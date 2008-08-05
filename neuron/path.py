# coding=UTF-8

import numpy as n
import unittest
import myutils
from myutils import pdbg
import pylab as p
import polygon
import os.path

myutils.debug = 'verbose'

nlocalspines = 4

# file stuff -----------------------------------------------------------------

testfile = '/home/amcmorl/lib/python-lib/neuron/test_cell.txt'


def read_points(fname):
    '''reads a points file and returns object array of lines'''
    f = open(fname, 'r')
    ids = []
    codes = []
    pos_s = []
    dis = []
    prnts = []
    line = None
    while not line == '':
        line = f.readline()
        if line <> '' and line[0] <> '#':
            bits = line.split()
            ids.append(int(bits[0]))
            codes.append(bits[1])
            pos_s.append(n.asarray([float(x) for x in bits[2:5]]))
            dis.append(float(bits[5]))
            prnts.append(int(bits[6]))

    ids = n.array(ids)
    codes = n.array(codes)
    pos_s = n.array(pos_s)
    dis = n.array(dis)
    prnts = n.array(prnts)
    return [codes, pos_s, dis, prnts]

def read_ldens(fname):
    lf = open(fname, 'r')
    lines = lf.readlines()
    lf.close()
    return n.array([float(x.strip()) for x in lines])

def read_resolution(fname):
    '''reads a resolution file in the same folder as fname'''
    cdir = os.path.split(fname)[0]
    f = open(cdir + '/info.txt')
    lines = f.readlines()
    x,y,z = [float(x) for x in lines[0].split()]
    f.close()
    return n.array((x,y,z))

def write_ldens(ldens):
    # should really check if it exists already
    f = open('ldens.txt', 'w')
    for i in xrange(ldens.shape[0]):
        f.write("%.4f\n" % ldens[i])
    f.close()

def sub_points(points, idx):
    '''select only a subtree of points - supports boolean indexing only'''
    most_points = [x[idx] for x in points]
    n_new_pts = most_points[0].shape[0]
    # if only it were this simple: return most_points

    # construct array of n_removeds
    parents = points[-1]
    n_pts = parents.shape[0]
    n_removeds = (1 - idx).cumsum()
    newparents = n.zeros(n_new_pts, dtype=int)
    j = -1
    for i in xrange(n_pts):
        if idx[i]:
            j += 1
            if parents[i] > 0:
                newparents[j] = parents[i] - n_removeds[parents[i]]
            elif parents[i] == -1:
                newparents[j] = -1
    most_points[-1] = newparents
    return most_points

def write_points(points, foutn):
    '''write points information to file   '''    
    # thread is in format x,y,z,di,code
    print "Writing file:", foutn
    fout = open(foutn, 'w')

    n_pts = points[0].shape[0]
    for i in xrange(n_pts):
        line = "%04d %s %7.3f %7.3f %7.3f %7.3f %4d\n" % \
               (i, points[0][i], \
                points[1][i,0], points[1][i,1], points[1][i,2], \
                points[2][i], points[3][i])
        fout.write(line)
    fout.close    

# basic path stuff -----------------------------------------------------------

def construct_children(points):
    # assumes each point has at most two children, and one of them is
    # directly below - doesn't allow for crossing of the soma
    parents = points[-1]
    npts = points[0].shape[0]
    children = n.zeros( npts, dtype=int )
    parofs = 0
    for i in xrange( npts ):
        #if i % 25 == 0: print i,
        #if parents[i] == -1: parofs = i
        if parents[i] <> 0: children[parents[i]] = i
    for i in xrange( npts ):
        #if i % 25 == 0: print i,
        if i == npts - 1: children[i] = -1
        # last point has no children
        elif parents[i+1] <> 0: children[i] = -1
    return children

hypot3d = lambda pt0, pt1: n.sqrt(((pt0 - pt1)**2).sum())

def cumulative_dist(points, children=None):
    if children == None:
            children = construct_children(points)
    npts = points[0].shape[0]
    parents = points[-1]
    pos_s = points[1]
    dists = n.zeros(npts, dtype=float)
    cdist = 0
    for i in xrange( npts ):
        if parents[i] == -1:
            dists[i] = 0
        else:
            if parents[i] == 0:
                par = i - 1
            elif parents[i] > 0:
                par = parents[i]
            delta = hypot3d(pos_s[par], pos_s[i])
            dists[i] = dists[par] + delta
    return dists

def reverse_dist(points, children=None, dists=None):
    '''calculate distance from terminals at each point'''
    if children == None:
        children = construct_children(points)
    if dists == None:
        dists = cumulative_dist(points, children)
    parents = points[-1]

    rdists = n.empty_like(dists)
    rdists[:] = 1e5 # shouldn't be any dendrites longer than this
    # need to parse properly through

    queue = []
    # queue terminal points
    tids = n.asarray(range(points[0].shape[0]))[children == -1]
    for tid in tids:
        queue.append( [tid, True, dists[tid]] ) # seed terminal points

    while queue:
        tpt, prox, refdist = queue.pop()
        tprnt = parents[tpt]
        #print "tpt: %3d" % tpt, 
        posdist = n.abs(refdist - dists[tpt])
        #print 'd = %6.2f' % posdist,
        # if this distance is smaller than currently in-place distance...
        if rdists[tpt] > posdist:
            rdists[tpt] = posdist

            # queue next point(s)
            if prox: # heading towards soma
                if tprnt == 0:
                    # next point is immediately above
                    queue.append( [tpt - 1, True, refdist] )
                    #print "add f", tpt -1
                    # need to check if parent has another child too
                    if children[tpt - 1] > 0:
                        newdist = dists[tpt - 1] - \
                                  (refdist - dists[tpt - 1])
                        queue.append( [ children[tpt - 1], False, newdist ])
                        
                elif tprnt > 0:
                    # queue parent
                    queue.append( [tprnt, True, refdist] )
                    #print "add f", tprnt
                    # queue parent's other child
                    #if parents[tprnt + 1] == 0:
                    newdist = dists[tprnt] - (refdist - dists[tprnt])
                    queue.append( [tprnt + 1, False, newdist] )
                    #print "add b", tprnt + 1
                else:
                    #print "soma"
                    continue # parent < 0 => reached soma, do next pt in queue
            else: # heading out
                if children[tpt] == -1:
                    #print "terminal"
                    continue
                else:
                    # queue any other children than might be present
                    if children[tpt] > 0:
                        #print "add", children[tpt]
                        queue.append( [children[tpt], prox, refdist] )
                    # queue next pt if its a child
                    if parents[tpt + 1] == 0:
                        #print "add", children[tpt]
                        queue.append( [tpt + 1, prox, refdist] )
        else:
            # we've already been closer to a terminal
            # go to next point in queue
            #print "too big"
            continue
    return rdists

def total_dist(points, children=None, dists=None):
    if children == None:
        children = construct_children(points)
    if dists == None:
        dists = cumulative_dist(points, children)
    parents = points[-1]
    totaldist = 0
    startdist = 0
    npts = points[0].shape[0]
    for i in xrange( npts ):
        if parents[i] > 0: startdist = dists[parents[i]]
        if parents[i] < 0: startdist = 0.
        # start dist = distance to this point
        if children[i] == -1: # end point
            thisend = dists[i]
            totaldist += (thisend - startdist)
    return totaldist

def dist_from_primary_mean(points, children=None, dists=None):
    if children == None:
        children = construct_children(points)
    if dists == None:
        dists = cumulative_dist(points, children)
    parents = points[-1]
    primaries = n.where(parents == -1)[0]
    #print "primary indexes:", primaries
    primeans = n.zeros(len(primaries))
    for i in xrange(len(primaries)):
        if i < len(primaries) - 1:
            tail = primaries[i+1]
        else:
            tail = None
        rng = dists[primaries[i]:tail]
        primeans[i] = rng.mean()
    #print "primary means:", primeans
    pridists = n.zeros_like(dists)
    for i in xrange(len(primaries)):
        if i < len(primaries) - 1:
            tail = primaries[i+1]
        else:
            tail = len(dists)
        for j in xrange(primaries[i],tail):
            pridists[j] = n.abs(dists[j] - primeans[i])
            #print "dist to soma %0.2f - mean dist %0.2f = %0.2f" % \
            #     (dists[j], primeans[i], pridists[j])
    return pridists

def deltas(points):
    npts = points[0].shape[0]
    parents = points[-1]
    pos_s = points[1]
    deltas = n.zeros( npts )
    for i in xrange( npts ):
        if parents[i] == -1:
            deltas[i] = 0
        elif parents[i] == 0:
            deltas[i] = hypot3d(pos_s[i - 1], pos_s[i])
        else:
            deltas[i] = hypot3d(pos_s[parents[i]], pos_s[i])
    return deltas

def soma_location(points):
    parents = points[-1]
    pos = points[1]
    soma_pts = n.where(parents == -1)
    return pos[soma_pts].mean(0)

# spine stuff ----------------------------------------------------------------

def average_density(points, totaldist=None, children=None, dists=None):
    '''returns average number of spines per micron across
    entire dataset'''
    if totaldist == None:
        if children == None:
            children = construct_children(points)
        if dists == None:
            dists = cumulative_dist(points, children)
        totaldist = total_dist(points, children, dists)
    nsps = (points[0] == 'p').sum() + \
          (points[0] == 'f').sum()
    #print "Number of spines", nsps
    spdens = nsps / totaldist
    return spdens

def local_density_distance(points, children=None, dists=None, count_dist=10.):
    if children == None:
        children = construct_children(points)
    if dists == None:
        dists = cumulative_dist(points, children)

    npts = points[0].shape[0]
    localdens = n.zeros( npts )
    parents = points[-1]
    # for each point
    for i in xrange( npts ):
        pdbg('terse', "Current point: %d" % i)
        # walk queues in each direction up to _dist_ distance
        # and count number of spines as you go
        # spine density = count / dist
        tcdist = dists[i] # this distance
        sp_count = 0

        # deal with case where current point is a spine
        if (points[0][i] == 'p') or (points[0][i] == 'f'):
            sp_count = 1
            pdbg('verbose', '#', '*')

        # set up queue
        queue = []
        # queue front and back
        if parents[i] <> -1:
            queue.append( [i, True, tcdist] ) # seed back direction
            pdbg('verbose',i-1, '*')
        if children[i] <> -1:
            queue.append( [i, False, tcdist] ) # seed forward direction
            pdbg('verbose',i, '*')
        pdbg('verbose','')

        while queue:
            tpt, prox, refdist = queue.pop()
            pdbg('verbose',"tpt:", tpt, '*')
            pdbg('verbose', 'distance = %0.2f' % \
                 n.abs( refdist - dists[tpt] ), '*')
            if n.abs( refdist - dists[tpt] ) < count_dist/2.:
                # if still close enough - test if is a spine
                if ((points[0][tpt] == 'p') or (points[0][tpt] == 'f')) \
                       and (i <> tpt):
                    # (don't count spines at central point
                    # because it's already been done
                    sp_count += 1
                    pdbg('verbose', '#','*')
                if prox:
                    # queue next points in same direction
                    pdbg('verbose', "-", '*')
                    # if this pt has a labelled parent...
                    if parents[tpt] > 0:
                        # queue it & it's other child
                        tprnt = parents[tpt]
                        queue.append( [tprnt, prox, refdist] )
                        pdbg('verbose', "add", tprnt, '*')
                        # queue parent's other child & change direction
                        if parents[tprnt + 1] == 0:
                            newdist = dists[tprnt] - (refdist - dists[tprnt])
                            queue.append( [tprnt + 1, not prox, newdist] )
                            pdbg('verbose',"add", tprnt + 1, '*')
                    elif parents[tpt] == -1:
                        # reached beginning
                        pass
                    else:
                        # parent is previous pt
                        queue.append( [tpt - 1, prox, refdist] )
                        pdbg('verbose',"add", tpt - 1, '*')
                        if (parents[tpt] == 0) and (children[tpt - 1] > 0):
                            # if parent has a labelled child, queue that
                            newdist = dists[tpt - 1] - \
                                      (refdist - dists[tpt - 1])
                            queue.append( [children[tpt - 1], \
                                           not prox, newdist] )
                            pdbg('verbose',"add", children[tpt - 1], \
                                 "refdist=", newdist, '*')
                else: # if not prox (==distal)
                    pdbg('verbose',"+", '*')
                    # walk forward one pt

                    if children[tpt] == -1:
                        # reached end of dendrite
                        pass
                    else:
                        if children[tpt] > 0:
                            pdbg('verbose', 'add', children[tpt], '*')
                            queue.append( [children[tpt], prox, refdist] )
                        if parents[tpt + 1] == 0:
                            # next pt is child
                            pdbg('verbose', 'add', tpt + 1, '*')
                            queue.append( [tpt+1, prox, refdist] )
                pdbg('verbose', '')
            else:
                # dist too big - stop
                pdbg('verbose', 'stopping: too big')
            
        pdbg('terse',"found %d spines" % sp_count)
        localdens[i] = sp_count / count_dist
        pdbg('terse', "dens here:", localdens[i], "\n")
    return localdens


def local_density_n_spines(points, children=None, dists=None):
    if children == None:
        children = construct_children(points)
    if dists == None:
        dists = cumulative_dist(points, children)

    npts = points[0].shape[0]
    localdens = n.zeros( npts )
    parents = points[-1]
    for i in xrange( npts ):
        pdbg('terse', "Current point: %d" % i)
        
        # count backwards four spines, collecting total_dist
        tcdist = dists[i] # this distance
        ntocount = nlocalspines
        locsps = n.ones( ntocount + 1 ) * 1e5
        # +1 to length gives room for one addition in [-1] slot
        # saves continual resizing
        
        # cdists of 8 spines flanking this pt
        # start with very high values and replace with lower ones

        pdbg('verbose', "i =", i, "seeds = ", '*')

        queue = [] # points/directions to check
        if parents[i] <> -1:
            queue.append( [i, -1, tcdist] ) # seed back direction
            pdbg('verbose',i-1, '*')
        if children[i] <> -1:
            queue.append( [i, 1, tcdist] ) # seed forward direction
            pdbg('verbose',i, '*')
        pdbg('verbose','')
        
        # while points in queue
        while queue:
            # check distance
            # check first point in queue
            tpt, way, refdist  = queue.pop()
            pdbg('verbose',"tpt:", tpt, '*')
            # walk to next pt - maybe
            if (n.abs( refdist - dists[tpt] ) < n.array(locsps).max()):
                # otherwise cdist is too big, can leave path

                # check if this point is a spine
                if (points[0][tpt] == 'p') or (points[0][tpt] == 'f'):
                    # is a spine or filopodium
                    tptpos = dists[tpt]
                    tptdist = n.abs(refdist - tptpos)
                    if tptdist == 0 and locsps[0] == 0.:
                        # don't count same self point in beth directions
                        pass
                    else:
                        locsps[-1] = tptdist
                        locsps.sort()
                    
                if way < 0:
                    pdbg('verbose', "-", '*')
                    # walk back one point
                    
                    #print '*', parents[tpt], children[tpt - 1]
                    if parents[tpt] > 0:
                        # if this pt has a labelled parent
                        # go to it & it's other child
                        tprnt = parents[tpt]
                        queue.append( [tprnt, way, refdist] )
                        pdbg('verbose', "add", tprnt, '*')
                        # queue parent's other child & change direction
                        if parents[tprnt + 1] == 0:
                            newdist = dists[tprnt] - (refdist - dists[tprnt])
                            queue.append( [tprnt + 1, way * -1, newdist] )
                            pdbg('verbose',"add", tprnt + 1, '*')
                    elif parents[tpt] == -1:
                        # reached beginning
                        pass
                    else:
                        # parent is previous pt
                        queue.append( [tpt - 1, way, refdist] )
                        pdbg('verbose',"add", tpt - 1, '*')
                        if (parents[tpt] == 0) and (children[tpt - 1] > 0):
                            #print refdist
                            #print dists[tpt - 1]
                            # if parent has a labelled child, queue that
                            newdist = dists[tpt - 1] - \
                                      (refdist - dists[tpt - 1])
                            queue.append( [children[tpt - 1], \
                                           way * -1, newdist] )
                            pdbg('verbose',"add", children[tpt - 1], \
                                  "refdist=", newdist, '*')

                else:
                    pdbg('verbose',"+", '*')
                    # walk forward one pt

                    if children[tpt] == -1:
                        # reached end of dendrite
                        pass
                    else:
                        if children[tpt] > 0:
                            pdbg('verbose', 'add', children[tpt], '*')
                            queue.append( [children[tpt], 1, refdist] )
                        if parents[tpt + 1] == 0:
                            # next pt is child
                            pdbg('verbose', 'add', tpt + 1, '*')
                            queue.append( [tpt+1, way, refdist] )
                pdbg('verbose', '')
            else:
                # dist too big - stop
                pdbg('verbose', 'stopping: too big')
            
        realdists = locsps[:-1][locsps[:-1] < 1e5]
        # -1 here adjusts for +1 in definition of locsps
        realdists.sort()
        pdbg('terse',"found dists", realdists)
        nrealdists = len(realdists)
        minidx = n.minimum(ntocount, nrealdists)
        pdbg('terse', "nsps found:", minidx)
        if minidx > 0:
            closestdists = realdists[0:minidx]
            localdens[i] = minidx / closestdists.sum()
        else:
            localdens[i] = 0
        pdbg('terse', "dens here:", localdens[i], "\n")
    return localdens

# branch stuff ---------------------------------------------------------------

def branch_order(points, children=None):
    if children == None:
        children = construct_children(points)
    npts = points[0].shape[0]
    orders = n.zeros( npts, dtype=n.uint8 )
    order = 1
    parents = points[-1]
    for i in xrange( npts ):
        if parents[i] == -1:
            order = 1
        if orders[i] <> 0:
            order = orders[i] # pick up this branch's order
        else:
            orders[i] = order # assign current order to this pt
        if children[i] > 0:
            order += 1 # step up an order
            orders[children[i]] = order
            # assign incremented order to other child
    return orders

   
def extents(points, children=None):
    '''calculates distance between soma point for each branch
    and each terminal point

    returns tuple:
      0 = xyz array of extents
      1 = sqrt(x^2 + y^2) = ''in-plane'' extent alone
    '''
    if children == None:
        children = construct_children(points)
    npts = points[0].shape[0]
    parents = points[-1]
    pos_s = points[1]
    soma = soma_location(points)
    endpts = n.where(children == -1)[0]
    n_endpts = endpts.shape[0]
    all_exts = (pos_s[endpts] - soma[n.newaxis,:]).transpose()
    all_xy_exts = all_exts[0:2]    
    z_exts = n.array((all_exts[2].max(), n.abs(all_exts[2].min())))
    xy_exts = n.sqrt((all_exts[0:2]**2).sum(0))
    xyz_exts = n.sqrt((all_exts**2).sum(0))
    return xyz_exts, xy_exts, z_exts


def branch_lengths(points, children=None, dists=None):
    if children == None:
        children = construct_children(points)
    if dists == None:
        dists = cumulative_dist(points, children)
    
    npts = points[0].shape[0]
    localdens = n.zeros( npts )
    parents = points[-1]
    branchlengths = []

    for i in xrange( npts ):
        # identify start dist
        if parents[i] == -1:
            # soma
            startdist = 0
        elif parents[i] > 0:
            # branches after a break
            startdist = dists[parents[i]]
        elif children[i-1] > 0:
            # follow on branches
            startdist = dists[i-1]
        # identify end pts
        if children[i] <> 0:
            branchlengths.append(dists[i] - startdist)

    return branchlengths

# wrappers -------------------------------------------------------------------

def cell_dend_stats(points, children=None, dists=None, orders=None, \
               tdist=None, text=False):
    if children == None:
        children = construct_children(points)
    if dists == None:
        dists = cumulative_dist(points, children)
    if orders == None:
        orders = branch_order(points, children)
    if tdist == None:
        tdist = total_dist(points, children, dists)

    data = {}

    # 1 maximum extent
    tot_exts, tv_exts, ax_exts = extents(points)
    data['max_extent_total'] = tot_exts.max()
    data['max_extent_transv'] = tv_exts.max()
    data['max_extent_axial'] = ax_exts.max()
    if text:
        #print "Maximum total extent: %.3f µm" % data['max_extent_total']
        print "Maximum in-plane extent: %.3f µm" % data['max_extent_inplane']
        print "Maximum axial extent: %.3f µm" % data['max_extent_axial']
    
    # 2 number of primary dendrites (parents == -1)
    parents = points[-1]
    data['num_primes'] = (parents == -1).sum()
    if text:
        print "Number of primary dendrites:", \
              data['num_primes']

    # 3 number of unbranched dendrites
    # ??

    # 4 number of terminals (children == -1)
    data['num_terms'] = (children == -1).sum()
    if text:
        print "Number of terminals:", data['num_terms']

    # 5 number of branch pts (children > 0)
    data['num_bpts'] = (children > 0).sum()
    if text:
        print "Number of branch pts:", data['num_bpts']

    data['num_branches'] = (children <> 0).sum()
    if text:
        print "Number of branches:", data['num_branches']

    # 6 combined dendritic length
    data['max_len'] = dists.max()
    data['tot_len'] = tdist
    if text:
        print "Maximum dendrite length: %.3f" % data['max_len']
        print "Total dendrite length: %.3f" % data['tot_len']
        
    return data

def branch_dend_stats(points, children=None, dists=None, orders=None, \
               tdist=None, text=False):
    ''' these ones need to be averaged across cells, so raw data
    should be passed (i.e. individual values) so that means can be
    calculated later
    '''
    if children == None:
        children = construct_children(points)
    if dists == None:
        dists = cumulative_dist(points, children)
    if orders == None:
        orders = branch_order(points, children)
    if tdist == None:
        tdist = total_dist(points, children, dists)

    data = {}

    # mean distance to branch points
    data['dist_bpts'] = dists[children > 0]
    #data['mean_dist_bpts'] = dists[children > 0].mean()
    if text:
        print "Mean distance to branch points: %.3f µm" % \
              data['dist_bpts'].mean()

    # mean distance to terminals
    data['dist_terms'] = dists[children == -1]
    if text:
        print "Mean distance to terminals: %.3f µm" % \
              data['dist_terms'].mean()

    # branch_lengths
    branchlens = branch_lengths(points)
    data['branch_lengths'] = branchlens
    if text:
        max_b_len = n.array(branchlens).max()
        mean_b_len = n.array(branchlens).mean()
        print "Maximum branch length: %.3f µm" % \
              max_b_len
        print "Mean branch length: %.3f µm" % \
              mean_b_len

    # average branch order
    # each branch sampled only once by children <> 0
    # (i.e. gets branch pts and terminals)
    data['b_order'] = orders[children <> 0]
    if text:
        print 'Mean branch order:   %.3f' % \
              data['b_order'].mean()

    # average terminal order
    end_orders = orders[children == -1]
    data['term_order'] = end_orders
    if text:
        print 'Mean terminal order: %.3f' % \
              data['term_order'].mean()

    return data

def orderwise_branch_stats(points):
    # number of branches at each branch order
    branch_orders = orders[children <> 0]
    branch_ord_fs, branch_ords = \
                   n.histogram(branch_orders, range(1, orders.max() + 1))
    if text:
        print 'Order    ' + ''.join(['%4d' % x for x in branch_ords])
        print 'Branches ' + ''.join(['%4d' % x for x in branch_ord_fs])

    # number of terminals at each branch order
    end_orders = orders[children == -1]
    end_ord_fs, end_ords = n.histogram(end_orders, range(1, orders.max() + 1))
    if text:
        print 'Terminals' + ''.join(['%4d' % x for x in end_ord_fs])

def spine_stats(points, children=None, dists=None, orders=None, \
               tdist=None, ldens=None, text=False):
    if children == None:
        children = construct_children(points)
    if dists == None:
        dists = cumulative_dist(points, children)
    if orders == None:
        orders = branch_order(points, children)
    if tdist == None:
        tdist = total_dist(points, children, dists)
    if ldens == None:
        ldens = local_density(points, children, dists)

    data = { }

    data['n_spines'] = (points[0] == 'p').sum() + \
                       (points[0] == 'f').sum()
    data['total_dist'] = total_dist(points, children, dists)
    #data['average density'] = average_density(points, tdist)
    if text:
        print "Mean spine density: %.3f" % (data['n_spines'] / \
                                            data['total_dist'])

    data['ldens'] = ldens
    if text:
        print "Maximum local spine density: %.3f" % ldens.max()
        print "Mean local spine density: %.3f" % ldens.mean()
        print "Minimum local spine density: %.3f" % ldens.min()

    spinepoints = n.logical_or(points[0] == 'p', points[0] == 'f')
    data['spords'] = orders[spinepoints]
    if text:
        print "Mean spine order: %.3f" % data['spords'].mean()
    data['spdists'] = dists[spinepoints]
    if text:
        print "Mean spine distance:  %.3f" % data['spdists'].mean()
    
    return data

def spine_distrib(points, children=None, dists=None, orders=None, \
                  text=False ):
    if children == None:
        children = construct_children(points)
    if dists == None:
        dists = cumulative_dist(points, children)
    if orders == None:
        orders = branch_order(points, children)

    data = {}

    # histogram of frequency vs distance from soma for spines
    data['dists'] = dists
    spinepoints = n.logical_or(points[0] == 'p', points[0] == 'f')
    data['spdists'] = dists[spinepoints]

    maxdist = dists.max()
    binsize = 25
    bins = n.arange(0, maxdist, binsize, dtype=float)

    # histogram of frequency vs order for spines
    data['orders'] = orders
    data['sp_orders'] = orders[spinepoints]
    
    if text:
        data['dhist'] = n.histogram(dists, bins=bins)
        data['sphist'] = n.histogram(spdists, bins=bins)
        data['normhist'] = sphist[0] / dhist[0].astype( float )
        
        f = p.figure(1)
        a0 = f.add_subplot(311)
        a1 = f.add_subplot(312)
        a2 = f.add_subplot(313)

        a0.bar(bins, sphist[0], width=0.9*bins[1])
        a1.bar(bins, dhist[0], width=0.9*bins[1])
        a2.bar(bins, normhist, width=0.9*bins[1])

    return data

def branch_density_profiles(points, children=None, dists=None, \
                                 ldens=None):
    # parse tree backwards, starting at terminals - build up profiles
    if children == None:
        children = construct_children(points)
    if dists == None:
        dists = cumulative_dist(points, children)
    if ldens == None:
        ldens = local_density_n_spines(points, children, dists)
    terms = n.where(children == -1)[0]
    parents = points[-1]
    densprofs = []
    distprofs = []
    for startpt in list(terms):
        tpt = startpt
        densprof = []
        distprof = []
        while tpt <> None:
            # start at point and walk backwards
            densprof.append(ldens[tpt])
            distprof.append(dists[tpt])
            if parents[tpt] > 0:
                tpt = parents[tpt]
            elif parents[tpt] == -1:
                tpt = None
            else:
                tpt -= 1
        densprofs.append(densprof)
        distprofs.append(distprof)
    return densprofs, distprofs
    
# subpoints stuff ------------------------------------------------------------

def find_sub_size(points, tlength):
    '''finds a size limit within which the total length of the remaining
points is close to tlength'''
    this_limit = points[1].max()
    this_length = total_dist(points)
    while this_length > tlength:
        this_limit -= 5 # increment in 5 micron steps
        subpoints = sub_points(points, n.all((points[1] < this_limit), axis=1))
        this_length = total_dist(subpoints)
    return this_limit, this_length, subpoints

# test stuff -----------------------------------------------------------------

class TestNeuronPathFunctions(unittest.TestCase):

    def setUp(self):
        testf = '/home/amcmorl/lib/python-lib/neuron/test_cell.txt'
        self.points = read_file(testf)
        self.children = construct_children(self.points)
        self.dists = cumulative_dist(self.points, self.children)
        self.tdist = total_dist(self.points, self.children, self.dists)
        self.deltas = deltas(self.points)
        

    def test_children(self):
        kids = construct_children(self.points)
        answer = n.array([ 0., 0., 0., 0., 0., 9., 0., 0., -1., \
                           0., 0., 0., 0., 0., -1., \
                           0., 0., -1])
        self.assertTrue( n.all(kids == answer) )


    def test_cumulative_dists(self):
        cdists = cumulative_dist(self.points, self.children)
        answer = n.array([ 0., 0.5, 1., 1.5, 2., 2.67082039, 2.87082039,  \
                           3.09442719, 3.19442719, 3.34164079, 4.06275104, \
                           4.28635784, 4.4277792, 4.62877671, 4.75683919, \
                           0., 0.5, 1.0])
        self.assertTrue( n.allclose( cdists, answer ) )


    def test_total_dist(self):
        tdist = total_dist(self.points, self.children, self.dists)
        answer = 6.2804459905
        self.assertTrue( n.allclose( tdist, answer ) )

    
    def test_average_density(self):
        avdens = average_density(self.points, self.tdist)
        answer = 5./6.2804459905
        self.assertTrue( n.allclose( avdens, answer ) )


    def test_local_density(self):
        avdens = local_density(self.points, self.children, self.dists)
        ds = self.deltas
        a0 = 4/n.array([ 1.5,  3.09442719,  4.28635784,  4.62877671]).sum()
        a4 = 4/n.array((ds[3:4].sum(), \
                        ds[5:8].sum(), \
                        ds[5:6].sum() + ds[9:12].sum(), \
                        ds[5:6].sum() + ds[9:14].sum())).sum()
        # need to also test a point with a spine and on second primary 
        self.assertTrue( n.allclose( avdens[0], a0 ) )


    def test_deltas(self):
        ds = deltas(self.points)
        answer = n.array([ 0., 0.5, 0.5, 0.5, 0.5, \
        0.67082039, 0.2, 0.2236068, 0.1, 0.67082039, \
        0.72111026, 0.2236068, 0.14142136, 0.20099751, 0.12806248, \
        0., 0.5, 0.5])
        self.assertTrue( n.allclose( ds, answer ) )


    def test_branch_lengths(self):
        bls = branch_lengths(self.points)
        ds = cumulative_dist(self.points)
        ans = [ ds[5] - ds[0],  \
                         ds[8] - ds[5],  \
                         ds[14] - ds[5], \
                         ds[17] - ds[15] ]
        self.assertTrue( n.allclose( bls, ans ) )

    
def test():
    suite = unittest.TestLoader().loadTestsFromTestCase( \
        TestNeuronPathFunctions)
    unittest.TextTestRunner(verbosity=2).run(suite)


