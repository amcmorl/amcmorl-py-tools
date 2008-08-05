# coding=UTF-8

import neuron.path as np
import numpy as n
from glob import glob
from scipy import stats
from copy import copy
import os.path
import pylab as p
import biometry
import numpy.core.ma as ma
import fit1d

def add_dicts(a, *args):
    '''silently overwrites any values of a with the same key in nexts'''
    c = copy(a)
    for next in args:
        for k,v in next.iteritems():
            c[k] = v
    return c

class vals_accumulator:
    def __init__(self):
        self.data = {'n_spines' : 0,
                     'total_dist' : 0}

    def append(self, one_stats):
        if one_stats.has_key('n_spines'):
            self.data['n_spines'] += one_stats.pop('n_spines')
        if one_stats.has_key('total_dist'):
            self.data['total_dist'] += one_stats.pop('total_dist')
        for k,v in one_stats.iteritems():
            v = n.asarray(v)
            if self.data.has_key(k):
                # subsequent additions
                self.data[k] = n.hstack((self.data[k], v))
            else:
                # first addition
                self.data[k] = v

class hist_stats:
    def __init__(self):
        self.data = { }
        self.outs = { }

    def append(self, one_stats):
        for k,v in one_stats.iteritems():
            if self.data.has_key(k):
                # subsequent additions
                self.data[k] = n.hstack((self.data[k], n.asarray(v)))
            else:
                # first addition
                self.data[k] = n.array(v)
        
    def calc(self):
        if self.data.has_key('dists'):
            dbinsize = 25

            # calculate distance histogram
            self.outs['dbins'] = n.arange(0, self.data['dists'].max(), \
                                          dbinsize, dtype=float)
            self.outs['dhist'] = n.histogram(self.data['dists'], \
                                    bins=self.outs['dbins'])[0]
            self.outs['spdhist'] = n.histogram(self.data['spdists'], \
                                    bins=self.outs['dbins'])[0]
            self.outs['normdhist'] = self.outs['spdhist'] / \
                        self.outs['dhist'].astype(float)

            # calculate order histogram
            self.outs['obins'] = n.arange(1, self.data['orders'].max(), \
                                          dtype=float)
            self.outs['ohist'] = n.histogram(self.data['orders'], \
                                             bins=self.outs['obins'])[0]
            self.outs['spohist'] = n.histogram(self.data['sp_orders'], \
                                               bins=self.outs['obins'])[0]
            self.outs['normohist'] = self.outs['spohist'] / \
                                    self.outs['ohist'].astype(float)

def collect_group_stats(wdir):
    # get cells
    ages = ['P0-2', 'P9-11', 'P19-21'] #glob('*')
    lists = [ ]
    hists = [ ]
    
    for age in ages:
        print "Processing group %s..." % age 
        cells = glob(age + '/*/global_coords.txt')
        #cells = glob(age + '/*/shortened_coords.txt')
        vals = vals_accumulator()
        sp_dists = hist_stats()
        
        for cell in cells:
            print "Processing %s..." % cell
            pts = np.read_points(cell)
            vals.append( np.cell_dend_stats(pts) )
            vals.append( np.branch_dend_stats(pts) )
            ldfile = os.path.split(cell)[0] + '/ldens.txt'
            ldens = np.read_ldens(ldfile)
            vals.append( np.spine_stats(pts, ldens=ldens) )
            sp_dists.append( np.spine_distrib(pts) )
            
        sp_dists.calc()
        lists.append(vals)
        hists.append( sp_dists.outs )
        if vals.data.has_key('ldens'):
            maxld = vals.data['ldens'].max()
            binsize = 0.05
            nbins = n.divide(maxld - n.mod(maxld, binsize), binsize)
            bins = n.arange(nbins) * binsize
            hists[-1]['ldhists'] = n.histogram( vals.data['ldens'], bins=bins )
    return ages, lists, hists

def sigmark(p, alpha=0.05, ntests=3):
    if p < alpha/ntests:
        return "*"
    else:
        return " "

def print_report(wdir):
    ages, lists, hists = collect_group_stats(wdir)

    p.rcParams['figure.facecolor'] = 'w'
    p.rcParams['text.usetex'] = False
    
    fig_dists = p.figure(1)
    fig_ldens_all = p.figure(4)
    maxxd = 0.
    maxxld = 0.
    maxxo = 0.

    print_order = [ # table 1
        'max_extent_total', \
        'max_extent_transv', \
        'max_extent_axial', \
        'max_len', \
        'num_branches', \
        'num_primes', \
        'num_terms', \
        'tot_len', \
        '', \
        # table 2
        'b_order', \
        'branch_lengths', \
        'term_order', \
        'dist_terms', \
        '', \
        # table 3
        'spine_density', \
        'ldens', \
        'spords', \
        'spdists']

    linetype=['b-o', 'r-o', 'y-o']

    barcolor=['blue', 'red', 'yellow']
                   
    for gp, age in enumerate(ages):
        print "Age %s" % age
        for k in print_order:
            if k == 'spine_density':
                if lists[gp].data['total_dist'] > 0:
                    print "%18s : %8.3f (%d µm)" % \
                          (k, lists[gp].data['n_spines'] / \
                           lists[gp].data['total_dist'], \
                           lists[gp].data['total_dist'])
            elif k == '':
                print ''
            else:
                if lists[gp].data.has_key(k):
                    print "%18s : %8.3f ± %8.5f (%d); <%.2f" % \
                          (k, \
                           lists[gp].data[k].mean(),
                           stats.sem(lists[gp].data[k]), \
                           lists[gp].data[k].size, \
                           lists[gp].data[k].max())
        
        # local density histogram
        if hists[gp].has_key('ldhists'):

            # all
            bins = hists[gp]['ldhists'][1]
            ax_ldens_all = fig_ldens_all.add_subplot(111)
            maxxld = n.maximum(bins.max() + 1. * bins[1], maxxld)
            heights = hists[gp]['ldhists'][0] + 1e-20
            ax_ldens_all.bar(bins, heights+1e-12, color=barcolor[gp], \
                         width=1.*(bins[1] - bins[0]), log=True, alpha=0.5 )

        if hists[gp].has_key('dbins'):
            # distance histogram
            subplot_num = [len(ages),1,gp + 1]
            ax = fig_dists.add_subplot(*subplot_num)
            dbins = hists[gp]['dbins']
            maxxd = n.maximum(dbins.max() + 1. * dbins[1], maxxd)
            ax.bar(dbins, hists[gp]['normdhist'], \
                   width = 1. * dbins[1])
            ax.set_xlabel('distance (um)')
            ax.set_ylabel('frequency')
            ax.set_title('Age group %d' % gp) 

    for ax in fig_dists.axes:
        ax.set_xlim((0, maxxd))
    ax_ldens_all.set_xlabel('local densities (spines/um)')
    ax_ldens_all.set_ylim(ymin=1.)
    ax_ldens_all.set_ylabel('frequency')
    #ax_ldens_all.legend()


    # stats testing
    print ""
    print "%18s   %16s %14s %14s %14s" % \
          ('P(diff)',  'Effect', \
           'P0-2*P9-11', \
           'P0-2*P19-21', \
           'P9-11*P19-21')
    for k in print_order:
        if lists[0].data.has_key(k) and \
               lists[1].data.has_key(k) and \
               lists[2].data.has_key(k):
            H, pe = stats.kruskal( lists[0].data[k], \
                                lists[1].data[k], \
                                lists[2].data[k] )
        else:
            H, pe = n.nan, n.nan
        if lists[0].data.has_key(k) and lists[1].data.has_key(k):
            t01, p01 = biometry.wilcoxon(lists[0].data[k], lists[1].data[k])
        else:
            t01, p01 = n.nan, n.nan
        if lists[0].data.has_key(k) and lists[2].data.has_key(k):
            t02, p02 = biometry.wilcoxon(lists[0].data[k], lists[2].data[k])
        else:
            t02, p02 = n.nan, n.nan
        if lists[1].data.has_key(k) and lists[2].data.has_key(k):
            t12, p12 = biometry.wilcoxon(lists[1].data[k], lists[2].data[k])
        else:
            t12, p12 = n.nan, n.nan
        print "%18s : %7.2f;%8.2e%1s %5.2f;%8.2e%1s" \
        "%5.2f;%8.2e%1s %5.2f;%8.2e%1s" % \
              (k, \
               H, pe, sigmark(pe),  \
               t01, p01, sigmark(p01), \
               t02, p02, sigmark(p02), \
               t12, p12, sigmark(p12))
            
def plot_branch_density_profiles(wdir):
    ages = ['P0-2', 'P9-11'] #glob('*')

    colors = ['r', 'y', 'g', 'c', 'b', 'm']
    
    # collect data
    max_num_cells = 0
    min_num_terms = 1e5 # not likely to be beaten
    alldensprofs = []
    alldistprofs = []
    for age in ages:
        print "Processing group %s..." % age 
        cells = glob(age + '/*/global_coords.txt')
        agedensprofs = []
        agedistprofs = []
        for i, cell in enumerate(cells):
            print "Processing %s..." % cell, i + 1
            if (i + 1) > max_num_cells:
                max_num_cells = (i + 1)
            pts = np.read_points(cell)
            ldfile = os.path.split(cell)[0] + '/ldens.txt'
            ldens = np.read_ldens(ldfile)
            celldensprofs, celldistprofs = \
                       np.branch_density_profiles(pts, ldens=ldens)
            if len(celldensprofs) < min_num_terms:
                min_num_terms = len(celldensprofs)
            agedensprofs.append(celldensprofs)
            agedistprofs.append(celldistprofs)
        alldensprofs.append(agedensprofs)
        alldistprofs.append(agedistprofs)
            
    # plot data
    p.rcParams['figure.facecolor'] = 'w'
    p.rcParams['text.usetex'] = False
    f = p.figure()
    i = 0
    print "Min number of terminals = %d" % min_num_terms
    h = 0
    maxdens = 0
    for agedensprofs, agedistprofs in zip(alldensprofs, alldistprofs):
        print "Plotting age %s..." % ages[h]
        for celldensprofs, celldistprofs in zip(agedensprofs, agedistprofs):
            print "Ploting cell %d which has %d terminals..." % \
                  (i, len(celldensprofs))
            ax = f.add_subplot(len(ages), max_num_cells, i+1)
            i += 1
            ax.set_title(ages[h])
            k = 0
            for j in xrange(min_num_terms):
                idx = n.random.uniform(high=len(celldensprofs))
                distprof = celldistprofs.pop(idx)
                densprof = celldensprofs.pop(idx)
                if max(densprof) > maxdens:
                    maxdens = max(densprof)
                ax.plot(distprof, densprof, \
                        colors[n.mod(k,len(colors))])
                ax.plot((distprof[0], ), (densprof[0],),
                        'o' + colors[n.mod(k,len(colors))])
                k += 1
        h += 1

    for ax in f.axes:
        ax.set_ylim(ymax=maxdens)

def fit_densities(wdir):
    ages, lists, hists = collect_group_stats(wdir)
    fig = p.figure()
    ax = fig.add_subplot(111)
    maxxld = 0.
    plotcolor=['blue', 'red', 'yellow']    
    for gp, age in enumerate(ages):
        print "Age %s" % age
        # local density histogram
        if lists[gp].data.has_key('ldens'):
            obs = lists[gp].data['ldens']
            obs.sort()
            n_ = obs.size
            i = n.arange(n_)
            F0_5 = biometry.F_delta(0.5, n_, i)
            #F0 = biometry.F_delta(0, n_, i)
            #F1 = biometry.F_delta(1, n_, i)
            ax.plot(obs, F0_5, linestyle='steps', \
                    color=plotcolor[gp], label = age + \
                    ' observed')

            # fit to obs cdf
            p_obs = fit1d.fit_exp_cdf(F0_5, r=obs)
            F_hat = fit1d.exp_cdf(p_obs, r=obs)
            ax.plot(obs, F_hat, linestyle='--', \
                    color=plotcolor[gp], label = age+' fitted')
            ax.legend()

            # test significance
            biometry.ks_samples(obs, F_hat)
            #print "test D %.3f; p %.3f" % (D, pval)
            
    ax.set_ylim(ymin=0)
    ax.set_ylabel("Cumulative relative frequency")    
    return i, F0_5, obs, p_obs

def fit_density_histograms(wdir):
    eps = 1e-20
    ages, lists, hists = collect_group_stats(wdir)
    fig = p.figure()
    ax = fig.add_subplot(111)
    maxxld = 0.
    plotcolor=['blue', 'red', 'yellow']    
    for gp, age in enumerate(ages):
        print "Age %s" % age
        if lists[gp].data.has_key('ldens'):
            maxld = lists[gp].data['ldens'].max()
            mean = lists[gp].data['ldens'].mean()
            binsize = 0.05
            nbins = n.divide(maxld - n.mod(maxld, binsize), binsize)
            bins = n.arange(nbins) * binsize
            hist = n.histogram( lists[gp].data['ldens'], bins=bins )
            
            #bins = hists[gp]['ldhists'][1]
            maxxld = n.maximum(bins.max() + 1. * bins[1], maxxld)
            f_obs = hist[0] + eps
            ax.bar(bins, f_obs + 1e-12, color=plotcolor[gp], \
                   width=1.*(bins[1] - bins[0]), log=True, alpha=0.5 )
            class_means = bins + (bins[1] - bins[0])/2.
            f_obs_m = ma.masked_array(f_obs, mask=f_obs < 1e-10)
            ax.plot(class_means, f_obs_m, 'o' + plotcolor[gp][0])
            regres = stats.linregress(class_means, n.log(f_obs_m))
            invtau,logA = regres[0:2]
            f_pred = n.exp(logA) * n.exp(class_means * invtau )
            ax.plot(class_means, f_pred, '-' + plotcolor[gp][0])

def fit_cell_ldens(wdir):
    p.rcParams['figure.facecolor'] = 'w'
    p.rcParams['text.usetex'] = False
    ages = ['P0-2', 'P9-11'] #glob('*')
    fig = p.figure()
    ax_cdfs = fig.add_subplot(211)
    ax_cdfs.set_ylabel('relative cumulative frequency')
    ax_cdfs.set_xlabel('spine density (spines/um)')
    ax_diffs = fig.add_subplot(212)
    ax_diffs.set_ylabel('cdf_obs. - cdf_pred.')
    ax_diffs.set_xlabel('spine density (spines/um)')
    maxxld = 0.
    plotcolor=['blue', 'red', 'yellow']
    plotsymbol=['o', '^']
    for j, age in enumerate(ages):
        print "Processing group %s..." % age 
        cells = glob(age + '/*/ldens_nspines.txt')
        sp_dists = hist_stats()
        for k, f in enumerate(cells):
            print "Processing cell %s..." % age
            ldens = np.read_ldens(f)
            tau = ldens.mean()
            print "Tau = %.3f" % tau
            ldens.sort()
            n_ = ldens.size
            i = n.arange(n_)
            F0_5 = biometry.F_delta(0.5, n_, i)
            ax_cdfs.plot(ldens, F0_5, '-' + plotcolor[j][0], \
                    label=age+' observed')
            F_hat = (1 - n.exp(-ldens/tau))
            ax_cdfs.plot(ldens, F_hat, '--' + plotcolor[j][0], \
                    label=age+' predicted')
            ax_cdfs.legend()
            ax_cdfs.set_ylim(ymin = 0)
            biometry.ks_samples(ldens, F_hat)
            ax_diffs.plot(ldens, F0_5 - F_hat, '-' + plotcolor[j][0], \
                          label=age+' difference')

def get_compartment_pts(pts, soma_dist, mid_dist, term_dist):
    dists = np.cumulative_dist(pts)
    rdists = np.reverse_dist(pts)
    # identify pts closest to soma
    apts = dists < soma_dist
    # identify mean distance from soma
    pdists = np.dist_from_primary_mean(pts)
    # identify pts within x/2 from mean
    bpts = n.abs(pdists) < mid_dist
    # identify pts closest to terminals
    # need to parse through points working backwards from termini
    cpts = rdists < term_dist
    Apts = apts & ~cpts
    Bpts = bpts & ~apts & ~cpts
    Cpts = cpts & ~apts
    return Apts, Bpts, Cpts

def combine_compartment_codes(Apts, Bpts, Cpts, vals):
    all = n.vstack((Apts, Bpts, Cpts)).astype(int)
    weighted = all * vals.reshape(3,1)
    return weighted.max(axis=0)

def visualise_compartments(pts, Apts=None, Bpts=None, Cpts=None):
    from sources_scatter_source import ScatterSource
    src = ScatterSource()
    locs = pts[1]
    src.points = locs[1:]
    drs = (locs - n.roll(locs, 1, axis=0))[1:]
    src.vector_data = drs
    return src

def write_region_dens_file(fname, age, cell, soma, mid, term):
    if not os.path.exists(fname):
        f = open(fname, 'w')
        f.write('%4s %3s %4s %6s %7s' % \
                ('', 'Age', 'Cell', 'Region', 'Density\n'))
        offset = 0
    else:
        f = open(fname, 'r')
        offset = int(f.readlines()[-1].split()[0])
        f.close()
        f = open(fname, 'a')
    idx = 0
    for i, reg in enumerate([soma, mid, term]):
        if i == 0: region = 's'
        elif i == 1: region = 'm'
        else: region = 't'
        for j, dens in enumerate(reg):
            f.write("%4d %3s %4d %6s %7.4f\n" % \
                    (offset + idx, age, cell, region, dens))
            idx += 1
    f.close()

def extract_compartment(wdir, x, y, z):
    '''divides dendritic tree points into three categories:
    (A) x micron from soma
    (B) x/2. micron from mean distance
    (C) x micron from terminals
    where each category is exclusive of each other.'''
    # need to identify cells & regions independently
    # loop over cells
    ages = ['P0-2', 'P9-11'] #glob('*')
    for j, age in enumerate(ages):
        if age == 'P0-2': symbol = 'D'
        else: symbol = 'o'
        print "Processing group %s..." % age
        cells = glob(age + '/*/global_coords.txt')
        for i, cell in enumerate(cells):
            # load cell pts and densities
            ldfile = os.path.split(cell)[0] + '/ldens_nspines.txt'
            pts = np.read_points(cell)
            if os.path.exists(ldfile):
                ldens = np.read_ldens(ldfile)
            Apts, Bpts, Cpts = get_compartment_pts(pts, x, y, z)
            somadens = ldens[Apts]
            middens  = ldens[Bpts]
            termdens = ldens[Cpts]

            outfname = wdir + 'dens20.txt'
            write_region_dens_file(outfname, age[1], i, \
                                   somadens, middens, termdens)
            
            offset = i + j * len(cells)
            p.plot(n.ones_like(somadens) * (offset - 0.05), \
                   somadens,'c' + symbol)
            p.plot(n.ones_like(middens) * (offset), \
                   middens,'y' + symbol)
            p.plot(n.ones_like(termdens) * (offset + 0.05), \
                   termdens,'r' + symbol)

def load_compartment_dens(file):
    f = open(file, 'r')
    lines = f.readlines()
    dens = n.zeros(len(lines) - 1)
    ages = n.zeros(len(lines) - 1)
    regn = n.zeros(len(lines) - 1)
    for i, line in enumerate(lines[1:]):
        regn_code = line.split()[3]
        if regn_code == 's': regn_num = 0
        elif regn_code == 'm': regn_num = 1
        else: regn_num = 2
        regn[i] = regn_num
        ages[i] = int(line.split()[1])
        dens[i] = float(line.split()[-1])
    return regn, dens, ages

def compare_regions(regn, dens, ages):
    neonates = ages == 0
    juveniles = ages == 9
    #1/0.
    assert(n.all(neonates == ~juveniles))

    ssn = dens[(regn == 0) & neonates]
    msn = dens[(regn == 1) & neonates]
    tsn = dens[(regn == 2) & neonates]
    ssnm, ssne = ssn.mean(), stats.sem(ssn)
    msnm, msne = msn.mean(), stats.sem(msn)
    tsnm, tsne = tsn.mean(), stats.sem(tsn)

    ssj = dens[(regn == 0) & juveniles]
    msj = dens[(regn == 1) & juveniles]
    tsj = dens[(regn == 2) & juveniles]
    ssjm, ssje = ssj.mean(), stats.sem(ssj)
    msjm, msje = msj.mean(), stats.sem(msj)
    tsjm, tsje = tsj.mean(), stats.sem(tsj)
    
    ax = p.figure().add_subplot(111)
    ax.bar(n.arange(3), (ssnm, msnm, tsnm), width = (0.4), \
           color=('#00ffff', '#ffff00', 'r'), yerr=(ssne, msne, tsne))
    ax.bar(n.arange(3) + 0.4, (ssjm, msjm, tsjm), width = (0.4), hatch='/', 
           color=('#00ffff', '#ffff00', 'r'), yerr=(ssje, msje, tsje))
    ax.set_xticks(n.arange(3) + 0.4)
    ax.set_xticklabels(['Proximal', 'Mid', 'Distal'])
    ax.set_xlim(0, 2.8)
    print "** Neonates **"
    print "Ns: proximal = %d, mid = %d, distal = %d" % \
          (((regn == 0) & neonates).sum(), \
           ((regn == 1) & neonates).sum(), \
           ((regn == 2) & neonates).sum())
    print "Means: proximal = %5.3f±%5.3e, " \
          "mid = %5.3f±%5.3e, "\
          "distal = %5.3f±%5.3e" % \
          (ssnm, ssne, msnm, msne, tsnm, tsne)
    print "Kruskal", stats.kruskal(ssn, msn, tsn)
    print "Wilcoxon - proximal x mid", biometry.wilcoxon(ssn, msn)
    print "Wilcoxon - mid x distal", biometry.wilcoxon(msn, tsn)
    print "Wilcoxon - proximal x distal", biometry.wilcoxon(ssn, tsn)

    print "** Juveniles **"
    print "Ns: proximal = %d, mid = %d, distal = %d" % \
          (((regn == 0) & juveniles).sum(), \
           ((regn == 1) & juveniles).sum(), \
           ((regn == 2) & juveniles).sum())
    print "Means: proximal = %5.3f±%5.3e, " \
          "mid = %5.3f±%5.3e, "\
          "distal = %5.3f±%5.3e" % \
          (ssjm, ssje, msjm, msje, tsjm, tsje)
    print "Kruskal", stats.kruskal(ssj, msj, tsj)
    print "Wilcoxon - proximal x mid", biometry.wilcoxon(ssj, msj)
    print "Wilcoxon - mid x distal", biometry.wilcoxon(msj, tsj)
    print "Wilcoxon - proximal x distal", biometry.wilcoxon(ssj, tsj)
    return ax

def better_density_distance_plot_cell(pts, binsize=25, step=1, \
                                      color='red', axes=None):
    ''' plot spine density as a function of distance from the
    soma for each cell'''
    # find max distance from soma
    dists = np.cumulative_dist(pts)
    maxdist = n.ceil(dists.max())
    xdists = n.arange(binsize/2., maxdist - binsize/2., step)
    spines = (pts[0] == 'p') | (pts[0] == 'f')
    print "Number of spines:", spines.sum()
    ydens = n.zeros(xdists.size)
    # loop over distances from soma
    for i, xd in enumerate(xdists):
        minx = xd - binsize/2.
        maxx = xd + binsize/2.
        vdists = (dists > minx) & (dists < maxx)
        ydens[i] = (spines & vdists).sum().astype(float) \
                   / binsize / vdists.sum()
    if axes == None:
        axes = p.figure().add_subplot(111)
    axes.plot(xdists, ydens, color=color)
    # define circle as location +/- 1/2 binsize = 25 micron
    # calculate number of spines in that circle
    # present density as number of spines / binsize
    return axes
    
def plot_better_density_distance_plot():
    bdir = '/home/amcmorl/working/cells/'
    files = ['P0-2/040623/global_coords.txt', \
             'P0-2/040624/global_coords.txt', \
             'P9-11/050822/global_coords.txt', \
             'P9-11/050824/global_coords.txt']
    colors = ['b','b','r','r']
    ax = None
    for f,c in zip(files,colors):
        pts = np.read_points(bdir + f)
        ax = better_density_distance_plot_cell(pts, binsize=50, \
                                               color=c, axes=ax)
    
def intercell_spine_densities(wdir):
    ages = ['P0-2', 'P9-11', 'P19-21'] #glob('*')
    for j, age in enumerate(ages):
        print "Processing group %s..." % age
        cells = glob(age + '/*/global_coords.txt')
        meandens = []
        for cell in cells:
            pts = np.read_points(cell)
            meandens.append(np.average_density(pts))
        meandens = n.asarray(meandens)
        meanmean = meandens.mean()
        meansem = stats.sem(meandens)
        print "Age: %s; density = %5.3f ± %5.3e spines/µm" % \
              (age, meanmean, meansem)

