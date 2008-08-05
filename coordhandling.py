from numpy import *
import scipy.interpolate
from vectors import unitvec
import numpy as n
import unittest

def identity3d( dims, iswhere=False ):
    '''creates a volume with 1s only along the 3D diagonal
    i.e. where x = y = z

    id3d = identity3d( dims )

    where dims is a tuple of the desired volume dimensions'''
    inds = n.indices( dims )
    arr = n.array((inds[0] == inds[1]) & (inds[0] == inds[2]))
    if not iswhere:
        return arr
    else:
        return n.where(arr)


def extract_box(data, centre, size):
    '''extracts a box from data, centered around "centre"
    (a tuple of co-ordinate) and with side of size "size"
    (a tuple of dimension sizes) or smaller if limited by
    edges of data

    usage: box = extract_box( data, centre, size )
    '''
    ndims = len(data.shape)
    try:
        assert ndims == len(centre)
        assert ndims == len(size)
    except AssertionError:
        print "data.shape, centre and size must all have same length"
        return None

    maxs = n.array( data.shape ) - 1
    mins = n.zeros( ndims )
    halfsize = n.array( size ) / 2
    diverr = n.array( size ) % 2
    lo = tuple( n.maximum( n.array( centre ) - halfsize, mins ) )
    hi = tuple( n.minimum( n.array( centre ) \
                              + halfsize + diverr, maxs ) )

    extr_str = ", ".join( [ 'lo[%d]:hi[%d]' % (i,i) for \
                            i in range( ndims ) ] )
    box = eval( 'data[' + extr_str + ']' )
    return box


def ninterpol(data, points, method='linear', mid=False):
    '''n-dimensional interpolation using sequential calls to
    scipy.interpolate.interp1d (I hope this is valid). Based
    loosely (signature only) on the routine of the same name
    in perldl.

    Usage:
    
        val = ninterpol( data, point )
        
    data = volume in which to interpolate
    point = co-ordinates of point to interpolate

    Probably has the same caveats as rebin.congrid (same engine)
    i.e. 1-D array must be promoted to shape (n,1).

    By definition the old values are defined at the lower bound of
    co-ordinate only (like points) rather than across the entire
    value (like pixels), which means that maximum offset _is_
    maxval - 1 so we don''t extrapolate
    '''
    ndims = len(data.shape)
    points = n.array( points )
    if mid:
        mid_ = 0.5
    else:
        mid_ = 0.

    # specify old dims
    olddims = [n.arange(i).astype(float) for i in list( data.shape )]
    
    # first interpolation - for ndims = any
    intrp = scipy.interpolate.interp1d( olddims[-1] + mid_, data, kind=method )
    newa = intrp( points[-1] )

    trorder = [ndims - 1] + range( ndims - 1 )
    for i in range( ndims - 2, -1, -1 ):
        newa = newa.transpose( trorder )
        intrp = scipy.interpolate.interp1d( olddims[i] + mid_, newa, \
                                            kind=method )
        newa = intrp( points[i] )
        
    if ndims > 1:
        # need one more transpose to return to original dimensions
        newa = newa.transpose( trorder )

    if ndims == 3:
        diag = identity3d( newa.shape, iswhere=True )
    elif ndims == 2:
        diag = n.where(identity( newa.shape[0] ))
    else:
        raise ValueError("ninterpol currently only supports 2 or 3 dimensions")
    return newa[diag]


def extract_line(data, ofs, vec, mid=True):
    '''extracts intensity values from a volume at points
    at unit length intervals along a line

    Usage: (vals, coords) = extract_line(data, ofs, vec)

    Inputs:
          data = volume to extract from
          ofs = point on line
          vec = direction of line
    Outputs:
          vals = returned values from volume
          coords = co-ordinates of extracted points'''

    if mid:
        mid_ = 0.5
    else:
        mid_ = 0.
    maxval = n.array( data.shape ) - 1

    # ensure inputs are numpy arrays
    ofs = n.asarray( ofs )
    vec = unitvec( vec )

    if n.alltrue( ofs <= maxval + mid_ ) and n.alltrue( ofs >= mid_ ):
        max_cnr = n.where(vec > 0, maxval, zeros(data.ndim)) + mid_
        min_cnr = n.where(vec < 0, maxval, zeros(data.ndim)) + mid_
        #print max_cnr
        #print min_cnr

        # work out how many steps before ofs
        presteps = (n.abs( min_cnr - ofs ) / n.abs( vec ) \
                    ).min().astype(int)

        # ... and how many after ofs
        poststeps = (n.abs( max_cnr - ofs ) / n.abs( vec ) \
                     ).min().astype(int)

        # construct list of steps ( in delta vecs )
        if presteps > 0:
            steps = [(presteps - i) * -1 \
                     for i in range( presteps + 1)]
                     # +1 to add 0 pt (at ofs)
        else:
            steps = [0]

        if poststeps > 0:
            steps += [(i + 1) for i in range( poststeps )]

        steps = n.array(steps)[newaxis,...]

        # construct array of actual pts
        pts = ofs[...,newaxis] + steps * vec[...,newaxis]
        #print pts
        val = ninterpol( data, pts, mid=mid )
        return val, pts
        
    else:
        raise ValueError("[extract_line] Offset must be within bounds of data.")
        return None


def ind2ax(ind, m):
    '''converts a 1-D index [ind] into an n-D
    index in an array of shape [m]'''
    if not type(m) == n.ndarray:
        m = n.array((m))
    mult = n.hstack((1, m[::-1].cumprod()[:-1] ))[::-1]
    cords = []
    for i in xrange(len(m)):
        this = ind / mult[i]
        cords.append(this)
        ind -= this * mult[i]
    return cords


def zmax(vol, pt):
    '''returns the index at which vol is maximum at the
    x-y co-ordinates determined by pt'''
    ind = n.where( vol[pt] == vol[pt].max().item() )
    #print ind
    if len(ind[0]) > 1:
        print "zmax Warning: More than one maximum found."
    return ind[0][0]


def mytest():
    print "Testing ind2ax...",
    a = zeros(1024)
    a[500] = 1
    a = a.reshape(8,16,8)
    if ind2ax(500, a.shape) == list(n.where(a == 1)):
        print 'passed'
    else:
        print 'failed'

    print "Testing zmax...",
    a = n.array([[[ 2, 10,  1, 10],
                [ 5,  2,  2,  1],
                [ 6,  8,  9,  3],
                [ 2, 10,  6,  6]],
               
               [[ 6,  2,  3,  2],
                [ 3,  8,  1, 10],
                [ 3,  9,  9,  1],
                [ 1,  9,  1,  1]],
               
               [[ 9,  1,  4,  0],
                [ 7,  8,  8,  8],
                [ 3,  3,  5,  7],
                [ 2,  2,  8,  5]],
               
               [[ 5,  8,  6, 10],
                [ 1,  6,  0,  7],
                [ 2,  9,  3,  8],
                [ 5,  1,  2,  9]]])
    if zmax(a, (0,2)) == 2:
        print 'passed'
    else:
        print 'failed'
        
    
def centre_of_mass(v):
    '''calculate centre-of-mass of volume

    \dot{c} = \dfrac{\sum{w \cdot \dot{p}_{i}}}
                    {\sum{w}}
    '''
    inds = n.indices(v.shape)
    vp = inds * v[n.newaxis,...]
    for i in range(n.rank(v)):
        vp = vp.sum(-1)
    return vp / v.sum()


def align_stacks_big_by_com(a, b):
    '''zero pad stacks to give compatible regions,
    aligned by centre of mass'''

    alc = centre_of_mass(a)
    blc = centre_of_mass(b)
    ga, gb = align_stacks_big_by_offs(a,alc,b,blc)
    
    return ga, gb


def align_stacks_big_by_offs(a, alc, b, blc):
    alc = n.asarray(alc)
    blc = n.asarray(blc)
    acu = n.asarray(a.shape) - alc
    bcu = n.asarray(b.shape) - blc

    glc = n.vstack((alc,blc)).max(0)
    gcu = n.vstack((acu,bcu)).max(0)

    ga = n.zeros(glc+gcu)
    gb = n.zeros(glc+gcu)
    al = glc - alc
    au = al + n.asarray(a.shape)
    ga_slice = [ slice(i,j) for (i,j) in zip( al, au ) ]
    bl = glc - blc
    bu = bl + n.asarray(b.shape)
    gb_slice = [ slice(i,j) for (i,j) in zip( bl, bu ) ]
    ga[ga_slice] = a
    gb[gb_slice] = b
    
    return ga, gb
    


class TestCoordhandlingFunctions(unittest.TestCase):

    def test_centre_of_mass(self):
        a = n.zeros((10,15)) + 1
        self.assertTrue( n.all(centre_of_mass(a) == n.array((4.5,7))) )

    def test_align_stacks_big_by_com(self):
        a = n.ones((11,15))
        # com_a = 5,7
        b = n.ones((15,21))
        b[:,0:3] = 9
        # com_b = 7,5.2
        ga, gb = align_stacks_big_by_com(a,b)
        res = n.zeros((15,22))
        # make like a
        res[2:13,:15] = 1
        # add b
        res[:,1:4] += 9
        res[:,4:] += 1
        self.assertTrue( n.all(res == ga+gb) )

    def test_align_stacks_big_by_offs(self):
        a = n.ones((11,15))
        ao = (5,7)
        b = n.ones((15,21))
        b[:,0:3] = 9
        bo = (7,5)
        ga, gb = align_stacks_big_by_com(a,b)
        res = n.zeros((15,22))
        # make like ga
        res[2:13,:15] = 1
        # add gb
        res[:,1:4] += 9
        res[:,4:] += 1
        self.assertTrue( n.all(res == ga+gb) )


def test():
    suite = unittest.TestLoader().loadTestsFromTestCase( \
        TestCoordhandlingFunctions)
    unittest.TextTestRunner(verbosity=2).run(suite)


import pylab as p

def graph_test_extract_line():
    v = (0.5,0.5,0.)
    a = (5,5,5)
    vol = n.indices((10,10,10))[0].astype(float)
    f = p.figure()
    ax1 = f.add_subplot(211)
    ax2 = f.add_subplot(212)
    vals, pts = extract_line(vol, a, v)
    ax1.imshow(vol[...,5].transpose(1,0))
    ax1.plot(pts[0], pts[1], 'ro')
    ax2.plot(vals, 'ro')
