from numpy import *
from numpy import linalg
from scipy.optimize import leastsq

def fit_plane_resids( p, pts ):
    '''returns the residual (distances to the plane for the
    points in pts) for a plane given by ax + by + cz + d = 0'''
    
    (a, b, c, d) = p
    Ds = (a * pts[:,0] + b* pts[:,1] + c * pts[:,2] + d) \
           / sqrt(a**2 + b**2 + c**2) 
    return Ds




def fit_plane( pts, p0=None ):
    '''fits a plane (in 3D) to an array of points (x,y,z co-ordinates)
    minimizing their orthogonal distance to the plane using leastsq'''
    if p0 == None:
        p0 = (1,1,1,1)
    plsq = leastsq( fit_plane_resids, p0, args=(pts,) )
    return plsq




def vol2coords(data, thr=0):
    '''converts a volume dataset into a list of
    co-ordinates and corresponding intensities
    at points where data > thr (default 0)'''
    
    cds = array( where(data > thr) )
    wi = data[tuple(cds)]
    return (cds, wi)




def dotover(a,b):
    '''Returns dot product of last two dimensions of two 2-D arrays,
    threaded over first dimension.'''
    try:
        assert a.shape[1] == b.shape[2]
        assert a.shape[0] == b.shape[0]
    except AssertionError:
        print "incorrect input shapes"
    res = zeros( (a.shape[0], a.shape[1], a.shape[1]), dtype=float )
    for i in range(a.shape[0]):
        res[i,...] = dot( a[i,...], b[i,...] )
    return res




def innerover(a,b):
    '''Returns inner product of last dimension of two 2-D arrays,
    threaded over first dimension.'''
    try:
        assert a.shape[0] == b.shape[0]
        assert a.shape[1] == b.shape[1]
    except AssertionError:
        print "incorrect input shapes"
    n = a.shape[0]
    res = zeros( (n, ), dtype=float )
    for i in range( n ):
        res[i] = inner( a[i,:], b[i,:] )
    return res




def fitlineNdw(x,w):
    '''fits a line (vector D through pt a)
    to co-ordinates x with intensities w'''

    N = x.shape[1]
    a = ((x * w**2).mean(1)) * N / (w**2).sum()
    y0 = w[:,newaxis] * (x.transpose() - a[newaxis,:])
    y = y0[:,newaxis,:]
    yt = y0[...,newaxis]
    M = -1 * dotover(yt, y).sum(0)
    M[arange( M.shape[0] ), arange( M.shape[1] )] += innerover( y0, y0 ).sum()
    (ev, e) = linalg.eig( M )
    D = e[:,where(ev == ev.min())].squeeze()
    return (a, D)
