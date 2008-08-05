import numpy as n
from numpy.linalg import norm
from coordhandling import unitvec
import time

def cylinder( vdim, A, D, r, inds=None):
    '''usage: cyl = cylinder( vdim, A, D, r )
    creates a volume of size vdim containing a binary (1s inside)
    cylinder running through A, in direction D, with radius r'''
    
    D = unitvec( n.asarray( D, dtype=float ) )
    A = n.asarray( A, dtype=float )
    if inds == None:
        Xi = n.indices( vdim )
    else:
        Xi = inds
    Y = Xi - A[:, n.newaxis, n.newaxis, n.newaxis]   # vectors from A to Xi
    Ytr = Y.transpose(1,2,3,0)
    di = n.inner( D, Ytr )        # length along D to point nearest Xi
    dims = di.shape
    tmp = Y - di[n.newaxis,...] * D[:, n.newaxis, n.newaxis, n.newaxis]
    comp = n.sqrt( (tmp**2).sum(0) )
    vol = r >= comp
    return vol

def clip_plane(vdim, p, norm, inds=None):
    ''' Creates a volume containing 0s and 1s separated by a defined
    plane. Useful for clipping volumes (by multiplying by result).

    Usage: res = clip_plane( norm, p, vdim )

           where vdim = dimensions of volume
                 norm = normal to plane
                 p    = point on plane
                 

    Equation of plane (in vector notation):

           norm.(r-r0) = 0

    which becomes:

           a(x-x0) + b(y-y0) + c(z-z0) = 0

    when norm = (a, b, c) & r = (x, y, z)'''

    norm = n.asarray(norm)
    p = n.asarray(p)

    if inds == None:
        co = n.indices( vdim )
    else:
        co = inds
    pl = norm[...,n.newaxis,n.newaxis,n.newaxis] * \
         (co - p[...,n.newaxis,n.newaxis,n.newaxis] )
    pl = pl.sum(0)
    return (pl >= 0)

def sphere( vdim, c, r, inds=None ):
    '''create_sphere - return a volume containing a binary
    sphere defined by a point and a radius

Usage:

 vol = create_sphere( vdims, c, r );

where vdim = dimensions of volume
      c     = centre point
      r     = radius'''
    c = n.asarray(c)

    if inds == None:
        Xi = n.indices( vdim )
    else:
        Xi = inds

    return r >= n.sqrt( ((Xi - c[...,n.newaxis,n.newaxis,n.newaxis\
                                 ])**2).sum(0) )

def identity3d( dims, iswhere=False, inds=None ):
    '''Creates a volume with 1s only along the 3D diagonal
    i.e. where x = y = z

    id3d = identity3d( dims )

    where dims is a tuple of the desired volume dimensions.
    Also supports 2d.'''
    if inds == None:
        inds = n.indices( dims )
        
    ndims = len( dims )
    arr = n.array( (inds[0] == inds[1]) )
    if ndims == 3:
        arr = arr & n.array( (inds[0] == inds[2]) )
    if not iswhere:
        return arr
    else:
        return n.where(arr)

def circlecoords( cx, cy, r, start=0, finish=2*n.pi, npts=100 ):
    a = n.linspace(start, finish, npts)
    xs = cx + r * n.sin( a )
    ys = cy + r * n.cos( a )
    return xs, ys

def centre_pad(arr, shape, pad_value=0):
    big = n.zeros(shape) + pad_value
    bigc = (n.array(big.shape)/2.)
    lc = bigc - (n.array(arr.shape)/2.)
    uc = bigc + (n.array(arr.shape)/2.)
    slices = [ slice(start,stop) for (start,stop) in zip( lc.astype(int), \
                                                          uc.astype(int) ) ]
    big[slices] = arr
    return big

def draw_circle(radius, size=None, centre=None):
    # not quite perfect - some slight shift from centre
    if (centre == None) and (size != None):
        raise ValueError("If centre is specified, size must be also.")
    radius = n.asarray(radius)
    if size == None:
        size = n.asarray([radius * 2.] * 2)
    else:
        size = n.asarray(size)
    if centre == None:
        centre = size / 2.
    else:
        centre = n.asarray(centre)
    #print "centre", centre
    #print n.indices(size).shape
    inds = n.indices(size) - centre[...,n.newaxis,n.newaxis]
    d = n.sqrt((inds**2).sum(0))
    #print d
    return d < radius
