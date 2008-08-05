#from numpy import *
#from tidbits import *
import numpy as n
import scipy.interpolate
from scipy import ndimage

def rebin_factor(a, scale_factor):
    '''wraps rebin_neighbour to allow a scale factor to be given'''
    newshape = (n.array(a.shape, dtype=float) * scale_factor).round()
    return rebin_neighbour( a, newshape )


def neighbour( a, newshape ):
    '''Rebin an array to a new shape using nearest_neighbour lookup.
    '''
    assert len(a.shape) == len(newshape)
    
    slices = [ slice(0,old, float(old)/new) \
               for old,new in zip(a.shape,newshape) ]
    coordinates = n.mgrid[slices]
    indices = coordinates.astype('i')
    #choose the biggest smaller integer index
    return a[tuple(indices)]


def neighbour_factor( a, newshape ):
    '''Rebin an array to a new shape.
    newshape must be a factor of a.shape
    Uses nearest neighbour lookup.
    '''
    assert len(a.shape) == len(newshape)
    assert not n.sometrue(n.mod( a.shape, newshape ))
    
    slices = [ slice(None,None, old/new) \
               for old,new in zip(a.shape,newshape) ]
    return a[slices]


def rebin_average(a, *args):
    '''rebin ndarray data into a smaller ndarray of the same rank
    whose dimensions are factors of the original dimensions.
    eg. An array with 6 columns and 4 rows can be reduced to
    have 6,3,2 or 1 columns and 4,2 or 1 rows.
    example usages:
    >>> a=rand(6,4); b=rebin(a,3,2)
    >>> a=rand(6); b=rebin(a,2)
    Returns averages of replaced elements.
    '''
    shape = a.shape
    lenShape = len(shape)
    assert len(args) == lenShape
    evList = ['a.reshape('] + \
             ['args[%d],shape[%d]/args[%d],'%(i,i,i) \
              for i in range(lenShape)] \
             + [')'] + ['.sum(%d)'%(i+1) for i in range(lenShape)] + \
             ['*args[%d]/shape[%d]'%(i,i) for i in range(lenShape)]
    return eval(''.join(evList))


def rebin_mean(a, *args):
    '''returns a float array with values taken from the mean of the original
    pixels'''
    shape = a.shape
    lenShape = len(shape)
    factor = n.asarray(shape)/n.asarray(args)
    evList = ['a.reshape('] + \
             ['args[%d],factor[%d],'%(i,i) for i in range(lenShape)] + \
             [')'] + ['.mean(%d)'%(i+1) for i in range(lenShape)]
    print ''.join(evList)
    return eval(''.join(evList))


def congrid(a, newdims, method='linear', centre=False, minusone=False):
    '''Arbitrary resampling of source array to new dimension sizes.
    Currently only supports maintaining the same number of dimensions.
    To use 1-D arrays, first promote them to shape (x,1).
    
    Uses the same parameters and creates the same co-ordinate lookup points
    as IDL''s congrid routine, which apparently originally came from a VAX/VMS
    routine of the same name.

    method:
    neighbour - closest value from original data
    nearest and linear - uses n x 1-D interpolations using
                         scipy.interpolate.interp1d
    (see Numerical Recipes for validity of use of n 1-D interpolations)
    spline - uses ndimage.map_coordinates

    centre:
    True - interpolation points are at the centres of the bins
    False - points are at the front edge of the bin

    minusone:
    For example- inarray.shape = (i,j) & new dimensions = (x,y)
    False - inarray is resampled by factors of (i/x) * (j/y)
    True - inarray is resampled by(i-1)/(x-1) * (j-1)/(y-1)
    This prevents extrapolation one element beyond bounds of input array.
    '''
    if not a.dtype in [n.float64, n.float32]:
        a = n.cast[float](a)
    
    m1 = n.cast[int](minusone)
    ofs = n.cast[int](centre) * 0.5
    old = n.array( a.shape )
    ndims = len( a.shape )
    if len( newdims ) != ndims:
        print "[congrid] dimensions error. " \
              "This routine currently only support " \
              "rebinning to the same number of dimensions."
        return None
    newdims = n.asarray( newdims, dtype=float )    
    dimlist = []

    if method == 'neighbour':
        for i in range( ndims ):
            base = n.indices(newdims)[i]
            dimlist.append( (old[i] - m1) / (newdims[i] - m1) \
                            * (base + ofs) - ofs )
        cd = n.array( dimlist ).round().astype(int)
        newa = a[list( cd )]
        return newa
    
    elif method in ['nearest','linear']:
        # calculate new dims
        for i in range( ndims ):
            base = n.arange( newdims[i] )
            dimlist.append( (old[i] - m1) / (newdims[i] - m1) \
                            * (base + ofs) - ofs )
        # specify old dims
        olddims = [n.arange(i, dtype = n.float) for i in list( a.shape )]

        # first interpolation - for ndims = any
        mint = scipy.interpolate.interp1d( olddims[-1], a, kind=method )
        newa = mint( dimlist[-1] )

        trorder = [ndims - 1] + range( ndims - 1 )
        for i in range( ndims - 2, -1, -1 ):
            newa = newa.transpose( trorder )

            mint = scipy.interpolate.interp1d( olddims[i], newa, kind=method )
            newa = mint( dimlist[i] )

        if ndims > 1:
            # need one more transpose to return to original dimensions
            newa = newa.transpose( trorder )

        return newa
    elif method in ['spline']:
        oslices = [ slice(0,j) for j in old ]
        oldcoords = n.ogrid[oslices]
        nslices = [ slice(0,j) for j in list(newdims) ]
        newcoords = n.mgrid[nslices]

        newcoords_dims = range(n.rank(newcoords))
        #make first index last
        newcoords_dims.append(newcoords_dims.pop(0))
        newcoords_tr = newcoords.transpose(newcoords_dims)
        # makes a view that affects newcoords

        newcoords_tr += ofs        

        deltas = (n.asarray(old) - m1) / (newdims - m1)
        newcoords_tr *= deltas

        newcoords_tr -= ofs

        newa = ndimage.map_coordinates(a, newcoords)
        return newa
    else:
        print "Congrid error: Unrecognized interpolation type.\n", \
              "Currently only \'neighbour\', \'nearest\',\'linear\',", \
              "and \'spline\' are supported."
        return None



from gaussian import gauss3d
from scipy import signal

def rebin_arb(data, cur_px_res, des_px_res):
    '''blurs then resamples data to mimick
    rebinning over any arbitrary dimensions

    doesn''t really belong here, since it''s
    application is too specific'''

    orishape = data.shape
    cur_px_res = array(cur_px_res)
    des_px_res = array(des_px_res)

    # first make a gaussian kernel with fwhms
    # 2 * des_px_res at cur_px_res spacing
    # to allow for correct Nyquist sampling
    fwhms = des_px_res * 2 / cur_px_res
    print "[rebin_arb] Creating kernel..."
    kern = gauss3d( *fwhms )
    # tick

    # next blur data by that amount
    print "[rebin_arb] Convolving data..."
    blurred = signal.fftconvolve( data, kern )
    # tick
    
    # and cut original size out of convolved data
    # don't want to change volume size by resampling
    # (except directly)
    ofs = ((array(kern.shape) - 1)/2)
    indstr = ', '.join([("ofs[%d]:-ofs[%d] - 1" % (i, i)) \
                        for i in range( len( orishape ) )])
    cblur = eval("blurred[" + indstr + "]")
    # tick - okay see now I could have used 'same' switch to fftconvolve
    
    # next resample (using linear interpolation is now okay)
    # at desired points in data
    print "[rebin_arb] Resampling data..."
    new_dims = (array( orishape ) * cur_px_res / des_px_res).round().astype(int)
    smpl_cblur = congrid(cblur, new_dims)
    
    return (smpl_cblur, cblur, kern)
