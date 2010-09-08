import numpy as np

def edge2cen(array, axis=-1):
    '''
    Converts firing rates at edges of bins to firing rates in bin centers by
    averaging neighbouring bins.

    Parameters
    ----------
    array : array_like, shape (..., n)
        values at bin edges

    Returns
    -------
    array : array_like, shape (..., n - 1)
        values at bin centers
    '''
    ndim = np.rank(array)

    if ndim == 1:
        v = np.ones(2)
        res = np.convolve(array, v, mode='valid') / v.sum()
    else:
        if not ((axis == -1) | (axis == ndim - 1)):
            if axis < 0:
                axis = ndim - axis
            transpose_list = range(ndim)
            transpose_list[axis], transpose_list[-1] = -1, axis
            temp_array = np.transpose(array, transpose_list)
        else:
            temp_array = array

        shape = list(temp_array.shape)
        preshape = np.asarray(shape[0:-1])
        v = np.ones(2)
        vw = temp_array.reshape(np.product(preshape), shape[-1])
        temp_res = np.apply_along_axis(np.convolve, 1, vw,
                                       *[v, 'valid']) / v.sum()
        shape[-1] -= 1
        temp_res.shape = shape
        
        if not ((axis == -1) | (axis == ndim - 1)):
            res = np.transpose(temp_res, transpose_list)
        else:
            res = temp_res
    return res

def clip_below(x, llim=0):
    '''
    Clip array at lower bound (i.e. values lower than bound are set to bound)

    Parameters
    ----------
    x : ndarray
      array to clip, is left intact
    llim : scalar
      lower bound

    Returns
    -------
    clipped : ndarray
      `x` with clipped values
    '''
    x = x.copy()
    x[x < llim] = llim
    return x
