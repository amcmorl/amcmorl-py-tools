import numpy as np

def arg_ind_fn_where(fn, arr, cond):
    '''
    Return the index/indices into `arr` given by `fn`, but only for elements of `arr` specified by `cond`.

    Parameters
    ----------
    fn : callable
        argX (X =min, max, sort) function, takes a 1-d argument
    arr : array
        1-d array
    cond : sequence
        index to arr

    Returns
    -------
    ind : array
    '''
    assert np.rank(arr) == 1
    sub = fn(arr[cond])
    return np.arange(arr.size)[cond][sub]

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

def rolling_window(a, window):
   """
   Make an ndarray with a rolling window of the last dimension

   Parameters
   ----------
   a : array_like
       Array to add rolling window to
   window : int
       Size of rolling window

   Returns
   -------
   Array that is a view of the original array with a added dimension
   of size w.

   Examples
   --------
   >>> x = np.arange(10).reshape((2,5))
   >>> rolling_window(x, 3)
   array([[[0, 1, 2], [1, 2, 3], [2, 3, 4]],
          [[5, 6, 7], [6, 7, 8], [7, 8, 9]]])

   Calculate rolling mean of last dimension:
   >>> np.mean(rolling_window(x, 3), -1)
   array([[ 1.,  2.,  3.],
          [ 6.,  7.,  8.]])
          
   """
   if window < 1:
       raise ValueError, "`window` must be at least 1."
   if window > a.shape[-1]:
       raise ValueError, "`window` is too long."
   shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
   strides = a.strides + (a.strides[-1],)
   return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)
