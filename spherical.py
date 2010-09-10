import numpy as np

def cart2pol(array, axis=-1):
    '''Convert an array containing xyz values into an array containing theta, phis

    Parameters
    ----------
    array : array_like, shape (3,)
      x,y,z components of vector

    Returns
    -------
    theta : scalar
    phi : scalar
    '''
    array = np.asarray(array)
    if array.shape[axis] != 3:
        raise(ValueError("size of dimension %d must be 3"  % (axis)))

    return np.apply_along_axis(_cart2pol, axis, array)

def pol2cart(array, axis=-1):
    '''Convert a point described by two angles to Cartesian co-ordinates.

    Parameters
    ----------
    array : array_like
      theta, phi values (differentiated along `axis`, which must be length 2)
      
    Returns
    -------
    vector : array, shape (3,)
      x,y,z co-ordinates of equivalent vector
    '''
    array = np.asarray(array)
    if array.shape[axis] != 2:
        raise(ValueError("size of dimension %d must be 2" % (axis)))

    return np.apply_along_axis(_pol2cart, axis, array)

def _pol2cart(tp):
    theta, phi = tp
    sin, cos = np.sin, np.cos
    x = sin(theta) * cos(phi)
    y = sin(theta) * sin(phi)
    z = cos(theta)
    return np.array((x,y,z))

def _cart2pol(xyz):
    x,y,z = xyz
    theta = np.arccos(z)
    phi = np.arctan2(y, x)
    if np.allclose(phi, 0):
       phi = 0 
    return np.hstack((theta, phi))

def pol2cart_seq(*args):
    '''Convert a point described by two angles to Cartesian co-ordinates.

    Parameters
    ----------
    args : ndarray or sequence
      if ndarray, first dimension must be length 2: theta, phi
      if sequence, must have two elements: scalars or arrays of theta, phi
      
    Returns
    -------
    vector : array, shape (3,)
      x,y,z co-ordinates of equivalent vector
    '''
    if len(args) == 2:
        # assume are theta, phi
        theta, phi = args
    elif (len(args) == 1) & (type(args[0]) == np.ndarray):
        # assume rows are theta, phi
        theta = args[0][0]
        phi = args[0][1]

    sin, cos = np.sin, np.cos
    x = sin(theta) * cos(phi)
    y = sin(theta) * sin(phi)
    z = cos(theta)
    return np.array((x,y,z))
    
def cart2pol_v2(array, axis=-1):
    '''Convert an array containing xyz values into an array containing theta, phis

    Parameters
    ----------
    array : array_like, shape (3,)
      x,y,z components of vector

    Returns
    -------
    theta : scalar
    phi : scalar
    '''
    if np.rank(array) == 1:
        doflat = True
    else:
        doflat = False
    array = np.atleast_2d(array)

    if array.shape[axis] != 3:
        raise(ValueError("size of `axis` dimension must be 3"))
    
    # make axis first
    if axis < 0:
        axis = np.rank(array) + axis
    tlist = range(np.rank(array))
    tlist.pop(tlist.index(axis))
    tlist.insert(0, axis)
    
    x,y,z = array.transpose(tlist)
    theta = np.arccos(z)
    phi = np.arctan2(y,x)
    phi[np.abs(theta) < 1e-8] = 0

    # restore shape
    result = np.concatenate((theta[None],phi[None]))
    rlist = range(1, np.rank(result))
    rlist.insert(axis, 0)
    result = result.transpose(rlist)
    if doflat:
        return result.squeeze()
    else:
        return result
