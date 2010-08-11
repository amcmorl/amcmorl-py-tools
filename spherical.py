import numpy as np

def cart2pol(array, axis=-1):
    '''Convert an array containing xyz values into an array containing theta, phis

    Parameters
    ----------
    mu : array_like, shape (3,)
      x,y,z components of vector

    Returns
    -------
    theta : scalar
    phi : scalar
    '''
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
    return result.transpose(rlist)

def cart2pol_old(vector):
    x,y,z = vector

    theta = np.arccos(z)
    if np.allclose(theta, 0):
        # theta should be of size 1, but allclose handles masked values better
        # vector points straight up, phi is irrelevant,
        #... automatically define as 0
        phi = 0
    else:
        phi = np.arctan2(y,x)
    return theta, phi

def pol2cart(*args):
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
