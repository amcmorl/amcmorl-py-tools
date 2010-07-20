import numpy as np

def cart2pol(vector):
    '''Convert a vector into two angles

    Parameters
    ----------
    mu : array_like, shape (3,)
      x,y,z components of vector

    Returns
    -------
    theta : scalar
    phi : scalar
    '''
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

def pol2cart(theta, phi):
    '''Convert a point described by two angles to Cartesian co-ordinates.

    Parameters
    ----------
    theta : scalar
    phi : scalar
           
    Returns
    -------
    vector : array, shape (3,)
      x,y,z co-ordinates of equivalent vector
    '''
    sin, cos = np.sin, np.cos
    x = sin(theta) * cos(phi)
    y = sin(theta) * sin(phi)
    z = cos(theta)
    return np.array((x,y,z))
