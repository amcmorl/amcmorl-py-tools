import numpy as np
from numpy.linalg import norm

def perpz(vec):
    '''returns the unit length vector perpendicular
    to vec and lying in the x-y plane (i.e. z = 0).

    Works for 2 and 3-D vectors.
    '''
    try:
        assert ((len(vec) >= 2) & (len(vec) <= 3))
    except AssertionError:
        print "Perpz is only defined for 2 and 3-D vectors."

    if len(vec) == 3:
        return unitvec( (vec[1], -vec[0], 0) )
    else:
        return unitvec( (vec[1], -vec[0]) )

def unitvec(vec):
    '''returns the unit length vector
    in the same direction as vec'''
    return vec / np.sqrt(np.sum(vec ** 2))

def angle_between(a, b):
    '''returns the angle (in rads) between 2 vectors'''

    costheta = np.dot( a, b ) / (norm( a ) * norm( b ))
    theta = np.arccos( costheta )
    return theta

def pt_nearest(pt, offs, dir):
    '''returns point ''new'' on line in direction ''dir'' through ''offs''
    nearest to point ''pt'' '''
    return offs + np.cross(dir, pt - offs) * dir

def rotate_about_centre(v, c, th):
    '''returns vector v after rotation about centre c
    angle th is given in rads'''
    v = np.asarray(v)
    c = np.asarray(c)
    
    rotation_matrix = np.array([(np.cos(th), np.sin(th)), \
                               (-np.sin(th), np.cos(th))])
    return np.dot(rotation_matrix, v - c) + c

def rotate_about_origin_3d(vector, normal, theta):
    '''rotates the vector v around the normal vector n through angle th'''
    cos = np.cos
    sin = np.sin
    x, y, z = vector
    u, v, w = normal
    dt = u*x + v*y + w*z
    lns = u**2 + v**2 + w**2
    ln = np.sqrt(lns)
    return np.array(( \
        (u * dt \
         + (x * (v**2 + w**2) - u * (v*y + w*z)) * cos(theta) \
         + ln * (-w*y + v*z) * sin(theta)) / lns,                      
        (v * dt \
         + (y * (u**2 + w**2) - v * (u*x + w*z)) * cos(theta) \
         + ln * (w*x - u*z) * sin(theta)) / lns,
        (w * dt \
         + (z * (u**2 + v**2) - w * (u*x + v*y)) * cos(theta) \
         + ln * (-v*x + u*y) * sin(theta)) / lns \
        ))

def rotate_by_angles(vector, theta, phi):
    '''Gives vector after rotation about theta, phi.

    Parameters
    ----------
    vector : array_like, shape (3,)
      axial components of vector
    theta : scalar
      angle of rotation to z axis
    phi : scalar
      angle of rotation to x and y axes

    Returns
    -------
    rotated vector : array_like, shape (3,)
      axial components of rotated vector
      '''
    cos, sin = np.cos, np.sin
    t, ph = theta, phi
    A = np.array(([cos(t) * cos(ph), 0,       sin(t) * cos(ph)], \
                  [cos(t) * sin(ph), cos(ph), sin(t) * sin(ph)], \
                  [-sin(t),          0,       cos(t)]))
    return np.dot(A, vector)

def convert_angles_to_components(theta, phi):
    '''Convert a point described by two angles to Cartesian co-ordinates.

    Parameters
    ----------
    theta : scalar
      angle of rotation to z axis
    phi : scalar
      angle of rotation to x and y axes

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

def convert_components_to_angles(vector):
    '''Convert a vector into two angles

    Parameters
    ----------
    mu : array_like, shape (3,)
      x,y,z components of vector

    Returns
    -------
    theta : scalar
      angle from z axis
    phi : scalar
      angle in x/y axes
    '''
    x,y,z = vector
    pi = np.pi

    theta = np.arccos(z)
    phi = np.arctan(y/x) # first pass

    # now sanitize phi according to convention
    # (anticlockwise rotation, angle +ve)
    if (x > 0) & (y > 0):
        assert (phi > 0) & (phi < pi/2.)
    elif (x < 0) & (y > 0):
        phi = pi/2. - phi
        assert (phi > pi/2.) & (phi < pi)
    elif (x < 0) & (y < 0):
        phi += pi/2.
        assert (phi > pi) & (phi < 3 * pi / 2.)
    elif (x > 0) & (y < 0):
        phi = 2 * pi - phi
        assert (phi > 3 * pi / 2.) & (pi < 2 * pi)
    elif (x == 0) and (y < 0):
        phi = 3 * pi / 2.
    elif (x == 0) and (y > 0):
        phi = pi / 2.
    elif (x > 0) and (y == 0):
        phi = 0.
    elif (x < 0) and (y == 0):
        phi = pi
    return theta, phi
