import numpy as np

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

def rotate_by_angles(vector, theta, phi, reverse_order=False, fixlen=False):
    '''Gives vector after rotation about theta, phi.

    Parameters
    ----------
    vector : array_like, shape (3,)
      axial components of vector
    theta : scalar
      angle of rotation to z axis
    phi : scalar
      angle of rotation in x-y plane, CCW from x axis
    reverse_order : bool
      perform the phi rotation first? Normally, theta rotation is first.

    Returns
    -------
    rotated vector : array_like, shape (3,)
      axial components of rotated vector

    Notes
    -----
    Theta and phi are relative to (0,0,1).
    This performs two rotations:
      theta, about the y-axis; followed by phi, about the z-axis.
      '''
    t, ph = theta, phi
    A = np.array(([cos(t) * cos(ph), cos(t) * sin(ph), -sin(t)],
                  [-sin(ph),         cos(ph),          0],
                  [sin(t) * cos(ph), sin(t) * sin(ph), cos(t)]))
    if not reverse_order:
        A = A.T
    return np.dot(A, vector)
