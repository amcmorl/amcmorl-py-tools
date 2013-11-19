import numpy as np
from numpy import cos, sin #, array, dot
#from numpy import arctan2, sqrt
from . import norm, tensor_product, cross_matrix, unitvec

def angle_between(a, b):
    '''returns the angle (in rads) between 2 vectors'''

    costheta = np.dot( a, b ) / (norm( a ) * norm( b ))
    theta = np.arccos( costheta )
    return theta

def pt_nearest(pt, offs, dir):
    '''returns point ''new'' on line in direction ''dir'' through ''offs''
    nearest to point ''pt'' '''
    return offs + np.cross(dir, pt - offs) * dir

def axis_angle2mat(axis, angle):
    '''
    Construct rotation matrix from axis of rotation and angle.
    
    Parameters
    ----------
    axis : array_like
      vector of axis of rotation
    angle : float
      amount to rotate in radians
    '''
    axis = unitvec(axis)
    if np.rank(axis) > 1:
        raise ValueError('axis should be 1-d only')
    nd = axis.shape[0]
    xm = cross_matrix(axis)
    tp = tensor_product(axis, axis)
    c, s = cos(angle), sin(angle)
    I = np.identity(nd)
    return I * c + s * xm + (1 - c) * tp
    
def rotmat_between_two_vecs(u, v):
    '''
    Calculate the rotation matrix to transform from one vector `u`
    to another 'v'.
    
    Parameters
    ----------
    u, v : array_like
      1-d vectors
      
    Returns
    -------
    r : ndarray
      rotation matrix
    '''
    u = unitvec(u)
    v = unitvec(v)
    # if vectors are parallel, return Indentity
    if np.all(u == v):
        return np.identity(u.shape[0])
    axis = np.cross(u, v)
    angle = np.arccos(np.dot(u,v))
    return axis_angle2mat(axis, angle)

#~ def Rx(theta):
    #~ '''
    #~ Construct rotation matrix for rotation about x.
#~ 
    #~ Parameters
    #~ ----------
    #~ theta : float
      #~ angle in radians
#~ 
    #~ Returns
    #~ -------
    #~ rotation_matrix : ndarray
      #~ 3x3 rotation matrix
    #~ '''
    #~ t = theta
    #~ return array([[1,      0,       0],
                  #~ [0, cos(t), -sin(t)],
                  #~ [0, sin(t),  cos(t)]])
#~ 
#~ def Ry(theta):
    #~ '''
    #~ Construct rotation matrix for rotation about y.
#~ 
    #~ Parameters
    #~ ----------
    #~ theta : float
      #~ angle in radians
#~ 
    #~ Returns
    #~ -------
    #~ rotation_matrix : ndarray
      #~ 3x3 rotation matrix
    #~ '''
    #~ t = theta
    #~ return array([[ cos(t), 0, sin(t)],
                  #~ [ 0,      1,      0],
                  #~ [-sin(t), 0, cos(t)]])
#~ 
#~ def Rz(theta):
    #~ '''
    #~ Construct rotation matrix for rotation about z.
#~ 
    #~ Parameters
    #~ ----------
    #~ theta : float
      #~ angle in radians
#~ 
    #~ Returns
    #~ -------
    #~ rotation_matrix : ndarray
      #~ 3x3 rotation matrix
    #~ '''
    #~ t = theta
    #~ return array([[cos(t), -sin(t), 0],
                  #~ [sin(t),  cos(t), 0],
                  #~ [0,            0, 1]])

#~ def ypr2mat(ypr):
    #~ '''
    #~ Construct rotation matrix from yaw, pitch, roll.
#~ 
    #~ Parameters
    #~ ----------
    #~ ypr : array_like
      #~ shape (3,),  yaw, pitch and roll angles
#~ 
    #~ Returns
    #~ -------
    #~ rotation_matrix : ndarray
      #~ shape (3,3) rotation matrix
    #~ '''
    #~ return dot(dot(Rx(ypr[0]), Ry(ypr[1])), Rz(ypr[2]))

#~ def rm2ypr(rm):
    #~ '''
    #~ Convert a rotation matrix (rm) to roll-pitch-yaw representation.
    #~ '''
    #~ r21 = r[1,0]
    #~ r11 = r[0,0]
    #~ r31 = r[2,0]
    #~ r32 = r[2,1]
    #~ r33 = r[2,2]
    #~ 
    #~ alpha = arctan2(r21, r11)
    #~ beta  = arctan2(-r31, sqrt(r32**2 + r33**2))
    #~ gamma = arctan2(r32, r33)
#~ 
    #~ return alpha, beta, gamma
