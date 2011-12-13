from numpy import arctan2, sqrt

from warnings import warn
warn("This module is deprecated. Use vecgeom package instead.")

def rm2ypr(rm):
    '''
    Convert a rotation matrix (rm) to roll-pitch-yaw representation.
    '''
    r21 = r[1,0]
    r11 = r[0,0]
    r31 = r[2,0]
    r32 = r[2,1]
    r33 = r[2,2]
    
    alpha = arctan2(r21, r11)
    beta  = arctan2(-r31, sqrt(r32**2 + r33**2))
    gamma = arctan2(r32, r33)

    return alpha, beta, gamma
