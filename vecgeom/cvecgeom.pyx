import numpy as np
cimport numpy as np
from libc.math cimport sqrt

def norm3d(np.ndarray[np.float_t, ndim=1] arr):
    cdef double nrm = sqrt(arr[0]**2 + arr[1]**2 + arr[2]**2)
    return nrm
