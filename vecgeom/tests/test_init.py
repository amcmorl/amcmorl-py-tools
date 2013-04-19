import numpy as np
from amcmorl_py_tools.vecgeom import unitvec, unitvec_f2d

def test_unitvec():
    a = np.array([[ 0.50654606,  0.05050327],
       [ 0.06780228, -0.10952565],
       [ 0.12116112, -0.14544285],
       [-0.0588865 , -0.14017103],
       [ 0.1167503 , -0.26414753],
       [-0.09625524,  0.07777135],
       [ 0.32561687,  0.08549398],
       [-0.16084578, -0.0788045 ],
       [ 0.37862188, -0.05553404],
       [-0.06879143, -0.15628546]])
    uv = unitvec(a, axis=1)

    norm = np.sqrt(np.sum(a**2, axis=1))
    other = a / norm[...,None]
    np.testing.assert_array_equal(uv, other)

def test_unitvec_f2d():
    a = np.array([[ 0.50654606,  0.05050327],
       [ 0.06780228, -0.10952565],
       [ 0.12116112, -0.14544285],
       [-0.0588865 , -0.14017103],
       [ 0.1167503 , -0.26414753],
       [-0.09625524,  0.07777135],
       [ 0.32561687,  0.08549398],
       [-0.16084578, -0.0788045 ],
       [ 0.37862188, -0.05553404],
       [-0.06879143, -0.15628546]])
    uvf = unitvec_f2d(a)
    uv = unitvec(a, axis=1)
    np.testing.assert_array_equal(uv, uvf)
