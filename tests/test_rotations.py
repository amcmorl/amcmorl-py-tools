from numpy.testing import *
import numpy as np
from amcmorl_py_tools.vecgeom.rotations import \
    rotate_about_origin_3d, rotate_by_angles

def test_rotate_about_origin3d():
    before = np.array((1,0,0))
    normal = np.array((0,1,0))
    theta = np.pi/2. # clockwise rotation
    after = rotate_about_origin_3d(before, normal, theta)
    assert_almost_equal(after, np.array((0,0,-1)))
    
def test_rotate_by_angles():
    before = np.array([1.,0.,0.])
    theta = -np.pi/3. # upwards
    phi = np.pi/3.
    after = rotate_by_angles(before, theta, phi)
    assert_almost_equal(after, np.array((0.25, np.sqrt(3)/4., np.sqrt(3)/2.)))
    
    before = np.array([0., 1., 0.])
    theta = np.pi/3.
    phi = np.pi/4.
    after = rotate_by_angles(before, theta, phi)
    assert_almost_equal(after, np.array((-np.sqrt(2.)/2., np.sqrt(2)/2., 0.)))

    # test stacking of vectors
    before = np.array([[1., 0., 0.], [0., 0., 1.]]).T # shape (2, 3)
    theta = -np.pi/3.
    phi = np.pi/3.
    after = rotate_by_angles(before, theta, phi)
    assert_almost_equal(after, 
        np.array([[0.25, np.sqrt(3)/4., np.sqrt(3)/2.],
                  [np.sin(-np.pi/3.) * np.cos(np.pi/3.), -0.75, .5]]).T)
