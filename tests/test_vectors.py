from numpy.testing import *
import numpy as np
import vectors

def test_rotate_about_origin3d():
    before = np.array((1,0,0))
    normal = np.array((0,1,0))
    theta = np.pi/2. # clockwise rotation
    after = vectors.rotate_about_origin_3d(before, normal, theta)
    assert_almost_equal(after, np.array((0,0,-1)))
    
def test_rotate_by_angles():
    before = np.array((1.,0.,0.))
    theta = -np.pi/3. # upwards
    phi = np.pi/3.
    after = vectors.rotate_by_angles(before, theta, phi)
    assert_almost_equal(after, np.array((0.25, np.sqrt(3)/4., np.sqrt(3)/2.)))
    
    before = np.array((0., 1., 0.))
    theta = np.pi/3.
    phi = np.pi/4.
    after = vectors.rotate_by_angles(before, theta, phi)
    assert_almost_equal(after, np.array((-np.sqrt(2.)/2., np.sqrt(2)/2., 0.)))
