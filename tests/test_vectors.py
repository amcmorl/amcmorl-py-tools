from numpy.testing import *
import numpy as np
import vectors

class TestRotateAboutOrigin3d(NumpyTestCase):
    before = np.array((1,0,0))
    normal = np.array((0,1,0))
    theta = np.pi/2. # clockwise rotation
    after = vectors.rotate_about_origin_3d(before, normal, theta)
    assert_almost_equal(after, np.array((0,0,-1)))
    
class TestRotateByAngles(NumpyTestCase):
    before = np.array((1,0,0))
    theta = np.pi/4.
    phi = np.pi/4.
    after = vectors.rotate_by_angles(before, theta, phi)
    assert_almost_equal(after, np.array((0.5, 0.5, -np.sqrt(2)/2.)))
    
class TestConvertAnglesToComponents(NumpyTestCase):
    theta = np.pi/3.
    phi = np.pi/3.
    x = 1/2. * np.sqrt(3/4.)
    y = 3/4.
    z = 1/2.
    res = vectors.convert_angles_to_components(theta, phi)
    assert_almost_equal(np.array((x,y,z)), res)
