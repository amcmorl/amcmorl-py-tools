import numpy as np
from spherical import pol2cart, cart2pol
from vectors import norm

def test_cart2pol():
    v = np.array([0., 0., 1.])
    np.testing.assert_equal(cart2pol(v), np.array([0.,0.]))

    v = np.array([0., 0., -1.])
    np.testing.assert_equal(cart2pol(v), np.array([np.pi, 0.]))

    v = np.array([1., 0., 0.])
    np.testing.assert_equal(cart2pol(v), np.array([np.pi/2., 0.]))

    v = np.array([1., 1., 0.])
    v /= norm(v)[...,None]
    np.testing.assert_equal(cart2pol(v), np.array([np.pi/2., np.pi/4.]))
    
    v = np.array([0., 1., 0.])
    v /= norm(v)[...,None]
    np.testing.assert_equal(cart2pol(v), np.array([np.pi/2., np.pi/2.]))

    v = np.array([[0., 0., -1.],[1., 0., 0.]]) # shape 2,3
    exp = np.array([[np.pi, 0.], [np.pi/2., 0.]])
    np.testing.assert_equal(cart2pol(v), exp)

    v = np.array([[0., 0., -1.],[1., 0., 0.]]).T # shape 2,3
    exp = np.array([[np.pi, 0.], [np.pi/2., 0.]]).T
    np.testing.assert_equal(cart2pol(v, axis=0), exp)

    
def test_pol2cart_cart2pol():
    v = np.random.random(size=(2, 3))
    v /= norm(v, axis=1)[...,None]
    np.testing.assert_almost_equal(v, pol2cart(cart2pol(v)))
