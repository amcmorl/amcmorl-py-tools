from numpy.testing import *
import numpy as np
from spherical_stats import *

class TestConvertPolarToCartesian(NumpyTestCase):
    theta = np.pi/3.
    phi = np.pi/3.
    x = 1/2. * np.sqrt(3/4.)
    y = 3/4.
    z = 1/2.
    res = convert_polar_to_cartesian(theta, phi)
    assert_almost_equal(np.array((x,y,z)), res)

class TestConvertCartesianToPolar(NumpyTestCase):
    theta = np.pi/3.
    phi = np.pi/3.
    x = 1/2. * np.sqrt(3/4.)
    y = 3/4.
    z = 1/2.
    res = convert_cartesian_to_polar(np.array((x,y,z)))
    assert_almost_equal(np.array((theta, phi)), np.array(res))

class TestCoordinateConversions(NumpyTestCase):
    n_tests = int(1e3)
    for i in xrange(n_tests):
        theta = np.random.uniform(low = 0., high = np.pi)
        phi = np.random.uniform(low = 0., high = 2 * np.pi)
        cart = convert_polar_to_cartesian(theta, phi)
        theta_prime, phi_prime = convert_cartesian_to_polar(cart)
        assert_almost_equal(np.array((theta, phi)),
                            np.array((theta_prime, phi_prime)))
    
class TestEstimateConfidAngle(NumpyTestCase):
    n_samples = int(1e3)
    n_pts = 1e3
    verbose = False
    alpha = 0.05
    # generate a vMF distribution with a random mu (mean) and kappa (spread)

    mu, kappa = random_pars(verbose=verbose)
    x,y,z = 0,1,2
    count = 0
    for i_sample in xrange(n_samples):
#         if not verbose:
#             if i_sample % 10 == 0:
#                 print ""
#                 print "%02d" % (i_sample),
                
        # generate a random distribution with mu and kappa
        P_i = vmf_rvs(mu, kappa, n_pts=n_pts)
        
        # calculate confidence cone
        theta_confid = estimate_confid_angle(P_i, alpha)
        if verbose:
            print 'theta_confid=%5f' % (theta_confid),

        # see if population mean is less than confidence
        # cone angle away from sample mean
        R, S = calc_R(P_i)
        sample_mean = S/R
        inter_mean_angle = np.arccos(np.dot(mu, sample_mean))
        if verbose:
            print 'inter-mean angle=%5f' % (inter_mean_angle),
        if inter_mean_angle < theta_confid:
            count += 1
            if verbose:
                print " o"
        else:
            if verbose:
                print ""
    P = count / float(n_samples)
    if verbose:
        print ""
        print "P = %5.3f" % (P)
        print "mu: %.3f %.3f %.3f, kappa %5.3f" % \
              (mu[x], mu[y], mu[z], kappa)
    assert np.abs(P - (1 - alpha)) < 0.05
