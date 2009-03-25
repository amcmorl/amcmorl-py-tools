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

def random_pars(mu=None, kappa=None, verbose=False):
    if mu == None:
        theta = np.random.uniform(low = -np.pi/2., high = np.pi/2.)
        phi = np.random.uniform(low = -np.pi, high = np.pi)
        if verbose: print "theta = %.3f, phi = %.3f" % (theta, phi)
        mu = convert_polar_to_cartesian(theta, phi)
    if kappa == None:
        kappa = np.random.uniform(low=0., high=5.)
    return mu, kappa

def assess_estimated_confid_angle(n_samples = int(1e3), n_pts=int(1e3),
				  verbose=False, alpha = 0.05,
				  kappa_range='any'):
    # generate a vMF distribution with a random mu (mean) and kappa (spread)
    x,y,z = 0,1,2
    mu, kappa = random_pars(verbose=verbose)
    if kappa_range == 'any':
        pass
    elif (kappa_range == '<5'):
        if (kappa > 5):
            kappa -= 5
    elif (kappa_range == '>5'):
        if (kappa < 5):
            kappa += 5
    else:
        try:
            kappa = float(kappa_range)
        except:
            print "Unknown kappa_range %s- ignored." % (kappa_range)
            
    count = 0
    for i_sample in xrange(n_samples):                
        # generate a random distribution with mu and kappa
        P_i = vmf_rvs(mu, kappa, n_pts=n_pts)
        
        # calculate confidence cone
        theta_confid = estimate_confid_angle(P_i, alpha, verbose=verbose)
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
        print "P = %5.3f" % (P)
        print "mu: %.3f %.3f %.3f, kappa %5.3f" % \
              (mu[x], mu[y], mu[z], kappa)
    return P

class TestVMFRVS(NumpyTestCase):
    mu, kappa = random_pars()
    n_pts = int(1e4)
    tolerance = 0.1
    #... acceptable angular deviation between set mean and sample mean
    P_i = vmf_rvs(mu, kappa, n_pts)
    sample_mean = mean_dir(P_i)
    angle_between = np.arccos(np.dot(mu, sample_mean))
    assert(angle_between < tolerance)
    
class TestEstimateKappa(NumpyTestCase):
    mu, kappa = random_pars()
    n_pts = int(1e4)
    tolerance = 0.1
    # acceptable deviation between set kappa (hopefully), and estimated kappa.
    # Could really do with a way to define a set of points with an established
    # kappa.
    P_i = vmf_rvs(mu, kappa, n_pts)
    khat = estimate_kappa(P_i, mu=mu)
    diff = np.abs(khat - kappa)
    assert(diff < tolerance)
    
# class TestEstimateConfidAngle(NumpyTestCase):
#     alpha = 0.05
#     acceptable_error = 0.05
#     n_samples = int(1e3)
#     if True:
#         # n = 100, kappa < 5
#         P = assess_estimated_confid_angle(n_samples=n_samples, n_pts=100,
#                             alpha=alpha, kappa_range='<5')
#         assert np.abs(P - (1 - alpha)) < acceptable_error
#         # n = 100, kappa > 5
#         P = assess_estimated_confid_angle(n_samples=n_samples, n_pts=100,
#                             alpha=alpha, kappa_range='>5')
#         assert np.abs(P - (1 - alpha)) < acceptable_error
#     if False:
#         P = assess_estimated_confid_angle(n_samples=n_samples, n_pts=50,
#                             alpha=alpha, kappa_range='<5')
#         assert np.abs(P - (1 - alpha)) < acceptable_error
#         # n = 50, kappa > 5
#         P = assess_estimated_confid_angle(n_samples=n_samples, n_pts=50,
#                             alpha=alpha, kappa_range='>5')
#         assert np.abs(P - (1 - alpha)) < acceptable_error
#     if True:
#         # n = 10, kappa < 5
#         P = assess_estimated_confid_angle(n_samples=n_samples, n_pts=10,
#                             alpha=alpha, kappa_range='<5')
#         assert np.abs(P - (1 - alpha)) < acceptable_error
#         # n = 10, kappa > 5
#         P = assess_estimated_confid_angle(n_samples=n_samples, n_pts=10,
#                             alpha=alpha, kappa_range='>5')
#         assert np.abs(P - (1 - alpha)) < acceptable_error
#     # exact kappa value
#     #P = assess_P_confid(n_samples=500, alpha=alpha, verbose=True,
#     #                    kappa_range='0.8084940961932108')
#     #assert np.abs(P - (1 - alpha)) < acceptable_error

def generate_Stephens1962_nomograms(resolution=3.):
    alphas = [0.1, 0.05, 0.01]
    n_alphas = len(alphas)
    Ns = range(3,17) + [18,20,25,30]
    n_Ns = len(Ns)
    coords = np.ma.empty((n_alphas, n_Ns, max(Ns) * resolution))
    coords.mask = ~np.isnan(coords)
    R0s = np.ma.empty((n_alphas, n_Ns, max(Ns) * resolution))
    R0s.mask = ~np.isnan(R0s)
    for i_alpha, alpha in enumerate(alphas):
        for i_N, N in enumerate(Ns):
            Rstars = np.linspace(0, N, N * resolution)
            for i_Rstar, Rstar in enumerate(Rstars):
                if Rstar == N:
                    Rstar -= 1e-1
                coords[i_alpha, i_N, i_Rstar] = Rstar
                R0 = calculate_R0(N, Rstar, alpha)
                R0s[i_alpha, i_N, i_Rstar] = R0
                R0s.mask[i_alpha, i_N, i_Rstar] = False
    return alphas, Ns, coords, R0s

def plot_Stephens1962_nomogram(alphas, Ns, Rstars, R0s,
                               figure=1, txt_dx=0.1, txt_dy=0):
    plt.subplots_adjust(top=0.95, hspace = 0.35)
    n_alphas, n_Ns, n_Rstars = Rstars.shape
    fig = plt.figure(num=figure)
    axes = []
    for i_alpha in xrange(n_alphas):
        # draw each plot - different alpha
        axes.append(fig.add_subplot(n_alphas, 1, i_alpha + 1))
        axes[-1].set_title('alpha = ' + str(alphas[i_alpha]))
        axes[-1].set_xlabel('R*')
        axes[-1].set_ylabel('R0')
        for i_N in xrange(n_Ns):
            # draw each line - different N
            axes[-1].plot(Rstars[i_alpha, i_N], R0s[i_alpha, i_N],
                          'k-o', markersize=2)
            Rstar_row = Rstars[i_alpha, i_N, ~Rstars.mask[i_alpha, i_N]]
            R0_row = R0s[i_alpha, i_N, ~Rstars.mask[i_alpha, i_N]]
            if Rstar_row.size != 0:
                last_Rstar = Rstar_row[-1]
                last_R0 = R0_row[-1]

            axes[-1].text(last_Rstar + txt_dx, last_R0 + txt_dy,
                          'N=' + str(Ns[i_N]), fontsize=7)
