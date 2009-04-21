import numpy as np
from enthought.mayavi import mlab
import scipy.optimize as opt
from scipy.integrate import quad
import vectors
import matplotlib.pyplot as plt
from scipy.io import read_array
import point_primitives

def convert_polar_to_cartesian(theta, phi):
    '''Convert a point described by two angles to Cartesian co-ordinates.

    Parameters
    ----------
    theta : scalar
      angle of rotation to z axis
    phi : scalar
      angle of rotation in x-y plane, counter-clockwise from x axis
      +ve y is defined to be next CCW from +x, as in:

            +y
            /\
            ||
      -x <--  --> +x
            ||
            \/
            -y
           
    Returns
    -------
    vector : array, shape (3,)
      x,y,z co-ordinates of equivalent vector
    '''
    sin, cos = np.sin, np.cos
    x = sin(theta) * cos(phi)
    y = sin(theta) * sin(phi)
    z = cos(theta)
    return np.array((x,y,z))

def convert_cartesian_to_polar(vector):
    '''Convert a vector into two angles

    Parameters
    ----------
    mu : array_like, shape (3,)
      x,y,z components of vector

    Returns
    -------
    theta : scalar
      angle from z axis
    phi : scalar
      angle in x-y plane, from counter-clockwise from x axis
      +ve y is defined to be next CCW from +x
    '''
    x,y,z = vector
    pi = np.pi

    theta = np.arccos(z)
    if theta == 0:
        # vector points straight up, phi is irrelevant,
        #... automatically define as 0
        phi = 0.
    else:
        phi = np.arctan(y/x) # first pass
        # now sanitize phi according to convention
        # (anticlockwise rotation, angle +ve)
        if (x > 0) & (y > 0):
            assert (phi > 0) & (phi < pi/2.)
        elif (x < 0) & (y > 0):
            phi += pi
            assert (phi > pi /2.) & (phi < pi)
        elif (x < 0) & (y < 0):
            phi += pi
            assert (phi > pi) & (phi < 3 * pi / 2.)
        elif (x > 0) & (y < 0):
            phi += 2 * pi
            assert (phi > 3 * pi / 2.) & (pi < 2 * pi)
        elif (x > 0) and (y == 0): # +x axis
            phi = 0.
        elif (x == 0) and (y > 0): # +y axis
            phi = pi / 2.
        elif (x < 0) and (y == 0): # -x axis
            phi = pi
        elif (x == 0) and (y < 0): # -y axis
            phi = 3 * pi / 2.
    return theta, phi

# -------------------------- von-Mises-Fisher dist ---------------------------

def calc_R(P_i):
    '''Returns R, the resultant vector, of a collection of vectors.'''
    S = np.sum(P_i, axis=0)   # 3.2
    R = np.sqrt(np.sum(S**2)) # 3.3
    return R, S

def mean_dir(P_i):
    '''Calculates the mean direction of a set of unit vectors on a sphere.

    Parameters
    ----------
    P_i : array_like, shape (n_pts, 3)
        n * x,y,z components of vectors

    Returns
    -------
    vector : array_like, shape (3,)
        x,y,z components of mean direction vector
    '''
    R, S = calc_R(P_i)
    return S/R

def estimate_kappa(P_i, mu=None):
    '''Returns the maximum likelihood estimate of the spread parameter (kappa)
    of a von-Mises-Fisher distribution.

    Parameters
    ----------
      pts : array_like, shape (n, 3)
        n unit vectors for which to find \kappa
      mu : array_like, shape (3,)
        vector of direction cosines for known population mean

    Returns
    -------
      kappa : scalar
        spread estimate for sample, assuming von-Mises-Fisher distribution
      
    Notes
    -----
    Code is based on equations in Statistical Analysis of Spherical Data,
    1987 by NI Fisher, T Lewis, and BJJ Embleton, section 5.3.2(iv).

    coth(\kappa) - 1/\kappa = R/n

    where R/n is the mean resultant vector.

    In the case where the true mean direction is known, an improvement on the
    above is to substitute R^* for R, where

    R^* = RC

    where C = \lambda\hat\lambda + \mu\hat\mu + v\hat{v}

    where \lambda, \mu, and v are the mean direction cosines and those with a
    hat are the sample statistic equivalents. Since this is the dot product
    of (\lambda, \mu, v) and (\hat\lamda, \hat\mu, \hatv), C is the cosine of
    the angle between the sample and population means.

    Numbers in the code comments refer to equation numbers in the text.'''
    
    n = P_i.shape[0]
    R, S = calc_R(P_i)
    #S = np.sum(P_i, axis=0)   # 3.2
    #R = np.sqrt(np.sum(S**2)) # 3.3
    if mu != None:
        sample_mean = S/R
        C = np.dot(mu, sample_mean)
        R *= C
    k0 = 1.
    lsq = opt.fmin_powell(to_min, k0, args=(R, n), disp=0)
    return lsq
    
def to_min(k, R, n):
    err = np.abs(1/np.tanh(k) - 1/k - R/n)
    #print 'k %.8f R %.4f n %d -> %.8f' % (k, R, n, err)
    return err

def C_F(kappa):
    '''From Statistical Analysis of Spherical Data,
    1987 by NI Fisher, T Lewis, and BJJ Embleton, equation 4.21'''
    return kappa / (2 * np.pi * (np.exp(kappa) - np.exp(-kappa)))

def vmf_pdf(theta, kappa):
    '''Returns the probability density function for the
    von-Mises-Fisher function, for the 3d case, at the angle theta.

    Equation is taken from Statistical Analysis of Spherical Data,
    1987 by NI Fisher, T Lewis, and BJJ Embleton, equation 4.23

    Parameters
    ----------
    theta : scalar
      angle for which to calculate f(theta)
    kappa : scalar
      spread parameter
    '''
    return C_F(kappa) * np.exp(kappa * np.cos(theta))

def vmf_cdf(theta, kappa):
    '''Returns the value of the cumulative density function of the
    von-Mises-Fisher function, 3d case, at angle theta.'''
    theta = np.asarray(theta)
    if theta.size > 1:
        if np.rank(theta) > 1:
            raise UserError("theta should be at most 1d.")
        cdf_list = [quad(vmf_pdf, 0, angle, args=(kappa))[0] for angle in theta]
        cdf = np.array(cdf_list)
        return cdf
    else:
        intgr, abserr = quad(vmf_pdf, 0, theta, args=(kappa))
        #print abserr
        return intgr

def vmf_rvs(mu, k, n_pts=1):
    '''Simulate a  Fisher distribution. From Statistical Analysis of
    Spherical Data, 1987, section 3.6.2'''
    size = (n_pts, 2)
    R = np.random.uniform(size=size)
    R_0 = R[...,0]
    R_1 = R[...,1]
    l = np.exp(-2 * k)
    colat = 2 * np.arcsin(np.sqrt(-np.log(R_0 * (1 - l) + l) / (2 * k)))
    longit = 2 * np.pi * R_1

    pts = convert_polar_to_cartesian(colat, longit)
    #print mean_dir(pts.T)
    alpha, beta = convert_cartesian_to_polar(mu)
    return vectors.rotate_by_angles(pts, alpha, beta).T

#---------------------------------------------------------------------------
# Stephens, 1962 Biometrika paper - derivation of table A.12 in Fisher, 1987
def factorial(n):
    assert (n >= 0)
    if n > 12:
        n = np.uint64(n)
    return np.product(np.arange(1, n+1))

def binomial_coef(n,k):
    return factorial(n) / (factorial(k) * factorial(n-k))

def lclip(a):
    return a if a >= 0 else 0

def P(N, R):
    '''from Stephens, 1962b, eqn 4'''
    ss = np.arange(0, N / 2.)
    rs = np.empty(ss.size)
    for i, s in enumerate(ss):
        rs[i] = (-1)**s * binomial_coef(N,s) * (lclip(N - R - 2 * s))**(N - 1)
    result = rs.sum()
    assert(~np.isnan(result))
    return rs.sum()

def calculate_R0(N, Rstar, alpha): #Ro_low, Ro_high):
    Ro_low = 0
    Ro_high = N
    return opt.fminbound(prob_fn, Ro_low, Ro_high,
                         args=(N, Rstar, alpha))#, disp=0)
    #return gradient_descent(prob_fn, start=3*Rstar/4., args=(N, Rstar, alpha))

def prob_fn(Ro, N, Rstar, alpha):
    return (P(N, Ro) / P(N, Rstar) - alpha)**2
#----------------------------------------------------------------------------

def estimate_confid_angle(P_i, alpha, mu=None, verbose=False):
    log = np.log
    sqrt = np.sqrt
    arccos = np.arccos
    z = 2
    n = P_i.shape[0]
    R, S = calc_R(P_i)
    if mu == None:
        mu = S/R
    #Rz = R * mu[z]
    Rz = R
    
    khat = estimate_kappa(P_i, mu=mu)
    #print "khat %3f" % (khat)
    if khat < 5:
        if n <= 30:
            R0 = calculate_R0(n, Rz, alpha)
            theta_alpha = arccos(R0 / R)
        elif (n > 8) & (n < 30):
            if (Rz > 0) & (Rz < n/4.):
                Ra1 = sqrt(Rz**2 - 2 / 3. * n * log( alpha ))
                theta_alpha = arccos(Ra1 / R)
            elif (Rz >= n/4.) and (Rz <= 3 * n / 5.):
                Ra1 = sqrt(Rz**2 - 2 / 3. * n * log( alpha ))
                astar = alpha**(1 / (n - 1))
                Ra2 = n * (1 - astar) + Rz * astar
                Ra = (Ra1 + Ra2) / 2.
                theta_alpha = arccos(Ra / R)
            else: # Rz > 3n/5
                assert (Rz > (3 * n / 5))
                astar = alpha**(1 / (n - 1))
                Ra2 = n * (1 - astar) + Rz * astar
                theta_alpha = arccos(Ra2 / R)
        else: # n > 30
            if verbose: print "-",
            assert(n >= 30)
            # !!!!!!!!!!!!!!!!!!!!!!!!
            theta_alpha = arccos(1 + log(alpha)/(khat * R))
    else: # khat >= 5
        if verbose: print "*",
        # !!!!!!!!!!!!!!!!!!!!!!!!            
        theta_alpha = arccos(1 - ((n - R) / R) *
                             ((1 / alpha)**(1/float(n-1)) - 1))
    return theta_alpha

def calc_confid_angle(P_i, alpha, mu=None):
    '''calculates the (1 - alpha) * 100% confidence interval for the given
    distribution of spherically distributed pts P_i

    Parameters
    ----------
    P_i : array, shape (n,3)
      co-ordinates of n pts (x,y,z) in distribution
    alpha : float
      probability of acceptable alpha error
    mu : array_like, shape (3,)
      mean direction of distribution, in event it is known a priori of
      variates

    Returns
    -------
    theta_alpha : float
      angle from the mean vector describing confidence cone'''
    z = 2

    n = P_i.shape[0]
    R, S = calc_R(P_i)
    if mu == None:
        mu = S/R
    Rz = R * mu[z]
    # n <= 8 find R^0_{1-\alpha}
    R0z = calculate_R0(n, Rz, alpha)
    print "Rz=%f, R0z=%f" % (Rz, R0z)
    theta_alpha = np.arccos(R0z / R)
    return theta_alpha

# dispersion measurement -----------------------------------------------------

def measure_percentile_angle_ex_kappa(kappa, percentile=0.95, n_pts=int(1e3)):
    '''
    Creates a random distribution with large `n` and dispersion `kappa`,
    and finds the angle that encompasses `confid` * 100% of the pts. 
    
    Parameters
    ----------
    kappa : float
        dispersion parameter of points
    percentile : float
        proportion of points to include inside angle
    n_pts : int
        number of points to create in distribution

    Returns
    -------
    theta_percentile : float
        angle that encompasses `percentile` * 100% of points `P_i`    
    
    Notes
    -----
    In this case I think theta_alpha does not depend on n, so it is best
    to estimate with a large n for better accuracy.
    '''
    perpvec = np.array((0., 0., 1.))
    P_i = vmf_rvs(perpvec, kappa, n_pts=n_pts)
    return measure_percentile_angle(P_i, percentile=percentile)

def measure_percentile_angle(P_i, percentile=0.95):
    '''
    Find the angle that encompasses percentile * 100% of the pts in `P_i`.
    
    Parameters
    ----------
    P_i : array_like, shape (n_pts, 3)
        x,y,z vectors of unit length points on a sphere
    percentile : float
        proportion of points to include inside angle

    Returns
    -------
    theta_percentile : float
        angle that encompasses `percentile` * 100% of points `P_i`
    '''
    # find angles between sample mean and each pt
    muhat = mean_dir(P_i)
    angles_to_mean = np.arccos(np.dot(muhat, P_i.T))
    
    # take percentile of angles
    sort_idxs = np.argsort(angles_to_mean)
    sorted_angles = angles_to_mean[sort_idxs]
    n_pts = P_i.shape[0]
    pc_pts = np.arange(1, n_pts+1) / float(n_pts)
    #... percentiles of each sorted pt
    pt_below = np.nonzero(pc_pts < percentile)[0][-1]
    pt_above = np.nonzero(pc_pts > percentile)[0][0]
    pc_below = pc_pts[pt_below]
    pc_above = pc_pts[pt_above]
    theta_below = sorted_angles[pt_below]
    theta_above = sorted_angles[pt_above]

    # now weight two angles by how close they are to confid
    #- i.e. linearly interpolate between them
    theta_percentile = ((percentile - pc_below) * theta_above + \
			(pc_above - percentile) * theta_below) / \
			(pc_above - pc_below)
    return theta_percentile

# ------------------ visualization ---------------------

def plot_pts(pts=None, mu=None):
    p3d = mlab.pipeline

    # plot unit sphere
    sphere_pts = point_primitives.sphere()
    sphere1 = mlab.mesh(*sphere_pts)
    sphere1.actor.mapper.scalar_visibility = False
    sphere1.actor.property.color = (241/255., 233/255., 199/255.)
    sphere1.actor.property.opacity = 0.25

    # plot pts
    if pts != None:
        x = pts[...,0]
        y = pts[...,1]
        z = pts[...,2]
        src = p3d.scalar_scatter(x,y,z)
        glyphs = p3d.glyph(src)
        glyphs.actor.property.color = (0.0, 0.0, 1.0)
        glyphs.glyph.glyph_source.glyph_source = \
            glyphs.glyph.glyph_source.glyph_list[4]
        glyphs.glyph.glyph.scale_factor = 0.1

    # plot mean
    if mu != None:
        zeros = np.zeros_like(mu[0])[..., np.newaxis]
        mlab.quiver3d(zeros, zeros, zeros,
                      mu[0, np.newaxis], mu[1, np.newaxis], mu[2, np.newaxis])

def plot_circle(mu, theta, scalars=None, scalar_max=None,
		color=None, radius=0.01,
		resolution=50.):
    x,y,z = 0,1,2
    phi_prime = np.linspace(0, 2 * np.pi, resolution)
    # start with defining cone around (0,0,1)
    origin = np.array((0.,0.,1.))
    P_j = np.array([vectors.rotate_by_angles(origin, theta, q)
                    for q in phi_prime]).T
    theta, phi = convert_cartesian_to_polar(mu)
    P_k = vectors.rotate_by_angles(P_j, theta, phi)
	    
    p3d = mlab.pipeline
    if scalars == None:
	tube = p3d.tube(p3d.line_source(P_k[x], P_k[y], P_k[z]))
    else:
	if type(scalars) == type(1) or type(scalars == type(1.)):
	    #'arraying %d' % (
	    scalars *= np.ones_like(P_k[0])
	    #'assigning scalars 
	    tube = p3d.tube(p3d.line_source(P_k[x],
					    P_k[y],
					    P_k[z], scalars))
	    tube.filter.radius = radius
	    surf = p3d.surface(tube)
    if color != None:
	surf.actor.actor.property.color = color
    if scalar_max != None:
	mm = surf.module_manager.scalar_lut_manager
	mm.data_range = np.array([0, scalar_max])

# def test_theta_P_single_dist(n_samples=1e5):
#     # generate a vMF distribution with a random mu (mean) and kappa (spread)
#     theta = np.random.uniform(low = -np.pi/2., high = np.pi/2.)
#     phi = np.random.uniform(low = -np.pi, high = np.pi)
#     mu = convert_polar_to_cartesian(theta, phi)
    
#     kappa = np.random.uniform(low=1., high=5.)
#     P_i = vmf_rvs(mu, kappa, n_pts=n_samples)
    
#     # estimate kappa from the distribution, and use the real mean
#     #k_hat = estimate_kappa(P_i, pop_mean=mu)
    
#     # calculate the confidence cone angle
#     theta_confid = estimate_confid_angle(P_i, 0.5, mu=mu)
    
#     # count how many of the variates lie inside the cone
#     angles = np.arccos(np.dot(P_i, mu))
#     return P_i, mu, theta_confid, (angles < theta_confid).sum()

# def make_sphere_pts(n=1000):
#     elevation = np.random.uniform(low=-1.0, high=1.0, size=1000) * np.pi * 2
#     azimuth = np.random.uniform(low=-1.0, high=1.0, size=1000) * np.pi * 2
#     pts = np.array([np.cos(azimuth) * np.cos(elevation),
#                     np.sin(azimuth) * np.cos(elevation),
#                     np.sin(elevation)])
#     return pts

# def points_on_sphere(n):
#     # some weird distribution - not uniform
#     n = float(n) # in casewe got an int which we surely got
#     pts = []
#     inc = np.pi * (3 - np.sqrt(5))
#     off = 2 / n
#     for k in range(int(n)):
#         y = k * off - 1 + (off / 2)
#         r = np.sqrt(1 - y * y)
#         phi = k * inc
#         pts.append([np.cos(phi) * r, y, np.sin(phi) * r]) 
#     return np.asarray(pts)

# def wrap_test_theta(n_repeats = 100):
#     total = 0
#     for i in xrange(n_repeats):
#         total += test_theta_P()
#     return total / float(n_repeats)

# def debug_Stephens1962_nomograms(resolution=3.):
#     alphas = [0.1,]
#     n_alphas = len(alphas)
#     Ns = range(3,17) + [18,20,25,30]
#     #Ns = [12,16]
#     n_Ns = len(Ns)
#     coords = np.ma.empty((n_alphas, n_Ns, max(Ns) * resolution))
#     coords.mask = ~np.isnan(coords)
#     R0s = np.ma.empty((n_alphas, n_Ns, max(Ns) * resolution))
#     R0s.mask = ~np.isnan(R0s)
#     for i_alpha, alpha in enumerate(alphas):
#         for i_N, N in enumerate(Ns):
#             print 'processing N=%d' % (N)
#             Rstars = np.linspace(0, N, N * resolution)
#             for i_Rstar, Rstar in enumerate(Rstars):
#                 if Rstar == N:
#                     Rstar -= 1e-1
#                     #print Rstar
#                 print 'processing Rstar=%5.3f' % (Rstar)
#                 coords[i_alpha, i_N, i_Rstar] = Rstar
#                 R0_low = 0
#                 R0_high = max(Ns)
#                 R0 = calculate_R0(N, Rstar, alpha)#, R0_low, R0_high)
#                 R0s[i_alpha, i_N, i_Rstar] = R0
#                 R0s.mask[i_alpha, i_N, i_Rstar] = False
#     return alphas, Ns, coords, R0s
    
# def debug2_Stephens1962_nomograms(N, Rstar, alpha, resolution=5., figure=2):
#     R0s = np.linspace(0, N, N*resolution)
#     fns = np.empty_like(R0s)
#     for i_R0, R0 in enumerate(R0s):
#         fns[i_R0] = prob_fn(R0, N, Rstar, alpha)
#     plt.figure(figure)
#     useful = fns < 1.
#     #plt.plot(R0s, fns, 'k.-')
#     slope = np.diff(fns[useful])
    
#     plt.plot(R0s[useful][:-1][slope < 0],
#              fns[useful][:-1][slope < 0], 'b')
#     plt.plot(R0s[useful][:-1][slope > 0],
#              fns[useful][:-1][slope > 0], 'r')
#     plt.draw()
#     return R0s, fns, slope

# def gradient_descent(fn, start=0, args=(), init_delta=1.,
#                      gamma=0.01, err=1e-8, full_output=False):
#     # let's start at 0 always - should be okay for this function
#     pt = start
#     delta = init_delta
#     while np.abs(delta) > err:
#         if full_output:
#             print 'pt: %10.8f' % (pt),
#         grad = (fn(pt + delta, *args) - fn(pt, *args)) / delta
#         if full_output:
#             print 'grad: %10.2g' % (grad),
#         delta = -gamma * grad / np.maximum(fn(pt, *args), 1.)
#         if full_output:
#             print '  delta: %5.2f' % (delta)
#         pt += delta
#     return pt

# def calc_confid_angle_Gidskehaug(kappa, alpha):
#     k = kappa
#     arccos = np.arccos
#     exp = np.exp
#     log = np.log
#     return arccos(1/k * log(exp(k) - (1-alpha) * (exp(k) - exp(-kappa))))

