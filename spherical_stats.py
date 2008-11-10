import numpy as np
from enthought.mayavi import mlab
from scipy.optimize import leastsq
from scipy.integrate import quad
import vectors

# -------------------------- von-Mises-Fisher dist ---------------------------

def estimate_kappa(P_i, mean_dir=None):
    '''Returns the maximum likelihood estimate of the spread parameter (kappa)
    of a von-Mises-Fisher distribution.

    Parameters
    ----------
      pts : array_like, shape (n, 3)
        n unit vectors for which to find \kappa
      mean_dir : array_like, shape (3,)
        vector of direction cosines for known population mean

    Returns
    -------
      kappa : scalar
        spread estimate for sample, assuming von-Mises-Fisher distribution
      
    Notes
    -----
    Code is based on equations in Statistical Analysis of Spherical Data,
    1987 by NI Fisher, T Lewis, and BJJ Embleton, section 5.3.2(iv).

    coth(\kappa) = 1/\kappa = R/n

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
    S = np.sum(P_i, axis=0)    # 3.2
    R = np.sqrt(np.sum(S**2)) # 3.3
    k0 = 1.
    lsq = leastsq(to_min, k0, args=(R, n))
    if lsq[-1] == 1:
        return lsq[0]
    
def to_min(k, R, n):
    return 1/np.tanh(k) - 1/k - R/n

def C_F(kappa):
    '''From Statistical Analysis of Spherical Data,
    1987 by NI Fisher, T Lewis, and BJJ Embleton, equation 4.21'''
    return kappa / (2 * np.pi * (np.exp(kappa) - np.exp(-kappa)))

def vmf_pdf(theta, kappa):
    '''Returns the value of the probability density function for the
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
    return quad(vmf_pdf, 0, theta, args=(kappa))

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

    pts = vectors.convert_angles_to_components(colat, longit)
    alpha, beta = vectors.convert_components_to_angles(mu)
    return vectors.rotate_by_angles(pts, alpha, beta)

# --------------------------- testing code ---------------------------

def make_sphere_pts(n=1000):
    elevation = np.random.uniform(low=-1.0, high=1.0, size=1000) * np.pi * 2
    azimuth = np.random.uniform(low=-1.0, high=1.0, size=1000) * np.pi * 2
    pts = np.array([np.cos(azimuth) * np.cos(elevation),
                    np.sin(azimuth) * np.cos(elevation),
                    np.sin(elevation)])
    return pts

def plot_pts(pts):
    x = pts[0]
    y = pts[1]
    z = pts[2]
    mlab.points3d(x,y,z)

def points_on_sphere(n):
    n = float(n) # in case we got an int which we surely got
    pts = []
    inc = np.pi * (3 - np.sqrt(5))
    off = 2 / n
    for k in range(int(n)):
        y = k * off - 1 + (off / 2)
        r = np.sqrt(1 - y * y)
        phi = k * inc
        pts.append([np.cos(phi) * r, y, np.sin(phi) * r]) 
    return np.asarray(pts)

