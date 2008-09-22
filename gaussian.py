from scipy.optimize import leastsq
import numpy as np
import unittest

# gaussians take the form
#
#  I(r) = Ae^(-r**2/k**2)
#
# where A is the peak
#       k determines the rate of curve such that:
#
#                     FWHM
#        k =  --------------------
#             (2 * sqrt( ln( 2 )))
#
#        or conversely:
#
#        FWHM = 2k * sqrt( ln( 2 ) )

# general service routines

def k2fwhm(k):
    '''converts k value (see above) to fwhm'''
    return 2 * k * np.sqrt( np.log( 2 ) )

def fwhm2k(fwhm):
    '''converts fwhm value to k (see above)'''
    return fwhm/(2 * np.sqrt( np.log( 2 ) ) )

def kerr2size(k,err):
    return 2 * k * np.sqrt( np.log( 1 / err ) )

# 1d radial gaussian generation and fitting

def erfr(p, I, r):
    '''gaussian error function
    returns residuals after gaussian approximation:
    data - gaussfit'''
    
    (A, k) = p
    return I - A * np.exp( -r**2 / k**2 )

def fitgauss1dr(data, p0=None, r=None):
    '''Fits a gaussian (using scipy.optimize.leastsq)
    to 1d data. Returns the tuple (A, fwhm)'''

    if r == None:
        r = np.arange( len( data ) )
    if p0 == None:
        p0 = (1,1)
    plsq = leastsq(erfr, p0, args=(data, r))
    #^used to be erf here - think it's wrong
    A = plsq[0][0]
    fwhm = k2fwhm(plsq[0][1])
    return A, fwhm

def gauss1dr(r, A, fwhm):
    '''returns the radial 1d gaussian given by A (peak) & fwhm
    '''
    return A * np.exp( -r**2 / fwhm2k(fwhm)**2 )

# 1d gaussian (non-radial) generation and fitting

def erf(p, I, r):
    '''gaussian error function
    returns residuals after gaussian approximation:
    data - gaussfit'''
    
    (A, k, c) = p
    return I - A * np.exp( -(r - c)**2 / k**2 )

def fitgauss1d(data, p0=None, r=None):
    '''Fits a gaussian (using scipy.optimize.leastsq)
    to 1d data.

    (A, fwhm, c) = fitgauss1d( data[, p0][, r] )

    where A = peak
          fwhm = full-width, half-max
          c = centre
          data = 1d profile to fit
          p0 = tuple of initial estimates
          r = positions of data
    '''    
    if r == None:
        r = np.arange( len( data ) )
    if p0 == None:
        p0 = (1,1,1)
    plsq = leastsq(erf, p0, args=(data, r))
    A = plsq[0][0]
    fwhm = k2fwhm(plsq[0][1])
    c = plsq[0][2]
    return A, fwhm, c

def gauss1d(r, A, fwhm, c):
    '''returns the 1d gaussian given by
    A (peak), fwhm, and c (centre)
    at positions given by r
    '''
    return A * np.exp( -(r-c)**2 / fwhm2k( fwhm )**2 )

# 1d with offset

def erf_const(p, I, r):
    '''gaussian error function
    returns residuals after gaussian approximation:
    data - gaussfit'''
    
    (A, k, c, d) = p
    return I - (A * np.exp( -(r - c)**2 / k**2 ) + d)

def fitgauss1d_const(data, p0=None, r=None):
    '''Fits a gaussian (using scipy.optimize.leastsq)
    to 1d data.

    (A, fwhm, c) = fitgauss1d( data[, p0][, r] )

    where A = peak
          fwhm = full-width, half-max
          c = centre
          data = 1d profile to fit
          p0 = tuple of initial estimates
          r = positions of data
    '''    
    if r == None:
        r = np.arange( len( data ) )
    if p0 == None:
        p0 = (1,1,1,1)
    plsq = leastsq(erf_const, p0, args=(data, r))
    A, fwhm_, c, d = plsq[0]
    fwhm = k2fwhm(fwhm_)
    return A, fwhm, c, d

def gauss1d_const(r, A, fwhm, c, d):
    '''returns the 1d gaussian given by
    A (peak), fwhm, and c (centre)
    at positions given by r
    '''
    return A * np.exp( -(r-c)**2 / fwhm2k( fwhm )**2 ) + d

# 2d gaussian construction

def gauss2d(fx, fy, err=0.01):
    '''Returns a volume containing a 2-D gaussian field
    Usage:
      g = gauss2d(fx, fy, err=allowed_error)

    Where:
      fx, fy are the fwhms in the x and y direction
      err is the maximum allowed error
        - this determines the size of the result
        (i.e. how small the tails have to get at the edges)
    AJC McMorland 13-9-2006
    '''
    ks = np.array( [fwhm2k(i) for i in (fx, fy)] )
    dims = np.array( [kerr2size(i,err) for i in list( ks )] ).round()
    dimvals = abs( np.indices( dims ) - \
                   (dims[..., np.newaxis, np.newaxis] - 1) / 2) \
                   / ks[..., np.newaxis, np.newaxis]
    Is = np.exp(-(dimvals**2).sum(0))
    return Is / Is.sum()

# 3d gaussian construction

def gauss3d(fx, fy, fz, err=0.01):
    '''Returns a volume containing a 3d gaussian field
    with the supplied fwhms fx, fy and fz.
    The size of volume is defined by err, which is the
    smallest value the volume should contain.'''

    ks = np.array( [fwhm2k(i) for i in (fx,fy,fz)] )
    dims = np.array( [kerr2size(i,err) for i in list( ks )] ).round()
    dimvals = abs(np.indices( dims ) - \
                  (dims[..., np.newaxis, np.newaxis, np.newaxis] - 1) / 2) \
                  / ks[..., np.newaxis, np.newaxis, np.newaxis]
    Is  = np.exp(-(dimvals**2).sum(0))
    return Is / Is.sum()


class TestGaussianFunctions(unittest.TestCase):

    def test_k2fwhm_fwhm2k(self):
        a = 1
        r = 2*a*np.sqrt( np.log( 2 ) )
        self.assertTrue( np.allclose(k2fwhm(a), r) )
        self.assertTrue( np.allclose(fwhm2k(r), a) )
        self.assertTrue( np.allclose(a, fwhm2k(k2fwhm(a))))

    def test_kerr2size(self):
        a, err = (1., 0.1)
        r = 2 * np.sqrt( np.log( 1. / 0.1 ) )
        self.assertTrue( np.allclose( kerr2size(a,err), r ) )


def test():
    suite = unittest.TestLoader().loadTestsFromTestCase( \
        TestGaussianFunctions)
    unittest.TextTestRunner(verbosity=2).run(suite)
