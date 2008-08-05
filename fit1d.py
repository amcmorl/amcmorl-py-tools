import numpy as _n
from scipy.optimize import leastsq
from scipy.interpolate import splprep, splev


# --- double_boltzman --------------------------------------

def double_boltzman_by_centre(centre, height, width, curve, \
                              length=None, r=None):
    '''wrapped for double_boltzman simplifying call for symmetric
    function with centre, multiplier (height iff function reaches 1),
    width between pivot points, curve and optionally x vals given by r
    or in a 1-D array of length.
    '''
    t1,t2 = (centre - width/2., centre + width/2.)
    return double_boltzman( t1, curve, t2, curve, height, \
                            length=length, r=r )


def double_boltzman(t1, tau1, t2, tau2, height=1, \
                    length=None, r=None):
    '''Returns a double boltzman function

    Usage:
      array = double_boltzman( t1, tau1, t2, tau2, height=1,
                               length=None, r=None )

    where:
      t1, t2     = pivot points
      tau1, tau2 = curviness on left and right respectively
      height     = multiplier (may not be actual height if
                   function doesn''t reach 1)
      length     = length of 1-D array to hold
       OR
      r          = x vals

    '''
    if r is None:
        if length is None:
            raise ValueError("One of r or length must be specified.")
        r = _n.arange(length)

    frac1denom = 1+_n.exp(-(r-t1)/tau1)
    frac2denom = 1+_n.exp((r-t2)/tau2)
    
    return height * (1/frac1denom) * (1/frac2denom)


def _double_boltzman_erf(p, d, r, tau):
    '''error function for Boltzman fitting
    returns residuals: data - double-boltzman

    fitting parameters (p) are:
      centre (c),
      height (h) and
      width between pivot points (w)

    other parameters are curvature (k) which is fixed prior to fitting
    (and doesn''t affect width)
    '''
    (c,h,w) = p
    t1,t2 = (c - w/2., c + w/2.)
    f = double_boltzman( t1, tau, t2, tau, h, r=r )
    return d - f


def fit_double_boltzman(data, p0=None, r=None, curve=1.):
    '''Fits a double boltzman function (using scipy.optimize.leastsq)
    to 1-D data. Returns p of best fit: centre, height, width'''

    if r == None:
        r = _n.arange( len( data ) )
    else:
        if len(r) != len(data):
            raise ValueError("Data and r (x-values) are not equal length.")

    if p0 == None:
        p0 = (1.,5.,1.)
    lsq = leastsq( _double_boltzman_erf, p0, args=(data, r, curve))
    if lsq[-1] == 1:
        return lsq[0]
    else:
        return None, None, None


def double_boltzman_fwhm(centre, mult, ipd, curve, length=None, r=None, ):
    '''Finds (numerically) the fwhm of the boltzman equation defined by
    parameters p, curve and optionally with x vals (r)

    Usage:
      fwhm = double_boltzman_fwhm(p, curve, r=None)

    Where:
      centre = centre of symmetry
      height = multiplier for function - not necessarily max point
      intra-pivot distance = distance between pivot points 
      curve  = curviness
      r      = x vals'''

    if r == None:
        if length == None:
            raise ValueError("One of r or length must be specified.")
        r = _n.arange(length)

    fvals = double_boltzman_by_centre(centre, mult, ipd, curve, \
                                      length=length, r=r )    
    # max is at centre, for symmetric function (curves the same)
    hmx = fvals.max() / 2.

    # quick & dirty method -> fwhm = 2 * dist( centre to 1 half-max )
    diffs = abs(fvals - hmx)
    return 2 * abs(centre - r[_n.where(diffs == diffs.min())])

# --- exponential decay ------------------------------------------

def exp_decay(a, c, length=None, r=None):
    ''' returns an exponential decay function in the form

    y = (1 - a)e^{\frac{x}{c}} + a

    which is a standard, normalized exponential from 1 to a

    length     = length of 1-D array to hold
       OR
    r          = x vals
    '''
    
    if r is None:
        if length is None:
            raise ValueError("One of r or length must be specified.")
        r = _n.arange(length)

    return (1 - a) * _n.exp(-r/c) + a

def _exp_decay_erf(p, d, r):
    '''error function for exponential decay fitting
    returns residuals: data - exp_decay

    fitting parameters (p) are:
      a - x-asymptote
      c - decay constant
    '''
    (a,c) = p
    f = exp_decay(a,c,r=r)
    return d - f

def fit_exp_decay(data, p0=None, r=None):
    '''Fits an exponential decay function (using scipy.optimize.leastsq) to
    to 1-D data. Returns p of best fit: asymptote, decay const.'''

    if r == None:
        r = _n.arange( len( data ) )
    else:
        if len(r) != len(data):
            raise ValueError("Data and r (x-values) are not equal length.")

    if p0 == None:
        p0 = (0., 1.)
    lsq = leastsq(_exp_decay_erf, p0, args=(data, r))
    if lsq[-1] == 1:
        return lsq[0]
    else:
        return None, None, None
    
# --- exponential ------------------------------------------

def exp_cdf(c, length=None, r=None):
    ''' returns an exponential decay function in the form

    y = e^{\frac{x}{c}}

    length     = length of 1-D array to hold
       OR
    r          = x vals
    '''
    if r is None:
        if length is None:
            raise ValueError("One of r or length must be specified.")
        r = _n.arange(length)
    return 1 - _n.exp(-r/c)

def _exp_cdf_erf(p, d, r):
    '''error function for exponential decay fitting
    returns residuals: data - exp_decay

    fitting parameters (p) are:
      a - x-asymptote
      c - decay constant
    '''
    c = p
    f = exp_cdf(c,r=r)
    return d - f

def fit_exp_cdf(data, p0=None, r=None):
    '''Fits an exponential cdf function (using scipy.optimize.leastsq) to
    1-D data. Returns p of best fit.
    '''
    if r == None:
        r = _n.arange( len( data ) )
    else:
        if len(r) != len(data):
            raise ValueError("Data and r (x-values) are not equal length.")
    if p0 == None:
        p0 = 1.
    lsq = leastsq(_exp_cdf_erf, p0, args=(data, r))
    if lsq[-1] == 1:
        return lsq[0]
    else:
        return None, None   
