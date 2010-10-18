import numpy as np

def _chkx(x):
    if type(x) == int:
        return np.arange(x)
    elif type(x) == np.ndarray:
        return x
    elif (type(x) == list) or (type(x) == dict):
        return np.asarray(x)
    else: raise ValueError("x must be int or array_like")

kerr2size = lambda k, err : 2 * k * np.sqrt(np.log(1. / err))
k2fwhm    = lambda k : 2 * k * np.sqrt(np.log(2.))
fwhm2k    = lambda fwhm : fwhm / (2 * np.sqrt(np.log(2.)))
gauss1d   = lambda x, A, fwhm, c : A * np.exp( -(_chkx(x) - c)**2 \
                                                    / fwhm2k(fwhm)**2 )

def gauss2d(fx, fy, err=0.01):
    '''Returns a volume containing a 2-D gaussian field.

    Parameters
    ----------
    fx, fy : float
      the fwhms in the x and y directions
    err : float
      the maximum allowed error, determines the size of the result
      (i.e. how small the tails have to get at the edges)

    Notes
    -----
    AJC McMorland 13-9-2006
    '''
    ks = np.array( [fwhm2k(i) for i in (fx, fy)] )
    dims = np.array( [kerr2size(i,err) for i in list( ks )] ).round()
    dimvals = abs( np.indices(dims) - \
                   (dims[..., None, None] - 1) / 2) \
                   / ks[..., None, None]
    Is = np.exp(-(dimvals**2).sum(0))
    return Is / Is.sum()

def gauss3d(fx, fy, fz, err=0.01):
    '''Returns a volume containing a 3d gaussian field.

    Parameters
    ----------
    fx, fy, fz : float
      the fwhms in the x and y directions
    err : float
      the maximum allowed error, determines the size of the result
      (i.e. how small the tails have to get at the edges)
    '''
    ks = np.array( [fwhm2k(i) for i in (fx,fy,fz)] )
    dims = np.array( [kerr2size(i,err) for i in list( ks )] ).round()
    dimvals = abs(np.indices( dims ) - \
                  (dims[..., None, None, None] - 1) / 2) \
                  / ks[..., None, None, None]
    Is  = np.exp(-(dimvals**2).sum(0))
    return Is / Is.sum()

def logistic(x, m, a, b):
    '''
    Logistic, s-curve, or squashing function.

    Parameters
    ----------
    x : array_like
      x values to calculate for
    m : float or int
      mid point of curve, i.e. x @ which f(x) = a/2.
    a : float or int
      amplitude, f(+inf) = a
    b : float or int
      baseline, f(0) = b

    Returns
    -------
    f(x) : ndarray
      s-curve of `x` data
    '''
    m = float(m)
    a = float(a)
    b = float(b)
    c = -m / np.log(1 / (a / b - 1))
    d = c * np.log(a / b - 1)
    return a / (1 + np.exp((-x + d) / c))

def squash(x, sqp):
    '''
    Alternative signature for squashing function, using a parameter dictionary.

    Parameters
    ----------
    x : array_like
      x values to calculate for
    sqp : dict
      parameters for squashing function, keys are:
      m - mid-point, x at which y=a/2
      a - amplitude, max y value, at x=+inf
      b - baseline, y value at x=0
      
    Returns
    -------
    f(x) - ndarray
      squashed `x` data
    '''
    return logistic(x, sqp['m'], sqp['a'], sqp['b'])

def set_squash_defaults(sqp, x):
    '''
    Set reasonable defaults on squashing function parameters. Works inplace.

    Parameters
    ----------
    sqp : dict
      parameters for squashing function, see `squash` definition
    x : array_like
      independent variables
    '''
    sqp.setdefault('m', (np.max(x) + np.min(x)) / 2.)
    sqp.setdefault('a', np.max(x))
    sqp.setdefault('b', 0.05 * np.max(x))

def dbl_boltzman(x, t1, tau1, t2, tau2, h):
    '''
    Parameters
    ----------
    x : ndarray
      independent variable
    t1, t2 : float
      times of inflection points
    tau1, tau2 : float
      curvatures
    h : float
      amplitude
    '''
    x = _chkx(x)
    denom1 = 1 + np.exp(-(x - t1) / tau1)
    denom2 = 1 + np.exp((x - t2) / tau2)
    return h * (1 / denom1) * (1 / denom2)
    
def dbl_boltzman_by_centre(x, cen, w, tau, h):
    t1, t2 = cen - w / 2., cen + w / 2.
    return dbl_boltzman(x, t1, tau, t2, tau, h)

def double_boltzman_fwhm(x, cen, A, ipd, tau):
    '''Finds (numerically) the fwhm of the boltzman equation defined by
    parameters p, curve and optionally with x vals (r)

    Parameters
    ----------
    x : array_like or int
        independent variable
    centre : float
        centre of symmetry
    intra-pivot distance : float
        distance between pivot points 
    tau : float
        curviness
    A : float
        multiplier for function, not necessarily max point
    '''
    x = _chkx(x)
    fvals = dbl_boltzman_by_centre(x, cen, ipd, tau, A)    
    # max is at centre, for symmetric function (curves the same)
    hmx = fvals.max() / 2.

    # quick & dirty method -> fwhm = 2 * dist( centre to 1 half-max )
    diffs = np.abs(fvals - hmx)
    return 2 * np.abs(cen - x[np.where(diffs == diffs.min())])

exp_decay = lambda x, a, c : (1 - a) * np.exp(-_chkx(x) / c) + a
exp_cdf   = lambda x, c : 1 - np.exp(-_chkx(x) / c)
cos_cdf   = lambda x, a, b, c : a * np.cos(b * _chkx(x)) + c
