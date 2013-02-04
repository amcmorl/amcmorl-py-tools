import numpy as np
from scipy import stats

def _get_z(alpha):
    return stats.distributions.norm.ppf(1 - alpha/2.)

def agresti_coull_interval(p, n, alpha):
    '''
    Returns the Agresti-Coull interval for approximating the binomial 
    confidence interval.
    
    From Wikipedia [1]:
    
    .. math:: 
    
        \tilde{n} = n + z^2_{1-\alpha/2}
        X = pn
        \tilde{p} = \frac{X + z^2_{1-\alpha/2}/2}{\tilde{n}}
        
    and the confidence interval is given by:
    
    .. math::
    
        \tilde{p} \pm z_{1-\alpha/2} \sqrt{\frac{\tilde{p}(1-\tilde{p})}{\tilde{n}}}
    
    References
    ----------
    
    [1] http://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
    '''
    z = _get_z(alpha)
    n_ = n + z**2
    X = p * n
    p_ = (X + z**2 / 2.) / float(n_)
    i = np.array([-1,1])
    return p_ + i * z * np.sqrt((p_ * (1 - p_)) / n_)
