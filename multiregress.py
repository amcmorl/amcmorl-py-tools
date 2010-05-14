import numpy as np

def multiregress(xs, y, with_err=True, add_const=True):
    '''Multiple linear regression with statistics.
    
    Calculates the partial regression co-efficients b_{0..k} in the
    equation Y = b_0 + b_1.x_1 + b_2.x_2 + ... + b_k.x_k

    Optionally returns the standard errors of the regression coefficients:
    s_{bY1}, s_{bY2}, as described in Biometry, Sokal and Rohlf, 3rd Ed.,
    box 16.2.

    Parameters
    ----------
    xs : array, shape (n, k)
        n is number of observations, k is number of variates, 
    y : array, shape (n,)
    with_err : bool
        calculate standard errors of regression coefficients
    add_const : bool
        include constant term in regression
    
    Returns
    -------
    b : array, shape (k + 1)
        partial regression co-efficients, first element is constant term
    errs : array, shape (k)
        standard errors of regression co-efficients in b,
        except constant term b_0
    '''
    # add ones for linear regression
    n, k = xs.shape
    if add_const:
        xs_c = np.concatenate((np.ones((n,1)), xs), axis=1)
    else:
        xs_c = xs
    b, resid, rank, s = np.linalg.lstsq(xs_c, y)
    if with_err:
        # resid is the 'unexplained sum-of-squares', a scalar
        s_y = np.sqrt(resid / (n - k - 1.)) # unexplained SD
        G_mult = np.linalg.inv( np.cov( xs.T ) * (n - 1) )
        se = s_y * np.sqrt(G_mult[np.eye(k, dtype=bool)])
        return b, se
    else:
        return b

def rsq(x, y, add_const=True):
    # calculate regression coefficients
    b = multiregress(x, y, with_err=False, add_const=add_const)
    b_nc = b[1:] # ignore constant term
    
    # standardize partial regression coefficients
    s_y = np.std(y, ddof=1)
    s_x = np.std(x, ddof=1, axis=0)    
    b_prime = b_nc * s_x / s_y

    # calculate correlation coefficients
    r_xy = np.corrcoef(y.T, x.T)[0,1:]
    rsq = np.sum(r_xy * b_prime)
    return rsq

def rsq_from_b(x, y, b):
    '''calculate r-squared value without recalculating regression.'''
    b_nc = b[1:]
    # standardize partial regression coefficients
    s_y = np.std(y, ddof=1)
    s_x = np.std(x, ddof=1, axis=0)    
    b_prime = b_nc * s_x / s_y

    # calculate correlation coefficients
    r_xy = np.corrcoef(y.T, x.T)[0,1:]
    rsq = np.sum(r_xy * b_prime)
    return rsq
    
