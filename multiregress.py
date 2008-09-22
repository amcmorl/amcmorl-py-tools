import numpy as np

def multiregress(xs, y):
    '''Multiple linear regression with statistics.
    
    Calculates the partial regression co-efficients b_{1..k} in the
    equation Y = b_0 + b_1.x_1 + b_2.x_2 + ... + b_k.x_k

    Returns the standard errors of the regression coefficients:
    s_{bY1}, s_{bY2}, as described in Biometry, Sokal and Rohlf, 3rd Ed.,
    box 16.2.

    Parameters
    ----------
    xs : array, shape (n, k)
        k is number of variates, n is number of observations
    y : array, shape (n,)

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
    xs_c = np.concatenate((np.ones((n,1)), xs), axis=1)
    b, resid, rank, s = np.linalg.lstsq(xs_c, y)
    s_y = np.sqrt(resid / (n - k - 1.)) # unexplained SD
    G_mult = np.linalg.inv( np.cov( xs.T ) * (n - 1) )
    se = s_y * np.sqrt(G_mult[np.eye(k, dtype=bool)])
    return b, se

