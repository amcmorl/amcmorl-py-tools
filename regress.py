# fcdf is in different places in different scipy versions
from scipy.stats import f
fcdf = f.cdf
import numpy as n
from LinearAlgebra import linear_least_squares, inverse

def regress(y,X):
    """
    Perform a multiple linear regression of y onto X.  X is a num
    observations by (num variables +1) array.  X should contain a
    column of ones to generate the intercept.  y is a num observations
    array of independent variables

    return value B, residuals, stats

    B: the regression coeffients;  Ypred = B.T * X.T
    residuals: array( y - Ypred.T)
    stats = Rsquared, F, p
    
    """

    # regression coeffs are given by (Xt*X)-1*Xt*y
    N = X.shape[0]
    y.shape = N, 1
    X = n.mat(X)
    Y = n.mat(y)
    Xt = X.T
    Xt_X_i = n.mat( n.linalg.inv(Xt * X) )
    B = Xt_X_i * Xt * Y

    Ypred = B.T * Xt
    residuals = n.array(Y - Ypred.T)
    CF = N * y.mean()**2     # correction factor

    SStotal = float(Y.T * Y - CF)
    SSregress =  float(B.T * Xt * Y - CF)
    SSerror =  SStotal - SSregress

    Rsquared = SSregress / SStotal

    dfTotal = N - 1
    dfRegress = len(B) - 1
    dfError = dfTotal - dfRegress

    F = SSregress / dfRegress / (SSerror / dfError)
    prob = 1 - fcdf(F, dfRegress, dfError)

    stats = Rsquared, F, prob
    return B, residuals, stats

if __name__ == "__main__":
    fhx = file('independent.dat', 'w')
    fhy = file('dependent.dat', 'w')

    # X is number observations by (number vars + constant)
    X = n.random.randn( 200,6 )    

    # The last column should be ones for the constant term
    X[:,5] = n.ones((200,), dtype=float)

    y = n.random.randn( 200 )   # the dependent var
    
    #save the data to file for comparison with matlab's regress
    for row in X:
        print >>fhx, ' '.join([str(val) for val in row])
        print >>fhy, '\n'.join([str(val) for val in y])
    fhx.close()
    fhy.close()

    # c are the coeffiecients
    # ss is the sum of squared residuals
    # sv are the singular values of X
    c, ss, rank, sv = n.linalg.lstsq(X,y)

    # below is just a demo to show how the results can be used to generate
    # a prediction
    p =  n.dot(X,c)
    ss1 = ((y - p)**2).sum() # should equal 'ss'

    B, resid, stats = regress(y,X)
    print stats  # validated against multiple_regress_demo.m
