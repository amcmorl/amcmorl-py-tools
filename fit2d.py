from numpy import array, arange, where, linspace

graphic = True
try:
    from pylab import figure, clf, plot
except:
    graphic = False
from scipy.interpolate import splprep, splev
    
def midcrossings(a, b = None, thresh = 1e-3, k = 5):
    '''usage res = midcrossings([x,] y)

    returns fwhm of y (a function with a single maximum value)
    from spline-interpolated midpoint crossings'''
    if b == None:
        y = array( a )
        x = arange(y)
    else:
        x = array( a )
        y = array( b )
    try:
        assert x.shape[0] == y.shape[0]
    except AssertionError:
        print "x and y must be same length"
        return None
    maxind = where( y==y.max() )[0].flatten()[0] # uses only first max pt

    (y1, y2) = y[:maxind], y[maxind:]
    (x1, x2) = x[:maxind], x[maxind:]

    s = 1.
    nest = -1
    interpvals = linspace( 0, 1, 251 )
      # 251 simply gives enough points to give one point close to 0.5

    print "thresholding to:", thresh
    # lower half
    nob1 = where( y1 > thresh )
    # ^ need to ignore baseline values when fitting splines
    y1tofit = y1[nob1] / y1[nob1].max()
    
    tckp, u = splprep( [y1tofit, x1[nob1]], s=s, k=k, nest=nest )
    y1p, x1p = splev( interpvals, tckp )
    dtohalf = abs(y1p - 0.5)
      # 0.5 because want width at _half_ max
    closest = where( dtohalf == dtohalf.min() )
    lowval = x1p[closest]

    # upper half
    nob2 = where( y2 > thresh )
    y2tofit = y2[nob2] / y2[nob2].max()
    tckp, u = splprep( [y2tofit, x2[nob2]], s=s, k=k, nest=nest )
    y2p, x2p = splev( interpvals, tckp )
    dtohalf = abs(y2p - 0.5)
    closest = where( dtohalf == dtohalf.min() )
    hival = x2p[closest]

    if graphic:
        fwhm = hival - lowval
        figure(2)
        clf()
        plot( x, y/y.max(), 'bx-', label='original')
        plot( x1p, y1p, 'r-', x2p, y2p, 'r-', label='spline fit' )
        plot( (lowval, hival), (0.5, 0.5), 'k-', label='width' )
    
    return fwhm
 
