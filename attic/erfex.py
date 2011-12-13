import numpy as np
import matplotlib.pyplot as plt
import plot_tools
from types import FunctionType

def _sample_erf(par, which, erf, *arg, **kwarg):
    '''
    Sample error function around a given set of fit parameters.

    Parameters
    ----------
    par : ndarray
      fit parameters
    erf : callable
      error function, must have signature (p, data, options)
    args : list
      arguments passed to `erf`
    nval : int
      number of values to evaluate in each cross-section

    Returns
    -------
    vals : ndarray
      independent values put into erf
    err : ndarray
      error function values returned by erf
    '''
    if 'nval' in kwarg:
        nval = kwarg['nval']
    else:
        nval = 20
    vals = np.array([np.linspace(x - x/2., x + x/2., nval) for x in par[which]])
    err = np.zeros_like(vals)
    for ipar in which:
        pcop = par.copy()
        for ival, val in enumerate(vals[ipar]):
            pcop[ipar] = val
            err[ipar, ival] = erf(pcop, *arg)
    return vals, err

class ErFEx(object):
    '''
    ErFEx - the Error Function Explorer

    For a given set of data, parameters and an error function,
    plot the orthogonal cross-sections through parameter space, centered
    at the value of the parameters.

    Parameters
    ----------
    erf : callable
      error function
    par : array_like
      fit parameters
    erf_args : list
      list of other arguments to `erf`
    '''
    def __init__(self, erf, par, erf_args, **kwarg):
        self.par = np.asarray(par)
        assert(type(erf) == FunctionType)
        self.erf = erf
        self.erf_args = erf_args

        # active subset of par
        if 'which' in kwarg.keys():
            self.which = np.asarray(kwarg['which'])
        else:
            # use all
            self.which = np.arange(self.pars.size)

        # make plot
        self.fig = plt.figure()
        if 'margins' in kwarg.keys():
            m = kwarg['margins']
        else:
            m = plot_tools.Margins(0.1, 0.1, 0.05, 0.1, 0.1, 0.075)
        self.grid = plot_tools.create_plot_grid(self.which.size, ncols=3,
                                                margin=m, fig=self.fig,
                                                direction='row', sharex='none',
                                                xspines='all', yspines='all')

        # labels for axes - default to 'p0', 'p1'...
        if 'label' in kwarg.keys():
            label = kwarg['label']
            if len(label) != self.which.size:
                raise(ValueError("number of labels must match number"
                                 "of selected parameters."))
        else:
            label = ['p%d' % (x) for x in self.which]

        # evaluate xsections
            vals, err = _sample_erf(self.par, self.which, self.erf, *self.erf_args)

        # plot xsections
        for i, g in enumerate(self.grid):
            g.set_title(label[i])
            g.plot(vals[i], err[i])
            g.axvline(self.par[i], color='r')
            yt = self.grid[i].get_ylim()
            self.grid[i].set_yticks([yt[0], yt[-1]])
            
        self.fig.canvas.mpl_connect('button_press_event', self.on_click)

    def on_click(self, evt):
        # if clicked in one of the grid axes
        if evt.inaxes in self.grid:
            pidx = self.grid.index(evt.inaxes)
            self.par[which[pidx]] = evt.xdata
            print "Changing" 
            print "Current error value: %f" % (self.get_current_erf())
            self.plot()
            
    def plot(self):
        val, err = _sample_erf(self.par, self.which, self.erf, *self.erf_args)
        for i, g in enumerate(self.grid):
            line = g.lines[0]
            line.set_xdata(val[i])
            line.set_ydata(err[i])
            g.set_ylim(np.min(err[i]), np.max(err[i]))
            g.set_yticks([np.min(err[i]), np.max(err[i])])
            g.set_xlim([np.min(val[i]), np.max(val[i])])
            line = g.lines[1]
            xval = self.par[i]
            line.set_xdata([xval, xval])
        plt.draw()

    def get_current_erf(self):
        return self.erf(self.par, *self.erf_args)
