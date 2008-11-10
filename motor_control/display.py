import numpy as np
import motor_control
import matplotlib.pyplot as plt

# center-out -----------------------------------------------------------------
def setup_direction_axes(figure=1):
    target = motor_control.target
    dirnames = [','.join(["%2s" % y for y in x]) for x in target]
    n_dirs = target.shape[0]
    fig = plt.figure(figure)
    fig.clf()
    axes = {}
    lmarg, offset = 0.075, 0.1
    w, midgap = 0.3, 0.1
    outleft  = lmarg
    outright = 0.5 + midgap/2. + offset
    inleft   = lmarg + offset
    inright  = 0.5 + midgap/2.
    bmarg, hgap, tmarg = 0.05, 0.1, 0.05
    h = (1 - 3 * hgap - bmarg - tmarg)/4.
    outhigh  = bmarg + 3 * (h + hgap)
    inhigh   = bmarg + 2 * (h + hgap)
    inlow    = bmarg + 1 * (h + hgap)
    outlow   = bmarg
    for i in xrange(n_dirs):
        # set up axes
        # lbwh - look at x,y
        if np.all(target[i,0:2] == [-1,-1]):
            l = outleft
        elif np.all(target[i,0:2] == [-1,1]):
            l = inleft
        elif np.all(target[i,0:2] == [1,-1]):
            l = outright
        elif np.all(target[i,0:2] == [1,1]):
            l = inright
        # look at y,z
        if np.all(target[i,1:3] == [-1,-1]):
            b = outlow
        elif np.all(target[i,1:3] == [-1, 1]):
            b = outhigh
        elif np.all(target[i,1:3] == [1,-1]):
            b = inlow
        elif np.all(target[i,1:3] == [1, 1]):
            b = inhigh            
        axrect = [l, b, w, h]
        axname = dirnames[i]
        axes[axname] = fig.add_axes(axrect)
        axes[axname].set_title(axname)
        axes[axname].axvline(0, linewidth=1)
    return fig, axes

def plot_raster_series(spike_series, ax):
    '''plots a series of rasters (in a list of spike times lists)
    in the given axes. Y-values are intergral from 0 to no. of spike trains'''
    n_trains = len(spike_series)
    for i, spikes in enumerate(spike_series):
        plot_raster_row(spikes, ax, i, i+1)

def plot_raster_row(times, ax, ymin, ymax):
    # assumes no signals are >75 seconds prior to onset
    shp = times.shape + (1,)
    spc = np.zeros(shp, dtype = np.int8) - 100
    dblx = np.concatenate((times[:,np.newaxis],
                         times[:,np.newaxis], spc), axis=1).flatten()
    ymins = np.zeros(shp) + ymin
    ymaxs = np.zeros(shp) + ymax
    dbly = np.concatenate((ymins, ymaxs, spc), axis=1).flatten()
    mskx = np.ma.masked_less(dblx, -75)
    msky = np.ma.masked_less(dbly, -75)
    ax.plot(mskx, msky, 'k-', solid_capstyle='butt', linewidth=1)
    return dblx, dbly, mskx, msky
