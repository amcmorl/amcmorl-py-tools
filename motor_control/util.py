import numpy as np
import scipy, scipy.stats #, scipy.io, scipy.signal

target = np.array([[-1, -1, -1], \
                   [-1, -1,  1], \
                   [-1,  1, -1], \
                   [-1,  1,  1], \
                   [ 1, -1, -1], \
                   [ 1, -1,  1], \
                   [ 1,  1, -1], \
                   [ 1,  1,  1]])

def getOldTargetDir(dircode):
    '''returns target direction and index of that direction in new target
    array'''
    if dircode == 1:
        dirn = (1, -1, -1)
    elif dircode == 2:
        dirn = (1, -1, 1)
    elif dircode == 3:
        dirn = (1, 1, -1)
    elif dircode == 4:
        dirn = (1, 1, 1)
    elif dircode == 5:
        dirn = (-1, -1, -1)
    elif dircode == 6:
        dirn = (-1, -1, 1)
    elif dircode == 7:
        dirn = (-1, 1, -1)
    elif dircode == 8:
        dirn = (-1, 1, 1)
    tview = target.view([('a', int), ('b', int), ('c', int)])
    dirnar = np.array(dirn, dtype=([('a', int), ('b', int), ('c', int)]))
    return dirn, np.nonzero(tview == dirnar)[0][0]

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

def norm_correlate(a, b, mode='valid'):
    a_ = (a - a.mean()) / scipy.stats.var(a)
    b_ = (b - b.mean()) / scipy.stats.var(b)
    return np.correlate(a_, b_, mode=mode)

def get_target_num(targetPos, targetArr):
    ''' returns index within targetArr which equals targetPos'''
    signedTargetPos = np.sign(targetPos)
    return np.nonzero((targetArr == signedTargetPos).min(1) == True)[0][0]
    
def getTaskType(eccentricity, illusion):
    '''return a code for the task type,
    based on the eccentricity and illusion values'''
    if (eccentricity != 1): # ellipse - visual
        if (illusion == 1):
            return 0 # ellipse
        else:
            return 1 # ellipse-illusion
    else:
        if (illusion == 1):
            return 2 # circle
        else:
            return 3 # circle-illusion

def getTaskName(eccentricity, illusion):
    '''return a code for the task type,
    based on the eccentricity and illusion values'''
    if (eccentricity != 1): # ellipse - visual
        if (illusion == 1):
            return "ellipse"
        else:
            return "ellipse-illusion"
    else:
        if (illusion == 1):
            return "circle"
        else:
            return "circle-illusion"

def nanPad(param):
    '''pad last element of last dimension with nans'''
    sh = list(param.shape)
    sh[-1] += 1
    newparam = np.empty(sh, dtype=param.dtype)
    newparam[...,0:-1] = param
    newparam[...,-1] = np.nan
    return newparam

def calcLag(pd, spike_hist, hand_vel, times, mode='first'):
    """Calculation of firing unit lag from an arbitrary time window.

    Calculates the lag between predicted firing rate of a unit based on its
    preferred direction and hand velocity and its actual firing rate.

    Parameters
    ----------
    pd   : array_like
        3-D vector of xyz co-ordinates of preferred direction of unit,
        as measured from center-out task.
    spike_hist : array_like
        Histogram of spike occurences, binned according to bin_edge_times
    hand_vel : array_like
        Hand velocity (xyz) vectors at time periods defined by bin_edge_times
    bin_edge_times : array_like
        n+1 time bins (including right-hand edge) into which spike_hists and
        hand_velocity have been placed

    Returns
    -------
    tau : float
        Lag time in s between measured and predicted firing rates
    """
    # calculate predicted firing rate, given pd and hand velocity
    # both this and actual firing rate are in bins given by bin_edge_times
    pfr = pd[0] + np.dot(pd[np.newaxis, 1:], hand_vel).squeeze()
    afr = spike_hist
    # replace negative predicted firing rates with 0
    pfr[pfr < 0] = 0.
    
    # replace nan values with geometric mean of other values
    pfr[np.isnan(pfr)] = scipy.stats.gmean(pfr[~np.isnan(pfr)])

    #x_cor = scipy.signal.correlate(pfr, afr, mode='valid')
    x_cor = norm_correlate(afr, pfr, mode='valid')
    lag_bin = int((x_cor.size - 1) / 2.) - x_cor.argmax()

    # if all bins are the same size, then next bit is easy
    if mode == 'first':
        lag_time = lag_bin * (times[1] - times[0])
    # but
    elif mode == 'mean':
        mean_time = np.diff(times).mean()
        lag_time = lag_bin * mean_time
    else:
        print "Time conversion mode is not supported.\n" \
              "Return value is bin number."
        return lag_bin

    return lag_time, lag_bin

def calc_pd(rates):
    '''calculates preferred directions of a unit from a center out task
    
    could be a method of a "unit" class
    which references the relevant data_file items

    Calculates the b_x, b_y, and b_z co-efficients described in eqn 1 of
    Schwartz et al, 1988.

    d(M) = b + b_x.m_x + b_y.m_y + b_z.m_z

    The algorithm regresses the firing rates (spike count / duration) against
    [1, cos theta_x, cos theta_y, cos theta_z], to give [b, b_x, b_y, b_z]
    '''
    n_directions = target.shape[0]
    cosInput = target / np.sqrt( ( target[0]**2 ).sum() )    
    xin = np.hstack((np.ones((n_directions, 1)), cosInput))
    coeffs = np.linalg.lstsq(xin, rates)[0]
    pd = np.array(coeffs)
    return pd
