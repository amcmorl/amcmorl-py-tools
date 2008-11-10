import numpy as np
import scipy #, scipy.io, scipy.signal
from multiregress import multiregress

target = np.array([[-1, -1, -1], \
                   [-1, -1,  1], \
                   [-1,  1, -1], \
                   [-1,  1,  1], \
                   [ 1, -1, -1], \
                   [ 1, -1,  1], \
                   [ 1,  1, -1], \
                   [ 1,  1,  1]])

def unbiased_histogram(spikes, bins=None, range=None, new=None):
    '''a drop in replacement for np.histogram that implements the unbiased
    firing rate metric employed by Zhanwu and Rob Kass'''
    vers = np.version.version
    if vers == '1.0.4':
        hist, bin_times = np.histogram(spikes,
                                       bins=bins,
                                       range=[range[0], range[1]])
        bin_times = np.hstack((bin_times, bin_times[-1] - bin_times[-2]))
    else:
        hist, bin_times = np.histogram(spikes,
                                       bins=bins,
                                       range=[range[0], range[1]],
                                       new=new)
    hist = hist.astype(float)
    hist += 1
    for i in xrange(hist.size):
        prespikes = spikes[spikes < bin_times[i]]
        if prespikes.size == 0:
            prespike = range[0]
        else:
            prespike = prespikes[-1]

        postspikes = spikes[spikes > bin_times[i+1]]
        if postspikes.size == 0:
            postspike = range[1]
        else:
            postspike = postspikes[0]
        time_interval = postspike - prespike
        hist[i] /= time_interval
    return hist, bin_times

def calc_pd(rates, with_err=False, n_samp=1000):
    '''calculates preferred directions of a unit from a center out task
    
    could be a method of a "unit" class
    which references the relevant data_file items

    Calculates the vector components of the preferred direction:
    c_x, c_y, and c_z from the
    b_x, b_y, and b_z partial regression co-efficients described in eqn 1 of
    Schwartz et al, 1988:

        d(M) = b + b_x.m_x + b_y.m_y + b_z.m_z,

    and where

        c_x = b_x / sqrt(b_x^2 + b_y^2 + b_z^2)
        c_y = b_y / sqrt(b_x^2 + b_y^2 + b_z^2)
        c_z = b_z / sqrt(b_x^2 + b_y^2 + b_z^2)

    The algorithm regresses the firing rates (spike count / duration) against
    [1, cos theta_x, cos theta_y, cos theta_z], to give [b, b_x, b_y, b_z]

    Parameters
    ----------
    rates : array, shape (ndir,)
        1d array of firing rates in each of the (usually eight) directions
    with_err : bool, default False
        whether to calculate and return std errors of pd components
    n_samp : integer
        when calculating error, is the number of samples from the
        theoretical distribution of PD to calculate
    Returns
    -------
    pd : array, shape (3,)
        preferred direction, unit-length
    err : array, shape (3, )
        std errors of components of pds
    '''
    cosInput = target / np.sqrt( ( target[0]**2 ).sum() )    
    if not with_err:
        n_directions = target.shape[0]
        xin = np.hstack((np.ones((n_directions, 1)), cosInput))
        b = np.linalg.lstsq(xin, rates)[0][1:]
        pd = b.copy()
        k = np.sqrt(np.sum(pd**2))
        # normalize to unit length
        pd /= k
        return pd, k    
    else:
        b, b_s = multiregress(cosInput, rates) # partial regression coeffs
        b = b[1:] # remove constant term for consideration of direction
        k = np.sqrt(np.sum(b**2))
        pd = b / k
        # b.shape (3,)
        
        # Now generate n_samp samples from normal populations
        # (mean: b_k, sd: err_k). Will have shape (n_samp, 3).
        b_rnd = np.random.standard_normal(size=(n_samp,3))
        # se.shape (3,)
        b_rnd *= b_s
        b_rnd += b
        ks = np.sqrt(np.sum(b_rnd**2, axis=1))
        #pds = b_rnd / ks[...,np.newaxis]
        #pd_s = np.std(pds, axis=0, ddof=1)
        k_s = np.std(ks, ddof=1)
        return pd, b_s, k, k_s

def calc_pp(rates, positions):
    '''Calculates positional gradient for a cell, based on the nine firing
    rates that correspond to eight corner target positions and the center
    position, supplied in the positions variable

    The regression is of the equation:

    fr = b_0 + b_x * x + b_y * y + b_z * z

    and b_0, c_x, c_y, and cz are returned

    where c_x,c_y,c_z = b_x,b_y,b_z / k; k = sqrt(b_x**2 + b_y**2 + b_z**2)

    Parameters
    ----------
    rates : array_like

    positions : array_like

    Returns
    -------
    pp : array_like, shape (3,)
      vector of positional gradient
    k : scalar
      modulation depth of positional encoding
    b0 : scalar
      baseline firing rate (b0 in regression
    '''
    n_rates = rates.size
    #assert False
    all_pos = positions
    #all_pos = np.vstack((np.zeros(3), target))
    xin = np.hstack((np.ones((n_rates, 1)), all_pos))
    res = np.linalg.lstsq(xin, rates)[0]
    b0, pp = res[0], res[1:]
    k = np.sqrt(np.sum(pp**2))
    pp /= k
    return pp, k, b0

def get_target_idx(vector):
    '''Return index of row in motor_control.util.target
    with same direction as vector.'''
    trg_dir = target / np.apply_along_axis(np.linalg.norm, 1,
                                                target)[..., np.newaxis]
    vec_dir = vector / np.linalg.norm(vector)    
    return np.nonzero( np.abs( trg_dir - vec_dir).sum(axis=1)
                          < 1e-8)[0].item()

# def get_target_num(targetPos, targetArr):
#     ''' returns index within targetArr which equals targetPos
#     Deprecated - use get_target_idx.'''
#     signedTargetPos = np.sign(targetPos)
#     return np.nonzero((targetArr == signedTargetPos).min(1) == True)[0][0]

# used ?????????????????????????????????????????????????????????????????????

# def getOldTargetDir(dircode):
#     '''returns target direction and index of that direction in new target
#     array'''
#     if dircode == 1:
#         dirn = (1, -1, -1)
#     elif dircode == 2:
#         dirn = (1, -1, 1)
#     elif dircode == 3:
#         dirn = (1, 1, -1)
#     elif dircode == 4:
#         dirn = (1, 1, 1)
#     elif dircode == 5:
#         dirn = (-1, -1, -1)
#     elif dircode == 6:
#         dirn = (-1, -1, 1)
#     elif dircode == 7:
#         dirn = (-1, 1, -1)
#     elif dircode == 8:
#         dirn = (-1, 1, 1)
#     tview = target.view([('a', int), ('b', int), ('c', int)])
#     dirnar = np.array(dirn, dtype=([('a', int), ('b', int), ('c', int)]))
#     return dirn, np.nonzero(tview == dirnar)[0][0]

# def norm_correlate(a, b, mode='valid'):
#     a_ = (a - a.mean()) / scipy.stats.var(a)
#     b_ = (b - b.mean()) / scipy.stats.var(b)
#     return np.correlate(a_, b_, mode=mode)
    
# def getTaskType(eccentricity, illusion):
#     '''return a code for the task type,
#     based on the eccentricity and illusion values'''
#     if (eccentricity != 1): # ellipse - visual
#         if (illusion == 1):
#             return 0 # ellipse
#         else:
#             return 1 # ellipse-illusion
#     else:
#         if (illusion == 1):
#             return 2 # circle
#         else:
#             return 3 # circle-illusion

# def getTaskName(eccentricity, illusion):
#     '''return a code for the task type,
#     based on the eccentricity and illusion values'''
#     if (eccentricity != 1): # ellipse - visual
#         if (illusion == 1):
#             return "ellipse"
#         else:
#             return "ellipse-illusion"
#     else:
#         if (illusion == 1):
#             return "circle"
#         else:
#             return "circle-illusion"

# def nanPad(param):
#     '''pad last element of last dimension with nans'''
#     sh = list(param.shape)
#     sh[-1] += 1
#     newparam = np.empty(sh, dtype=param.dtype)
#     newparam[...,0:-1] = param
#     newparam[...,-1] = np.nan
#     return newparam

# def calcLag(pd, spike_hist, hand_vel, times, mode='first'):
#     """Calculation of firing unit lag from an arbitrary time window.

#     Calculates the lag between predicted firing rate of a unit based on its
#     preferred direction and hand velocity and its actual firing rate.

#     Parameters
#     ----------
#     pd   : array_like
#         3-D vector of xyz co-ordinates of preferred direction of unit,
#         as measured from center-out task.
#     spike_hist : array_like
#         Histogram of spike occurences, binned according to bin_edge_times
#     hand_vel : array_like
#         Hand velocity (xyz) vectors at time periods defined by bin_edge_times
#     bin_edge_times : array_like
#         n+1 time bins (including right-hand edge) into which spike_hists and
#         hand_velocity have been placed

#     Returns
#     -------
#     tau : float
#         Lag time in s between measured and predicted firing rates
#     """
#     # calculate predicted firing rate, given pd and hand velocity
#     # both this and actual firing rate are in bins given by bin_edge_times
#     pfr = pd[0] + np.dot(pd[np.newaxis, 1:], hand_vel).squeeze()
#     afr = spike_hist
#     # replace negative predicted firing rates with 0
#     pfr[pfr < 0] = 0.
    
#     # replace nan values with geometric mean of other values
#     pfr[np.isnan(pfr)] = scipy.stats.gmean(pfr[~np.isnan(pfr)])

#     #x_cor = scipy.signal.correlate(pfr, afr, mode='valid')
#     x_cor = norm_correlate(afr, pfr, mode='valid')
#     lag_bin = int((x_cor.size - 1) / 2.) - x_cor.argmax()

#     # if all bins are the same size, then next bit is easy
#     if mode == 'first':
#         lag_time = lag_bin * (times[1] - times[0])
#     # but
#     elif mode == 'mean':
#         mean_time = np.diff(times).mean()
#         lag_time = lag_bin * mean_time
#     else:
#         print "Time conversion mode is not supported.\n" \
#               "Return value is bin number."
#         return lag_bin

#     return lag_time, lag_bin

