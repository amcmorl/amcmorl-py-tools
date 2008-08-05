from scipy.stats.stats import _chk_asarray
import numpy as np

def dparse(default_values, dic):
    '''edits a dictionary to provide default values
    if not otherwise supplied'''
    for k, v in default_values.iteritems():
        if not dic.has_key(k):
            dic[k] = v
    return dic


debug = 'error'
debug_levels = ['verbose','terse','warning','error']


def pdbg(level, *args):
    if debug_levels.index(level) >= debug_levels.index(debug) and args:
        if args[-1] == '*':
            args = list(args)
            args.pop()
            print " ".join([str(x) for x in args]), 
        else:
            print " ".join([str(x) for x in args])


def printif( condition, *args ):
    if condition:
        print " ".join([str(x) for x in args])

import os.path

def find_unique_filename(fname):
    i = 0
    if not os.path.exists(fname):
        print "Tested %s - not present." % fname
        return fname
    else:
        name, ext = os.path.splitext(fname)
        newfname = "%s_%d%s" % (name, i, ext)
        while os.path.exists(newfname):
            print "Tested %s - present" % (newfname)
            i += 1
            newfname = "%s_%d%s" % (name, i, ext)
        return "%s_%d" % (name, i)

def nanstd(x, axis=0, bias=False):
    """Compute the standard deviation over the given axis ignoring nans
 	
    :Parameters:
        x : ndarray
            input array
        axis : int
            axis along which the standard deviation is computed.
        bias : boolean
            If true, the biased (normalized by N) definition is used. If false,
            the unbiased is used (the default).
    
    :Results:
        s : float
            the standard deviation.
    Taken from scipy trac 5/9/07"""
    x, axis = _chk_asarray(x,axis)
    x = x.copy()
    Norig = x.shape[axis]
    
    Nnan = np.sum(np.isnan(x),axis)*1.0
    n = Norig - Nnan
    
    x[np.isnan(x)] = 0.
    m1 = np.sum(x,axis)/n
 	
    # Kludge to subtract m1 from the correct axis
    if axis!=0:
        shape = np.arange(x.ndim).tolist()
        shape.remove(axis)
        shape.insert(0,axis)
        x = x.transpose(tuple(shape))
        d = (x-m1)**2.0
        shape = tuple(np.array(shape).argsort())
        d = d.transpose(shape)
    else:
        d = (x-m1)**2.0
    m2 = np.sum(d,axis)-(m1*m1)*Nnan
    if bias:
        m2c = m2 / n
    else:
        m2c = m2 / (n - 1.)
    return np.sqrt(m2c)

def nansem(x, axis=0, bias=False):
    """Compute the standard error of the mean of array x along axis axis"""

    x, axis = _chk_asarray(x,axis)
    Norig = x.shape[axis]
    
    Nnan = np.sum(np.isnan(x),axis)*1.0
    n = Norig - Nnan

    stds = nanstd(x, axis, bias)
    return stds / np.sqrt(n)
    
