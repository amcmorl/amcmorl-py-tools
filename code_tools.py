def dparse(default_values, dic):
    '''edits a dictionary to provide default values
    if not otherwise supplied'''
    for k, v in default_values.iteritems():
        if not dic.has_key(k):
            dic[k] = v
    return dic

import numpy as np
anynans = lambda x : np.any(np.isnan(x))

# -----------------------------------------------------------------------------
# interactive tools
# -----------------------------------------------------------------------------

def explore_npz(file_name):
    '''
    Print a useful summary of npz files.
    '''
    f = np.load(file_name)
    for k, v in f.iteritems():
        print k, type(v), v.shape, v[0:2]
        print "\n"
