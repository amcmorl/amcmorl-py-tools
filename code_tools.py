def dparse(default_values, dic):
    '''edits a dictionary to provide default values
    if not otherwise supplied'''
    for k, v in default_values.iteritems():
        if not dic.has_key(k):
            dic[k] = v
    return dic

import numpy as np
anynans = lambda x : np.any(np.isnan(x))

