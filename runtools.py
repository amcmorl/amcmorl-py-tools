import os
import numpy as np

def rundisp(fn):
    def _(*args, **kwargs):
        try:
            res = fn(*args, **kwargs)
            os.system("notify-send 'Done!'")
            return res
        except:
            os.system("notify-send 'Oops!'")
            raise
    return _

def loras(filename, fn, *args, **kwargs):
    '''Load Or Run And Save'''
    if os.path.exists(filename):
        print "Loading from %s" % (filename)
        npzfile = np.load(filename)
        if len(npzfile.files) == 1:
            return npzfile[npzfile.files[0]]
        else:
            return tuple([npzfile[x] for x in npzfile.files])
    else:
        res = fn(*args, **kwargs)
        if type(res) == list:
            # deal with multiple returned values
            np.savez(filename, *res)
        else:
            # deal with a single returned value
            np.savez(filename, res)
        return res

def run_n_save(fn, filename):
    def _(*args, **kwargs):
        res = fn(*args, **kwargs)
        np.savez(filename, *res)
        return res
    return _

