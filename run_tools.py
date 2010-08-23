import os
import numpy as np
import os.path

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
