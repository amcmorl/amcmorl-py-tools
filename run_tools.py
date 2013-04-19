import os
import sys
import numpy as np
import os.path

def rundisp(fn):
    def _(*args, **kwargs):
        try:
            res = fn(*args, **kwargs)
            os.system("notify-send 'Done!'")
            return res
        except:
            typ, val, tbk = sys.exc_info()
            os.system("notify-send '%s: %s'" % (typ.__name__, str(val)))
            raise
    return _

def _chk_ar(inp, rank=None, shape=None):
    inp = np.asarray(inp)
    assert np.rank(inp) == rank
    return inp

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

def shc(s, c):
    '''
    SHell Colour: add a colour `c` specification to string `s`.

    s : string
      string to colourize
    c : string
      character code of colour to use, ala matplotlib,
      currently one of 'r', 'g'
    '''
    col = {'r' : 31,
           'g' : 32}
    if c not in col.keys():
        raise ValueError("`c` must be one of %s" % (str(col.keys())))
    return "\033[1;%2dm%s\033[1;m" % (col[c], s)


#!/usr/bin/env python
# encoding: utf-8

import sys

def switch_on_call_tracing(base):
    def trace_calls(frame, event, arg):
        if event != 'call':
            return
        co = frame.f_code
        func_name = co.co_name
        if func_name == 'write':
            # Ignore write() calls from print statements
            return
        func_line_no = frame.f_lineno
        func_filename = co.co_filename
        if not base in func_filename:
            return
        caller = frame.f_back
        caller_line_no = caller.f_lineno
        caller_filename = caller.f_code.co_filename
        print '%20s:%2s of %30s from %30s:%2s' % \
            (func_name, func_line_no, func_filename,
             caller_filename, caller_line_no)
        return

    sys.settrace(trace_calls)
