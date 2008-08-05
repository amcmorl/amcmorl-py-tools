debug = 'debug verbose'
debug_levels = ['debug verbose','debug terse','warning','error']

def pdbg(level, *args):
    if debug_levels.index(level) >= debug_levels.index(debug) and args:
        if args[-1] == '*':
            args = list(args)
            args.pop()
            print " ".join([str(x) for x in args]),
        else:
            print " ".join([str(x) for x in args])
