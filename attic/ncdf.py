from numpy import asarray, float64, ndarray, array, uint8
import numpy as n
import os.path

posmods = ['pynetcdf', 'Scientific.IO.NetCDF', 'Scientific_NetCDF',
           'pupynere', 'scipy.io.netcdf']
tried = 0
for mod in posmods:
    try:
        exec("from %s import NetCDFFile" % mod)
        break
    except ImportError:
        tried += 1
if tried == len(posmods):
    raise ImportError("No module containing NetCDFFile found.")

def r(fname, varname='var0',compat=False, \
          start=None, size=None, smallest=True):
    '''Reads in a netCDF file into a numpy array.

    Usage:
      data = rncdf(fname, varname=''var0'', compat=False,
        start=None, size=None)
        
    where
      fname = filename
      varname = variable name
      compat = whether to switch axes around to support
        deconvolution output format
      start = co-ordinates of lower corner
      stop = size of block to extract'''
    if start <> None:
        start = n.asarray(start)
    if size <> None:
        size  = n.asarray(size)
    if not os.path.isfile(fname):
        print fname, "is not a valid file name"
    else:
        ncfile = NetCDFFile( fname, 'r' )
        var = ncfile.variables[varname]
        if not start is None and not size is None:
            if start.size != size.size:
                raise ValueError("start & size must have same no. of elements.")
            if compat:
                start = start[::-1]
                size = size[::-1]
            stop = start+size
            index_string = ",".join(['start[%s]:stop[%s]' % \
                                     (i,i) for i in range(len(start))])
            exec('data = var[' + index_string + ']')
        else:
            data = n.asarray(var[:])

        if data.dtype.kind == 'i' and smallest:
            data = n.cast[n.uint8](data)
        if compat:
            return data.transpose(2,1,0)
        else:
            return data


def w(fname, data, name='var0'):
    ncfile = NetCDFFile( fname, 'w' )
    dimslist = []
    for i in range( len( data.shape ) ):
        dimname = 'dim%s' % i
        ncfile.createDimension( dimname, data.shape[i] )
        dimslist.append(dimname)
    if data.dtype.type == uint8:
        typecode = 'i'
    else:
        typecode = 'd'
    var = ncfile.createVariable( name, typecode, tuple( dimslist ) )
    var[:] = data
    ncfile.close()


def wncdf(fname, data, name='var0'):
    print "wncdf is now deprecated - use w instead"
    w(fname, data, name=name)


def rncdf(fname, varname='var0', compat=False, \
          start=None, size=None, smallest=True):
    print "rncdf is now deprecated - use r instead"
    return r(fname, varname=varname, compat=compat, \
      start=start, size=size, smallest=smallest)

def get_size(fname, varname='var0', compat=False):
    ncfile = NetCDFFile(fname, 'r')
    var = ncfile.variables[varname]
    if compat == False:
        tx = var[:,0,0]
        ty = var[0,:,0]
        tz = var[0,0,:]
    else:
        tx = var[0,0,:]
        ty = var[0,:,0]
        tz = var[:,0,0]
    ncfile.close()
    return n.array((tx.shape[0], ty.shape[0], tz.shape[0]))
