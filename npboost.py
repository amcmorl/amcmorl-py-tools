import numpy as n
from scipy import weave

def indices3d(vol_size):
    '''Returns the equivalent of n.indices for 3-D only.
    Gives a ~6x speed up.
    '''
    xsz = vol_size[0]
    ysz = vol_size[1]
    zsz = vol_size[2]
    assert( len(vol_size) == 3 )
    assert( type(xsz) == int )
    assert( type(ysz) == int )
    assert( type(zsz) == int )
    code = '''
    #line 14
    int i,j,k;
    npy_int dims[4] = {3,xsz,ysz,zsz};
    long *curpos;
    PyArrayObject* ar =
      (PyArrayObject *)PyArray_ZEROS(4, &dims[0], NPY_LONG, 0);
    for (i = 0; i<xsz; i++)
        for (j = 0; j<ysz; j++)
            for (k = 0; k<zsz; k++)
            {
                curpos = (long *)PyArray_GETPTR4(ar, 0, i, j, k);
                *curpos = i;
                curpos = (long *)PyArray_GETPTR4(ar, 1, i, j, k);
                *curpos = j;
                curpos = (long *)PyArray_GETPTR4(ar, 2, i, j, k);
                *curpos = k;
                
            }
    return_val = PyArray_Return(ar);
    '''
    return weave.inline( code, ['xsz','ysz','zsz'] ) 

