import numpy as n
from string import strip
import os.path

# ---------------------------------------------------------
# service routines to complete, write and read spine
# parameter data
# ---------------------------------------------------------


model_pxres = 0.015
model_pxres_ar = n.array((0.015,0.015,0.015))
real_pxres = n.array((0.156,0.156,0.4))


def complete_params( p ):
    '''fills in any gaps in params p with nan'''
    shudve = ['orig_file',   \
              'coords',      \
              'dend_diam',   \
              'angle_roundz',\
              'angle_toz',   \
              'angle_den',   \
              'neck_diam',   \
              'neck_length', \
              'head_diam',   \
              'head_int',    \
              'neck_int']
    for item in shudve:
        if not p.has_key( item ):
            p[item] = n.nan
    return p

def write_params( fname, spname, p ):
    '''appends a line containing spine parameters p to file fname'''
    if not os.path.exists(fname):
        setup_params_file(fname)
    f = open( fname, 'a' )
    f.write( "%16s|" % spname )
    f.write( "%10s|" % str( p['orig_file'] ) )
    f.write( "%10s|" % str( p['coords'] ) )
    f.write( "%6.4f|" % p['dend_diam'] )
    f.write( "%6.4f|" % p['angle_roundz'] )
    f.write( "%6.4f|" % p['angle_toz'] )
    f.write( "%6.4f|" % p['angle_den'] )
    f.write( "%6.4f|" % p['neck_diam'] )
    f.write( "%6.4f|" % p['neck_length'] )
    f.write( "%6.4f|" % p['head_diam'] )
    f.write( "%6.4f|" % p['head_int'] )
    f.write( "%6.4f\n" % p['neck_int'] )
    f.close()
    return None


def read_params( fname, spname ):
    '''reads spine spname from file fname, returning params
    Use spname == \'*\' to return all spines, in which case
    the function returns a list of spine param dictionaries'''
    f = open( fname )
    lines = f.readlines()
    ps = None
    if spname == '*':
        ps = []
    for line in lines:
        if line[0] != '#' and line.strip() != '':
            # ignore commented & empty lines
            bits = line.split('|')
            if spname == '*':
                lineps = {'spname'       : strip( bits[ 0] ), \
                          'orig_file'    : strip( bits[ 1] ), \
                          'coords'       : strip( bits[ 2] ), \
                          'dend_diam'    : float( bits[ 3] ), \
                          'angle_roundz' : float( bits[ 4] ), \
                          'angle_toz'    : float( bits[ 5] ), \
                          'angle_den'    : float( bits[ 6] ), \
                          'neck_diam'    : float( bits[ 7] ), \
                          'neck_length'  : float( bits[ 8] ), \
                          'head_diam'    : float( bits[ 9] ), \
                          'head_int'     : float( bits[10] ), \
                          'neck_int'     : float( bits[11] )}
                ps.append(lineps)
            elif spname == strip( bits[0] ):
                ps = {'spname'       : strip( bits[ 0] ), \
                      'orig_file'    : strip( bits[ 1] ), \
                      'coords'       : strip( bits[ 2] ), \
                      'dend_diam'    : float( bits[ 3] ), \
                      'angle_roundz' : float( bits[ 4] ), \
                      'angle_toz'    : float( bits[ 5] ), \
                      'angle_den'    : float( bits[ 6] ), \
                      'neck_diam'    : float( bits[ 7] ), \
                      'neck_length'  : float( bits[ 8] ), \
                      'head_diam'    : float( bits[ 9] ), \
                      'head_int'     : float( bits[10] ), \
                      'neck_int'     : float( bits[11] )}
    return ps


def read_clist(fname, spname):
    ''' Reads user parameters titled name from file '''
    f = open( fname )
    lines = f.readlines()
    ps = None
    if spname == '*':
        ps = []
    for line in lines:
        if line[0] != '#' and line.strip() != '':
            # ignore commented & empty lines
            bits = line.split()
            if spname == '*':
                lineps = [n.array(bits[1:4]).astype(int), \
                          n.array(bits[7:10]).astype(int), \
                          n.array(bits[10:14]).astype(int)]
                ps.append(lineps)
            elif spname == strip( bits[0] ):
                ps = [n.array(bits[1:4]).astype(int), \
                      n.array(bits[7:10]).astype(int), \
                      n.array(bits[10:14]).astype(int)]
    return ps


def pxres():
    '''returns xyz pixel resolution of models (in microns)'''
    return 0.015


def setup_params_file(fname):
    if os.path.isfile(fname):
        raise ValueError("File already exists.")
    f = open(fname, 'w')
    f.write( "#" )
    f.write( "%15s|" % 'Spine name' )
    f.write( "%10s|" % 'Orig. file' )
    f.write( "%10s|" % 'Co-ords' )
    f.write( "%6s|" % 'Dend d' )
    f.write( "%6s|" % 'Ang rz' )
    f.write( "%6s|" % 'Ang tz' )
    f.write( "%6s|" % 'Ang dz' )
    f.write( "%6s|" % 'Neck d' )
    f.write( "%6s|" % 'Neck l' )
    f.write( "%6s|" % 'Head d' )
    f.write( "%6s|" % 'Head I' )
    f.write( "%6s\n" % 'Neck I' )
    f.write( "#" )
    f.write( "%s\n" % ''.join(['-' for x in xrange(99)]) )
    f.close()
    return None
    

def test_params():
    return {'origin_x': 133,
            'origin_y': 100,
            'origin_z': 80,
            'angle_den': 0.0,
            'angle_roundz': 0.0,
            'angle_toz': 0.0,
            'dend_diam': 1.0,
            'head_diam': 0.5,
            'neck_diam': 0.25,
            'neck_length': 0.75}

# # ---------------------------------------------------------
# # read spine data from file, fit model to it, append results
# # to data file and read file to gather multi-spine data
# # ---------------------------------------------------------

# def plot_data( params ):

#     dd = array( data['dend_diams'] )
#     nd = array( data['neck_diams'] )
#     hd = array( data['head_diams'] )
#     nl = array( data['neck_lengths'] )

#     # dendrite diameter vs neck diameter
#     p1 = subplot(221)
#     p1.get_xaxis().label.set_text(r'$\rm{dendrite\ diameter (\mu m )}$')
#     p1.get_yaxis().label.set_text(r'$\rm{neck\ diameter (\mu m)}$')
#     p1.plot( dd, nd, 'b+' )

#     # dendrite diameter vs head diameter
#     p2 = subplot(222)
#     p2.get_xaxis().label.set_text(r'$\rm{dendrite\ diameter (\mu m )}$')
#     p2.get_yaxis().label.set_text(r'$\rm{head\ diameter (\mu m)}$')
#     p2.plot( dd, hd, 'b+' )
    
#     p3 = subplot(223)
#     p3.get_xaxis().label.set_text(r'$\rm{neck\ diameter (\mu m )}$')
#     p3.get_yaxis().label.set_text(r'$\rm{head\ diameter (\mu m)}$')
#     p3.plot( nd, hd, 'b+' )
    
#     p4 = subplot(224)
#     p4.get_xaxis().label.set_text(r'$\rm{neck\ length (\mu m )}$')
#     p4.get_yaxis().label.set_text(r'$\rm{neck\ diameter (\mu m)}$')
#     p4.plot( nl, nd, 'b+' )


def init_spine_view(sp):
    from enthought.mayavi.sources.array_source import ArraySource
    from enthought.mayavi.modules.iso_surface import IsoSurface
    from enthought.mayavi.modules.outline import Outline
    from enthought.mayavi.tools import mlab
    mv = mlab.figure()
    src = ArraySource()
    src.scalar_data = sp
    mv.add_source(src)
    o = Outline()
    isos = IsoSurface()
    mv.add_module(isos)
    mv.add_module(o)
    return src


def later_spine_view(src, sp):
    src.scalar_data = sp
