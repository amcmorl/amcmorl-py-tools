from numpy import *
import numpy as n
from spine.model import spine_model, spine_image
from spine.utils import pxres, read_params, write_params, complete_params, \
     model_pxres
from string import strip
from rebin import congrid
import os

if os.uname()[4] == 'i686':
    from scipy import ndimage

def random_in_range( lims ):
    ''' returns a random number within the limits given by
    lims (a length-2 tuple), centered around the mid-point
    and with the tails clipped'''

    try:
        assert len(lims) == 2
        assert lims[0] < lims[1]
    except AssertionError:
        print "lims must have 2 values - (lower, upper)"
        
    midpoint = (lims[1] + lims[0]) / 2.
    half_range = (lims[1] - lims[0]) / 2.
    var = half_range**2 / 4 # gives (hopefully 2 s.d. to limiting values)
    # which should keep 95% of random results within appropriate range
    val = array( n.random.standard_normal() * var + midpoint )
    return val.clip( *lims )
    

def convert_coord(cd, sh1, sh2):
    '''correct point co-ords for convolve''s padding'''
    diff = abs( array( sh1 ) - array( sh2 ) )
    return cd + diff/2

    
def random_pars():
    '''creates a random set of spine parameters for use in random_spine etc

    split out of random spine - should stand alone'''
    
    ranges = { 'dend_diam'    : (0.45,  1.5),  \
               'neck_diam'    : (0.05, 0.75), \
               'head_diam'    : (0.05, 0.9), \
               'neck_length'  : (0.5,  1.5),  \
               'angle_roundz' : (0,    pi/3), \
               'angle_toz'    : (0,    pi/3), \
               'angle_den'    : (0,    pi/2) }
    noisiness = 3.
    pick = {}

    for k, v in ranges.iteritems():
        val = random_in_range( v )
        if k in ['angle_roundz','angle_toz', 'angle_den']:
            pick[k] = val
        else:
            # maximum accuracy is pixel spacing * 2 to meet Nyquists
            pick[k] = val - val % (model_pxres * 2)
        print k, pick[k]
        
    # head diam must be >= neck diam
    val = random_in_range( (pick['neck_diam'], 0.75) )
    pick['head_diam'] = val - val % (pxres() * 2)
    #print 'head_diam', pick['head_diam']
    return pick

    
def random_spine(fname, spname):
    '''creates a randomly generated spine named spname,
    and writes to file fname

    need to check if need convert_coord stuff still - don''t think I do
    since I use fftconvolve2 now, which preserves size'''
    
    pick = random_pars()

    (sp, orig, neck, head, skel) = \
         spine_model(px_res = pxres(), \
                     dend_diam = pick['dend_diam'], \
                     head_diam = pick['head_diam'], \
                     neck_diam = pick['neck_diam'], \
                     neck_length = pick['neck_length'], \
                     angle_roundz = pick['angle_roundz'], \
                     angle_toz = pick['angle_toz'])
    
    (spi, psf) = spine_image(sp, pxr = pxres())


    iorig = convert_coord( orig, sp.shape, spi.shape )
    ihead = convert_coord( head, sp.shape, spi.shape )
    ineck = convert_coord( neck, sp.shape, spi.shape )
    spi_hint = spi[tuple( ihead )]
    spi_nint = spi[tuple( ineck )]
    spi_orig = spi[tuple( iorig )]
    print "Head intensity", spi_hint / spi_orig
    print "Neck intensity", spi_nint / spi_orig

    pick['head_int'] = spi_hint / spi_orig
    pick['neck_int'] = spi_nint / spi_orig
    write_params( fname, spname, complete_params( pick ) )

    # need to change noise model to Possonian!!!
    noise = abs( randn(*(spi.shape)) * noisiness )
    spi *= ((255 - noise.max()) / spi.max())
    # scale spi to accomodate addition of noise
    spin = spi + noise

    #if os.uname()[4] != 'x86_64':
    rot_xy_angle = rand() * 90 # ndimage.rotate works in degrees !?!?
    spinner = ndimage.rotate( spin, rot_xy_angle, axes = (0, 1) )
    #else:
    #spinner = spin
    #   print "Warning: rotation requires ndimage," \
    # + "so is not supported on 64-bit systems."
    
    # resample to realistic pixel sizes
    startpxres = array( (pxres(), pxres(), pxres()) )
    # ^ defined first in spine_model
    despxres = array( (0.156, 0.156, 0.4) )
    # ^ would be acquired experimentally
    scl_factor = array(startpxres / despxres)
    newdims = tuple( (array( spinner.shape ) * \
                      scl_factor).round(0).astype(int) )
    spinners = congrid(spinner, newdims, method='linear' )

    shapestr = str(spinners.shape).replace(', ', '_').strip( ' ()' )
    spinners.tofile( spname + shapestr + '.pybin' )
    return (skel, sp, spi, spinners, psf)


def create_variants( filename, spname ):
    '''creates a family of parameters based on the given parameters,
    varying each parameter in turn by up to +/-50%'''

    par = read_params( filename, spname )

    dimvals = [-50, -40, -30, -20, -10, -5, 5, 10, 20, 30, 40, 50] 
    varatrs = {'angle_roundz':'ar', \
               'angle_toz'   :'at', \
               'angle_den'   :'ad', \
               'dend_diam'   :'dd', \
               'neck_length' :'nl', \
               'head_int'    :'hi', \
               'neck_int'    :'ni'}
    angles = ('ad', 'ar', 'at')
    sizes = ('dd', 'nl')
    intns = ('hi', 'ni')
    for (fulln, shortn) in varatrs.iteritems():
        if shortn in angles:
            for tval in range(1, len(dimvals) + 1):
                thesepars = par.copy()
                thesepars[fulln] = pi/3. * tval/float(len(dimvals))
                ncode = 'p' + str(tval)
                newspn = spname + '_' + ncode + '_' + shortn
                write_params( filename, newspn, thesepars )
        else:
            for tval in dimvals:
                thesepars = par.copy()
                thesepars[fulln] += thesepars[fulln] * tval / 100.0
                if shortn in sizes:
                    thesepars[fulln] -= thesepars[fulln] % pxres()
                ncode = str(tval).replace('-','m')
                if tval > 0:
                    ncode = 'p' + ncode
                newspn = spname + '_' + ncode + '_' + shortn
                write_params( filename, newspn, thesepars )
 

def recon_spine( filename, spname ):
    '''creates a spine model with the given parameters'''

    par = read_params( filename, spname )

    return spine_model( px_res = model_pxres, \
                        dend_diam = par['dend_diam'], \
                        head_diam = par['head_diam'], \
                        neck_diam = par['neck_diam'], \
                        neck_length = par['neck_length'], \
                        angle_roundz = par['angle_roundz'], \
                        angle_toz = par['angle_toz'], \
                        angle_den = par['angle_den'] )

