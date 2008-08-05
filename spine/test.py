# -*- coding: utf-8 -*-

import numpy as n
import spine
import spine.fit, spine.utils, spine.model
from rebin import congrid
import time
from gaussian import gauss3d
from scipy.stats import poisson, stats
from vectors import unitvec
from glob import glob
from string import letters
from myutils import nansem
import pylab as p
import matplotlib as mpl

p.rcParams['text.usetex'] = False

# single test spine stuff here ------------------------------------------------

default_pars = { 'vol_size' : (222,201,160), \
              'd_diam'   : 1.0, \
              'space'    : 40, \
              'user_pts' : [(6,12,2),(20,14,3),(20,26,4)]
              }
default_pars['origin'] = [ default_pars['vol_size'][0] - \
                           default_pars['space'] - \
                           default_pars['d_diam'] / spine.utils.model_pxres, \
                           default_pars['vol_size'][1]/2., \
                           default_pars['vol_size'][2]/2. ]

test_pars = { 'origin_x'     : default_pars['origin'][0], \
              'origin_y'     : default_pars['origin'][1], \
              'origin_z'     : default_pars['origin'][2], \
              'head_diam'    : 0.5,  \
              'neck_diam'    : 0.25, \
              'neck_length'  : 0.75, \
              'angle_toz'    : 0.2, \
              'angle_roundz' : 0.2, \
              'angle_den'    : 0.2, \
              'dend_diam'    : default_pars['d_diam'], \
              'orient_den'   : 0.
              }

def get_vals(par, ns=None):
    if   par == 'hd' : vals = n.arange(0.05, 1.55, 0.05)
    elif par == 'hd+': vals = n.arange(0.7, 0.8, 0.005)
    elif par == 'nd' : vals = n.arange(0.05, 1.55, 0.05)
    elif par == 'nd+': vals = n.arange(0.7, 0.8, 0.005)
    elif par == 'nl' : vals = n.arange(0.1, 3.1, 0.1)
    elif par == 'dd' : vals = n.arange(0.1, 2.1, 0.1)
    elif par == 'arz':
        if ns == None: ns = 29
        vals = n.linspace(-2*n.pi/5., 2*n.pi/5., ns)
    elif par == 'arz+':
        if ns == None: ns = 49
        vals = n.linspace(-2*n.pi/5., 2*n.pi/5., ns)        
    elif par == 'atz':
        if ns == None: ns = 29
        vals = n.linspace(-2*n.pi/5., 2*n.pi/5., ns)
    elif par == 'A'  : raise NotImplemented("A is not currently supported.")
    elif par == 't'  : vals = n.arange(-2., 2., 0.1)
    return vals

def make_sp(model_pars='rand', other_pars=None, noise=0, fit_pars=True, \
            imagepars=None, name='anon', neck_length=None):
    '''Tests if spines can be reconstructed with variable levels of noise.

Creates a spine based on pars with a variable level of noise and tests whether
its parameters can be recovered after "imaging" i.e. blurring and resampling

:Parameters:
    model_pars : dict | "rand" | "test" (default = "rand")
        Spine parameters (as required by spine.model_weave)
        OR "rand" specifies that random spine parameters generated from EM
          measurements of CA1 spines should be used
        OR "test" signifies that standard test parameters should be used
    other_pars : dict
        Other parameters required for simulated imaging, and post-processing
        (vol_size, d_diam, space, user_pts)
    noise      : int (default = 0)
        Equivalent to maximum number of photons per pixel from spine
    fit_pars   : bool  | "test" (default = True)
        If TRUE, estimates parameters from image, otherwise uses original input
        parameters for reconstruction
        If "test", uses pre-supplied points - only works with test params and
          normal neck length
    imagepars  : dict
        A dictionary containing at least the items "fwhms" and "pxspc", which
        contain the FWHMs of the PSF and the pixel spacing of the experimental
        data in shape=(3,) ndarrays
    neck_length : float (default = None)
        optional specifier for neck length, _only_ used when
        model_pars = "test", to set neck length for geometry tests

:Returns:
    pars       : dictionary of parameters returned by reconstruction algorithm
    '''
    rand = False
    if imagepars == None:
        imagepars = { 'fwhms' : n.array((0.21,0.21,0.68)), \
                      'pxspc' : n.array((0.104,0.104,0.4)) }
    if other_pars   == None:
        other_pars = default_pars
    if model_pars   == 'test':
        model_pars = test_pars
        volsize, origin = None, None
        if type(neck_length) == float:
            model_pars['neck_length'] = neck_length
            volsize, origin = calc_size_and_origin(model_pars, \
                                                   spine.utils.model_pxres, \
                                                   space=other_pars['space'])
            model_pars['origin_x'] = float(origin[0])
            model_pars['origin_y'] = float(origin[1])
            model_pars['origin_z'] = float(origin[2])
    elif model_pars == 'rand':
        rand = True
        model_pars = create_parameters()
        volsize, origin = calc_size_and_origin(model_pars, \
                                               spine.utils.model_pxres, \
                                               space=other_pars['space'])
        model_pars['origin_x'] = float(origin[0])
        model_pars['origin_y'] = float(origin[1])
        model_pars['origin_z'] = float(origin[2])
    if volsize == None:
        if not other_pars.has_key('vol_size'):
            volsize, origin = calc_size_and_origin(model_pars, \
                                                   spine.utils.model_pxres)
        elif other_pars['vol_size'] <> n.nan and other_pars['vol_size'] <> None:
            volsize = other_pars['vol_size']
            print("Using default volume size and origin.")
            origin  = other_pars['origin']
        else:
            volsize, origin = calc_size_and_origin(model_pars, \
                                                   spine.utils.model_pxres)
    model = spine.model.build( \
         px_res   = spine.utils.model_pxres, \
         pars     = model_pars, \
         vol_size = volsize )
    mid_pxspc    = n.zeros(3) + spine.utils.model_pxres * 4.
    mid_fwhms    = mid_pxspc * 2.
    mid_fwhms_px = mid_fwhms / spine.utils.model_pxres
    blur1 = gauss3d(mid_fwhms_px[0], mid_fwhms_px[1], \
                    mid_fwhms_px[2], 1e-3)
    second_fwhms = n.sqrt(imagepars['fwhms']**2 - mid_fwhms**2) / mid_pxspc
    blur2 = gauss3d(second_fwhms[0], second_fwhms[1], \
                    second_fwhms[2], 1e-3)
    ires = spine.model.image_steps( \
                model, spine.utils.model_pxres, \
                mid_pxspc, imagepars['pxspc'], imagepars['fwhms'], \
                blur1 = blur1, blur2 = blur2 )
    del model
    ires -= ires.min()
    if noise <> 0:
        sz = ires.shape
        tmp = ires.flatten() / ires.max() * noise + 1e-12
        ires = poisson.rvs( tmp, size=ires.size ).reshape(sz).astype(float)
    sp = spine.Spine("%s%03d" % (name, noise), imagepars['pxspc'], \
                     imagepars['fwhms'], source_data = ires)
    if fit_pars == False:
        sp.params = model_pars
        sp.rotate_data()
    elif fit_pars == 'test':
        sp.load_user_pts(other_pars['user_pts'])
    else: # random parameters
        sp.params = model_pars
        print "Pick user points before proceeding."
    return sp

# residual function characterisation - execution ------------------------------

def test_pair(sp, xname, yname):
    self = sp
    # check params are valid
    if n.isnan(self.params['head_diam']):
        self.params['head_diam'] = 0.5
    if n.isnan(self.params['neck_diam']):
        self.params['neck_diam'] = 0.25
        
    # need to get dvec, dvec_a again...
    dvec_a, dvec = spine.fit.fit_dendrite_vector(self.rotated_data)
    print "Dendrite vector:", dvec, "\n & point:", dvec_a
    dvec *= self.pxspc / spine.utils.model_pxres
    dvec_a *= self.pxspc / spine.utils.model_pxres

    # pre-calculate some constant stuff for use in the loops
    model_size = (n.asarray(self.rotated_data.shape) * \
                  self.pxspc / spine.utils.model_pxres).round().astype(int)

    # fit angles and origin
    t0 = 0.
    arz0 = self.params['angle_roundz']
    atz0 = self.params['angle_toz']

    # fit diameters and neck length
    mid_pxspc = n.zeros(3) + spine.utils.model_pxres * 4
    mid_fwhms = mid_pxspc * 2
    mid_fwhms_px = mid_fwhms / spine.utils.model_pxres
    blur1 = gauss3d(mid_fwhms_px[0], mid_fwhms_px[1], \
                    mid_fwhms_px[2], 1e-3)
    second_fwhms = n.sqrt(self.fwhm**2 - mid_fwhms**2) / mid_pxspc
    blur2 = gauss3d(second_fwhms[0], second_fwhms[1], \
                    second_fwhms[2], 1e-3)

    A0 = self.rotated_data.sum()
    nd0 = self.params['neck_diam']
    hd0 = self.params['head_diam']
    nl0 = self.params['neck_length']
    dd0 = self.params['dend_diam']

    fout = open('/home/amcmorl/working/local_diam/resid_fn/save_%s_%s.txt' \
                % (xname, yname),'w')

    xs = get_vals(xname)
    ys = get_vals(yname)
    xname = xname.rstrip('+')
    yname = yname.rstrip('+')

    resids = n.zeros((xs.size,ys.size))
    for i in xrange(xs.size):      # arz
        for j in xrange(ys.size):  # dd
            exec('%s0 = xs[i]' % xname)
            exec('%s0 = ys[j]' % yname)
            p0 = (hd0, nd0, nl0, dd0, arz0, atz0, A0, t0)
            print "Processing %s=%6.3f, %s=%6.3f" % \
                  (xname, xs[i], yname, ys[j]),
            resids[i,j] = spine.fit.residual( p0, dvec_a, dvec, \
                                      model_size, blur1, blur2, \
                                      self, False )
            print "resid=%6.3f" % resids[i,j]
            fout.write('%s = %6.3f, %s = %6.3f, resid = %8.3f\n' % \
                       (xname, xs[i], yname, ys[j], resids[i,j]))
    fout.close()
    return xs, ys, resids

# residual function characterisation - visualisation --------------------------

def to_ind(val, mn, step):
    return int(n.round((val - mn) / step))

def collate_pair(fname, xplus=False, yplus=False):
    f = open(fname,'r')
    lines = f.readlines()
    xname, yname = lines[0].split()[0], lines[0].split()[3]
    if xplus:
        xname += '+'
    if yplus:
        yname += '+'
    xvals = spine.test.get_vals(xname)
    xmin, xstep = xvals[0], xvals[1] - xvals[0]
    yvals = spine.test.get_vals(yname)
    ymin, ystep = yvals[0], yvals[1] - yvals[0]
    dat = n.zeros((xvals.size, yvals.size))
    for line in lines:
        bits = line.split()
        xval = float(bits[2].rstrip(','))
        yval = float(bits[5].rstrip(','))
        resid = float(bits[8])
        xind = to_ind(xval, xmin, xstep)
        yind = to_ind(yval, ymin, ystep)
        dat[xind, yind] = resid
        if xind == 0:
            print line, "xind %d, yind %d\n" % (xind, yind)
        if yind == 0:
            print line, "xind %d, yind %d\n" % (xind, yind)
    return dat, xvals, yvals

def plot2D_surface_mesh(ar, xs, ys):
    from enthought.mayavi.tools import mlab
    from enthought.mayavi.sources.api import ArraySource
    from enthought.mayavi.modules.api import Axes, Surface
    from enthought.mayavi.filters.api import WarpScalar
    src = ArraySource()
    src.scalar_data = ar
    sc = mlab.figure()
    sc.add_child(src)
    surf = Surface()
    src.add_child(surf)
    surf.actor.mapper.interpolate_scalars_before_mapping = True
    warpsc = WarpScalar()
    src.add_child(warpsc)
    aspect = 4.
    height = (ar.max() - ar.min())
    width = xs.size
    cal = 0.2
    warpsc.filter.scale_factor = aspect * width / height * cal
    warpsurf = Surface()
    warpsc.add_child(warpsurf)
    warpsurf.actor.actor.property.representation = 'wireframe'
    warpsurf.actor.actor.property.line_width = 2.0
    ax = Axes()
    warpsc.add_child(ax)
    ax.axes.axis_label_text_property.bold = False
    ax.axes.axis_label_text_property.italic = False
    ax.axes.axis_label_text_property.shadow = True
    ax.axes.x_label = ""
    ax.axes.y_label = ""
    ax.axes.z_label = ""
    ax.axes.use_ranges = True
    ax.axes.ranges = [xs[0], xs[-1], ys[0], ys[-1], \
                      ar.min(), ar.max()]
    ax.property.line_width = 2.0

# snr relationship ============================================================

# generation of random spine parameters ---------------------------------------

def ca1_stats(verbose=False):
    # from Harris et al, 1989
    # dendritic diameters in microns, weighted by length in microns
    ddf = n.array((5.9,  9.7,  7.4,  8.5, 10.7, 11.2, 12.4))
    ddx = n.array((0.56, 0.47, 0.71, 0.65, 0.59, 0.41, 0.51))
    
    ddm = 1.8 * n.array(ddf * ddx).sum() / ddf.sum()
    #ddm = 0.7
    # ^ sum(f.Y)/sum(f)
    dds = n.sqrt(n.array(ddf * (ddx - ddm)**2).sum() / ddf.sum())
    # ^ sum(f.y**2)/sum(f)

    hvm = 0.051 # microns-cubed
    hdm = vol_to_diam(hvm)
    # hopefully sd_diam = 2 / pi / hdm**2 * hvs
    hvs = 0.07
    hds = 2 / n.pi / hdm**2 * hvs
    ndm, nds = 0.15, 0.06
    nlm, nls = 0.45, 0.29

    if verbose:
        print "%20s %0.3f ± %0.4f" % ('Dendritic diameter', ddm, dds)
        print "%20s %0.3f ± %0.4f" % ('Head diameter', hdm, hds)
        print "%20s %0.3f ± %0.4f" % ('Neck diameter', ndm, nds)
        print "%20s %0.3f ± %0.4f" % ('Neck length', nlm, nls)
        print "%20s %0.3f ± %0.4f" % ('Angle with dendrite', ndm, nds)

    stats = {'dend_diam'   : [ddm, dds], \
             'head_diam'   : [hdm, hds], \
             'neck_diam'   : [ndm, nds], \
             'neck_length' : [nlm, nls], \
             'angle_roundz': [-n.pi/3., n.pi/3.], \
             'angle_toz'   : [-n.pi/2., n.pi/2.], \
             'angle_den'   : [-n.pi/2., n.pi/2.]}
    return stats

def vol_to_diam(vol):
    '''head_vol = 4/3 * pi * (head_diam/2.)**3
    this does reverse calculation'''
    return (vol / (4/3.) / n.pi)**(1/3.) * 2

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
    val = n.array( n.random.standard_normal() * var + midpoint )
    return val.clip( *lims )

def create_parameters():
    ''' constructs n spine model variants with randomly varied parameters,
    based on mean parameters returned by mean_params'''
    stats = ca1_stats()
    norm_dists = ('head_diam','neck_diam','neck_length','dend_diam')
    uniform_dists = ('angle_toz', 'angle_den')
    limit_dists = ('angle_roundz')
    pars = spine.utils.complete_params( {} )
    for k in pars.keys():
        if k in norm_dists:
            tmp = n.random.standard_normal() * stats[k][1] \
                  + stats[k][0]
            #pars[k] = tmp - n.mod(tmp, spine.utils.model_pxres)
            pars[k] = tmp
        elif k in uniform_dists:
            pars[k] = n.random.random_sample() * \
                      (stats[k][0] - stats[k][1]) + (stats[k][1])
        elif k in limit_dists:
            pars[k] = random_in_range( tuple(stats[k]) )
        else:
            pass
    pars['orient_den'] = 0.
    return pars

def calc_size_and_origin(pars, pxres, space=25, dend_len=0):
    dpx  = pars['dend_diam']   / pxres
    nwpx = pars['neck_diam']   / pxres
    hpx  = pars['head_diam']   / pxres
    nlpx = pars['neck_length'] / pxres
    dlpx = dend_len            / pxres
    arz  = pars['angle_roundz']
    atz  = pars['angle_toz']
    dend_vec = n.array([0, n.cos(pars['angle_den']), n.cos(pars['angle_den'])])
    nec_vec_dir = n.array( (n.cos( arz ) * n.cos( atz ), \
                            # x component
                            n.sin( arz ) * n.cos( atz ), \
                            # y component
                            n.sin( atz )
                            # z component
                            ) )
    nec_in_den = (dpx/2.) / n.sqrt( 1 - n.dot(nec_vec_dir, dend_vec)**2 )
    nec_vec_len = nlpx + nec_in_den + hpx/2.
    nec_vec = nec_vec_len * nec_vec_dir
    vol_width = dpx / 2. + abs( nec_vec[0] ) + hpx / 2. + space * 2
    vol_height = n.array( (dpx, nwpx, hpx) ).max() \
                 + abs( nec_vec[2] ) + space * 2
    vol_depth = n.array( (dpx, nwpx, hpx, dlpx) ).max() \
                + abs( nec_vec[1] ) + space * 2
    sz = n.array( (vol_width, vol_depth, vol_height) ).round().astype(int)
    
    o = n.array( (vol_width - (space + dpx / 2.), \
                    vol_depth / 2. + nec_vec[1] / 2., \
                    vol_height / 2. + nec_vec[2] / 2.) \
                   ).round(0).astype(n.int)
    #print o
    return sz, o

# snr visualisation -----------------------------------------------------------

def plot_snr(pad=0.5):
    # assumes you mean from files in wdir
    # get list of files
    allfiles = glob('*.txt')
    noiselevels = [10,25,50,100,200,400]
    ncols = len(noiselevels)
    filear = []
    for noiselevel in noiselevels:
        filear.append([])
    # separate into list at each noise level
    nrows = 0
    for thisfile in allfiles:
        thisnoiselevel = int( \
            thisfile.lstrip(letters).rstrip(letters + '.').split('_')[0])
        col = noiselevels.index(thisnoiselevel)
        filear[col].append(thisfile)
        if len(filear[col]) > nrows: nrows = len(filear[col])
        
    # for each file at each noise level calculate absolute difference between
    #   construction and fitted head and neck diameters
    headdiffs = n.zeros((ncols, nrows)) + n.nan
    neckdiffs = n.zeros((ncols, nrows)) + n.nan

    for i, nl in enumerate(noiselevels):
        for j, thisfile in enumerate(filear[i]):
            print "%d,%d" % (i, j),
            headdiffs[i,j], neckdiffs[i,j], resid = read_optimisation_pars(thisfile)
    # for each noise level calculate mean and standard dev. head
    # and neck differences plot differences vs sqrt(noise level)
    snr = n.sqrt(n.asarray(noiselevels))
    headmeans, headsems = stats.nanmean(headdiffs, axis=1), \
                          nansem(headdiffs, axis=1)
    neckmeans, necksems = stats.nanmean(neckdiffs, axis=1), \
                          nansem(neckdiffs, axis=1)
    f = p.figure(1)
    f.clf()
    axhead = f.add_subplot(211)
    axhead.set_xlabel('SNR')
    axhead.set_ylabel('Head Relative Error')
    axhead.errorbar(snr, headmeans, headsems, label='Head diameter', color='b')
    axhead.set_xlim(snr.min() - pad, snr.max() + pad)
    axneck = f.add_subplot(212)
    axneck.set_xlabel('SNR')
    axneck.set_ylabel('Neck Relative Error')
    axneck.errorbar(snr, neckmeans, necksems, label='Neck diameter', color='g')
    axneck.set_xlim(snr.min() - pad, snr.max() + pad)
    #axhead.legend()
    
def read_optimisation_pars(fname):
    f = open(fname,'r')
    lines = f.readlines()
    if len(lines) >= 45:
        # for random-param files
        # asserts check that format is what we think it is
        assert(lines[3].split()[0] == 'head_diam')
        chead = n.abs(float(lines[3].split()[2]))
        assert(lines[13].split()[0] == 'neck_diam')
        cneck = n.abs(float(lines[13].split()[2]))
        assert(lines[34].split()[0] == 'head_diam')
        fhead = n.abs(float(lines[34].split()[2]))
        assert(lines[41].split()[0] == 'neck_diam')
        fneck = n.abs(float(lines[41].split()[2]))
        if len(lines) >= 46:
            if len(lines[45].split()) > 1:
                assert(lines[45].split()[0] == 'Residual')
                resid = float(lines[45].split()[2])
            else:
                resid = n.nan
        else:
            resid = n.nan
        headdiff = n.abs(chead - fhead) / chead
        neckdiff = n.abs(cneck - fneck) / cneck
        print "%s %5f %5f" % (fname, headdiff, neckdiff)
        return headdiff, neckdiff, resid
    elif len(lines) >= 43:
        # for test_param files
        # asserts check that format is what we think it is
        assert(lines[3].split()[0] == 'head_diam')
        chead = n.abs(float(lines[3].split()[2]))
        assert(lines[9].split()[0] == 'neck_diam')
        cneck = n.abs(float(lines[9].split()[2]))
        assert(lines[30].split()[0] == 'head_diam')
        fhead = n.abs(float(lines[30].split()[2]))
        assert(lines[37].split()[0] == 'neck_diam')
        fneck = n.abs(float(lines[37].split()[2]))
        assert(lines[42].split()[0] == 'Residual')
        resid = float(lines[42].split()[2])
        headdiff = n.abs(chead - fhead) / chead
        neckdiff = n.abs(cneck - fneck) / cneck
        print "%s %5f %5f" % (fname, headdiff, neckdiff)
        return headdiff, neckdiff, resid
    else:
        print "%s - -" % fname 
        return None, None, None

# geometry relationship =======================================================

def plot_geom(pad=0.5):
    allfiles = glob('*.txt')
    necklengths = [0.0, 1.5]
    ncols = len(necklengths)
    filear = []
    for necklength in necklengths:
        filear.append([])
    # separate into list at each noise level
    nrows = 0
    for thisfile in allfiles:
        thisnecklength = float( \
            thisfile.lstrip(letters).rstrip(letters + '.').split('_')[0]) / 10.
        try:
            col = necklengths.index(thisnecklength)
            filear[col].append(thisfile)
            if len(filear[col]) > nrows: nrows = len(filear[col])
        except ValueError:
            # if not listed, ignore file
            pass
        
    # for each file at each noise level calculate absolute difference between
    #   construction and fitted head and neck diameters
    headdiffs = n.zeros((ncols, nrows)) + n.nan
    neckdiffs = n.zeros((ncols, nrows)) + n.nan

    for i, nl in enumerate(necklengths):
        for j, thisfile in enumerate(filear[i]):
            print "%d,%d" % (i, j),
            headdiffs[i,j], neckdiffs[i,j], resid = read_optimisation_pars(thisfile)
    # for each noise level calculate mean and standard dev. head
    # and neck differences plot differences vs sqrt(noise level)
    snr = n.sqrt(n.asarray(necklengths))
    headmeans, headsems = stats.nanmean(headdiffs, axis=1), \
                          nansem(headdiffs, axis=1)
    neckmeans, necksems = stats.nanmean(neckdiffs, axis=1), \
                          nansem(neckdiffs, axis=1)
    f = p.figure(1)
    f.clf()
    axhead = f.add_subplot(111)
    axhead.set_xlabel('Neck length (micron)')
    axhead.set_ylabel('Head Relative Error')
    axhead.errorbar(snr, headmeans, headsems, label='Head diameter', color='b')
    axhead.set_xlim(snr.min() - pad, snr.max() + pad)
    
#     axneck = f.add_subplot(212)
#     axneck.set_xlabel('Neck length (micron)')
#     axneck.set_ylabel('Neck Relative Error')
#     axneck.errorbar(snr, neckmeans, necksems, label='Neck diameter', color='g')
#     axneck.set_xlim(snr.min() - pad, snr.max() + pad)
    #axhead.legend()

# residuals-accuracy relationship =============================================

def plot_residuals(pad=0.5):
    # assumes you mean from files in wdir
    # get list of files
    allfiles = glob('*.txt')
    noiselevels = [10,25,50,100,200,400]
    nnoises = len(noiselevels)
    filear = []
    for noiselevel in noiselevels:
        filear.append([])
    # separate into list at each noise level
    nrows = 0
    for thisfile in allfiles:
        thisnoiselevel = int( \
            thisfile.lstrip(letters).rstrip(letters + '.').split('_')[0])
        col = noiselevels.index(thisnoiselevel)
        filear[col].append(thisfile)
        if len(filear[col]) > nrows: nrows = len(filear[col])
        
    # for each file at each noise level calculate absolute difference between
    #   construction and fitted head and neck diameters
    headdiffs = n.zeros((nnoises, nrows)) + n.nan
    neckdiffs = n.zeros((nnoises, nrows)) + n.nan
    resids    = n.zeros((nnoises, nrows)) + n.nan

    for i, nl in enumerate(noiselevels):
        for j, thisfile in enumerate(filear[i]):
            print "%d,%d" % (i, j),
            headdiffs[i,j], neckdiffs[i,j], resids[i,j] \
                            = read_optimisation_pars(thisfile)

    colors = mpl.cm.Blues(n.arange(nnoises)/float(nnoises))
    f = p.figure(1, facecolor='w')
    f.clf()
    left = 0.075
    bottom = 0.075
    top = 0.05
    vgap = 0.075
    right = 0.025
    barwidth = 0.2
    hgap = vgap
    width = 1 - left - right - barwidth - hgap
    height = (1 - bottom - top - vgap)/2.
    axhead = f.add_axes([left, bottom, width, height])
    axhead.set_xlabel('Head Relative Error')
    axhead.set_ylabel('Residual Value')
    axneck = f.add_axes([left, bottom+height+vgap, width, height])
    axneck.set_xlabel('Neck Relative Error')
    axneck.set_ylabel('Residual Value')
    axbar = f.add_axes([left+width+hgap, bottom, barwidth, 1-top-bottom])
    axbar.imshow(colors[n.newaxis,...].transpose(1,0,2))
    # for each noise level plot head and neck accuracies vs residual
    # labeling points with snr
    plotlines = []
    labels = []
    for i, nl in enumerate(noiselevels):
        labels.append(str(n.sqrt(nl)))
        valids = n.equal(n.isnan(resids[i,:]), False)
        plotlines.append( axhead.plot(headdiffs[i,:][valids], \
                                      resids[i,:][valids], 'o', \
                                      color=colors[i]) )
        axneck.plot(neckdiffs[i,:][valids], resids[i,:][valids], 'o', \
                    color=colors[i])
    axhead.set_xlim(xmax = 1.01 * n.nanmax(headdiffs[n.equal(n.isnan(resids), \
                                                             False)]))
    axneck.set_xlim(xmax = 1.01 * n.nanmax(neckdiffs[n.equal(n.isnan(resids), \
                                                             False)]))
