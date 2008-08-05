import numpy
from vectors import perpz, pt_nearest
from numpy.linalg import norm
from fit3d import vol2coords, fitlineNdw
from fit1d import fit_double_boltzman, double_boltzman_fwhm
from gaussian import fitgauss1d
from coordhandling import extract_line, extract_box
from neuron import pdbg
import neuron

neuron.debug = 'warning'

def kwparse(default_values, kwargs):
    '''edits a dictionary to provide default values
    if not otherwise supplied'''
    for k, v in default_values.iteritems():
        if not kwargs.has_key(k):
            kwargs[k] = v
    return kwargs


def bbox(pts):
    '''returns the upper and lower corners describing the
    bounding box that encompasses all points in ''pts'' '''
    return pts.min(0), pts.max(0)


def get_cnrs(pt, vec, fordsz, sidesz, mincnr, maxcnr):
    '''returns upper and lower (n-D) corners of a box defined
    by a central point ''pt'', a forward vector ''vec'' and
    width & length scaling values ''fordsz'' and ''sidesz'',
    and within the bounds of ''lclip'' and ''hclip'' '''
    perpnotz = numpy.cross(vec, perpz(vec))
    perpnotz /= norm(numpy.cross(vec, perpz(vec)))
    cnrs = numpy.array( bbox( pt + numpy.array([ fordsz * vec, \
                                           fordsz * -vec, \
                                           sidesz * perpz( vec ), \
                                           sidesz * -1 * perpz( vec ), \
                                           sidesz * perpnotz, \
                                           sidesz * -1 * perpnotz ]) ) )
    #cnrs[1] += 1
    cnrs = cnrs.clip( min=mincnr, max=maxcnr )
    return cnrs[0], cnrs[1]


def fwhm_diam(block, a, pD, curve=1.2, precision=0.005):
    '''Measures dendritic diameter from double-boltzman function fitted
    to in-plane line perpendicular to dendrite. Units are in pixels'''
    crossprof, profpts = extract_line( block.astype(int), a, pD )
    if crossprof.shape[0] < 3:
        return None, None
    pdbg('debug verbose', "dbg:", crossprof.shape)
    hpk = crossprof.max() / 2.
    xs = numpy.arange( len(crossprof) )

    c = numpy.where( crossprof == crossprof.max() )[0][0]
    w = 2
    A = crossprof.max()

    (c, A, w) = fit_double_boltzman( crossprof, p0=(c,A,w), \
                                     r=xs.astype(float), curve=curve )
    if c is None:
        return None, None
    else:
        interpxs = numpy.linspace( 0, xs.max(), xs.max() / precision + 1)
        dend_diam = double_boltzman_fwhm( c, A, w, curve, r=interpxs )
        # calcs fwhm, accounting for case
        # where peak of either one curve is not reached
        return dend_diam, hpk


def trace_one(img, a, D, di, hpk, **kwargs):
    '''performs one iteration of the tracing loop for
    tracing a structure through a confocal stack'''
    pdbg('debug verbose', )
    pdbg('debug terse', "start pt ", a, "\n     di", di)
    lclip, hclip = (0,0,0), img.shape
    pdbg('debug verbose', "clip vals %s,\n          %s" % \
          (lclip.__str__(), hclip.__str__()))
    kwargs = kwparse({'fitbox_length' : 8, \
                      'fitbox_width'  : di, \
                      'angthresh'     : numpy.pi / 5., \
                      'ofsthresh'     : 6, \
                      'dithresh'      : 4, \
                      'diam_width'    : di + 5, \
                      'hpkthresh'     : 1/2., \
                      'diam_cube_size': 2, \
                      'npt_sp'        : 8}, kwargs)
    ok = False
    while ((not ok) and (kwargs['npt_sp'] != 0)):
        pdbg('debug terse','->', kwargs['npt_sp'], '*')
        ok = True
        pt = numpy.round( a + kwargs['npt_sp'] * D )
        if neuron.debug == 'debug verbose':
            kwargs['frame'].plot_pt(pt, 'rD')
        pdbg('debug verbose', "predicting point", pt)
        lc, uc = get_cnrs(pt, D, kwargs['fitbox_length'], \
                        kwargs['fitbox_width'], lclip, hclip)
        pdbg('debug verbose',"cnrs", lc, ",\n      ", uc)
        bk = img[lc[0]:uc[0], lc[1]:uc[1], lc[2]:uc[2]]
        cds, wi = vol2coords(bk, 0)
        nDofs, nD = fitlineNdw(cds, wi)
        pdbg('debug verbose', "old D  %s\nnew D  %s" % (D, nD))
        # 
        if (numpy.inner(nD, D) < 0): # i.e. angle between vectors > pi/2.
            nD *= -1
            pdbg('debug verbose', "flipping D to %s" % (nD))

        # appropriateness checking
        angl = numpy.inner(nD, D)
        if angl < numpy.cos(kwargs['angthresh']):
            ok = False
            pdbg('warning',"ERROR: angle (%0.2f) too divergent!" % (angl))
            kwargs['npt_sp'] -= 1
            continue
        
        na = pt_nearest(pt, nDofs + numpy.array(lc), nD)
        pdbg('debug verbose', "candidate pt ", na)
        # check fitted point is close to predicted point
        if (abs(pt - na).max() > kwargs['ofsthresh']):
            ok = False
            pdbg('warning', "     ERROR: fitted point to far from predicted point!")
            kwargs['npt_sp'] -= 1
            continue

        perpnDz = numpy.cross( nD, perpz(nD) )
        dlc, duc = get_cnrs( na, nD, kwargs['diam_cube_size'], \
                             kwargs['diam_width'], lclip, hclip )
        pdbg('debug verbose',"diameter corners  %s,\n                  %s" \
              % (dlc.__str__(), duc.__str__()))
        diam_bk = img[dlc[0]:duc[0], dlc[1]:duc[1], dlc[2]:duc[2]]
        ndi, nhpk = fwhm_diam(diam_bk, na-dlc, perpz(nD))
        pdbg('debug verbose',"new diameter  ", ndi)
        # check new diameter is defined
        if ndi is None:
            ok = False
            pdbg('warning',"     ERROR: unable to calculate diameter!")
            kwargs['npt_sp'] -= 1
            continue
        # check new diameter is similar to old one
        if (di != 0) & (ndi < di / kwargs['dithresh']) \
            | (ndi > di * kwargs['dithresh']):
            ok = False
            pdbg('warning',"     ERROR: diameter (", di, ":", \
                 ndi, ") varied too quickly!")
            kwargs['npt_sp'] -= 1
            continue
        # check new peak is > some fraction of prev peak
        thr = hpk * kwargs['hpkthresh']
        if (nhpk < thr):
            ok = False
            pdbg('warning',"     ERROR: peak value too low! (%f < %f)" \
                 % (nhpk, thr))
            kwargs['npt_sp'] -= 1
            continue

    if (not ok) | (kwargs['npt_sp'] == -1):
        pdbg('debug terse', "Exiting")
        return None, None, None, None

    pdbg('debug terse', "     Okay")
    return na, nD, ndi, nhpk
