import rebin
import os
import numpy as n
import spine.fit
from numpy import array, arange, pi, where, \
     matrix, linspace, fromfile, float64
from coordhandling import extract_box, extract_line
from vectors import perpz, angle_between, unitvec, rotate_about_centre
from gaussian import gauss1d, fitgauss1d
from fit1d import fit_double_boltzman, double_boltzman_by_centre, \
     double_boltzman_fwhm
from numpy.linalg import norm
from os import path

graphic = True
try:
    import pylab as p
except:
    graphic = False
    
from primitives import cylinder
from fit3d import vol2coords, fitlineNdw
from fit2d import midcrossings
from spine.utils import write_params, complete_params
import ncdf
from scipy.interpolate import splprep, splev
from myutils import printif
# requires data file binarydiam_vs_blurreddiam.nc
#          in '/home/amcmorl/working/forward_model

def sp_extract_dend_diam( den, pdvec, sp, res, psf_fwhm, \
                          precision=0.005, verbose=False ):
    ''' measure dend_diam
    - extract in-plane line perpendicular to dendrite (=profile)'''
    (crossprof, profpts) = extract_line( sp, den, pdvec )
    xs = arange( len(crossprof) ) * res[0]
    #dend_diam = midcrossings( xs, crossprof, k = 3, \
    # thresh = 0.01 * crossprof.max())
    if graphic and verbose:
        p.figure(1)
        p.clf()
        p.subplot(211)
        p.plot(xs, crossprof)
    (A, w, c) = fitgauss1d( crossprof, r=xs )
    printif( verbose, "Awc:", A, w, c )
    #A = crossprof.max()
    #w = crossprof
    A, w, c = abs(A), abs(w), abs(c)
    # ^ corrects for weirdness in fitting routine
    # which sometimes returns -ve number
    curve = psf_fwhm/2.5 # arbitrary guess - seems to fit okay
    (c, A, w) = fit_double_boltzman( crossprof, p0=(c,A,w), \
                                     r=xs, curve=curve )
    interpxs = n.linspace( 0, xs.max(), xs.max() / precision + 1)

    dend_diam = double_boltzman_fwhm( c, A, w, curve, r=interpxs )
    #profit = gauss1d( xs, A, dend_diam, c )
    profit = double_boltzman_by_centre(c, A, w, curve=curve, \
                                       r=interpxs )
    if graphic and verbose:
        p.plot(interpxs, profit)
    printif( verbose, "pre-correction", dend_diam )
    if graphic and verbose:
        p.subplot(212)
    dend_diam = correct_diam_for_blurring( dend_diam )
    printif( verbose, "dendrite width =", dend_diam, "microns" )
    return dend_diam

def sp_extract_neck_length( spvec, dvec, hd, dd, verbose=True ):
    '''calculate neck_length
     - defined as from edge of dendrite to edge of head'''
    printif(verbose, "Spine", spvec, "Dendrite", dvec)
    pttopt = norm(spvec)
    printif(verbose, "Total length", pttopt)
    hd_rad = hd/2.
    printif(verbose, "Head radius", hd_rad)
    in_dend = dd/(2. * n.sqrt(1-(n.dot(unitvec(dvec), unitvec(spvec)))**2))
    printif(verbose, "in dendrite", in_dend)
            
    nl = pttopt - hd_rad - in_dend
    printif( verbose, "neck length =", nl, "microns" )
    return n.maximum(nl, 0.)

def extract_orient_den(dvec, spvec, pxspc, verbose=True):
    # calculates angle relative to +ve y axis
    # - subject to flip depending on orientation of spine
    orient_den = n.arctan( dvec[0] / float(dvec[1]) )
    printif(verbose, "Angle to y (for round z):", orient_den)
    spvecr = rotate_about_centre(spvec[0:2], (0,0), -orient_den)
    if spvecr[0] > 0:
        orient_den = (orient_den + n.pi) % (2*n.pi)
    return orient_den
        
def sp_extract_angle_roundz( spvec, pdvec, orient_den, verbose=False):
    '''calculate angle_roundz
    - angle between perp-to-dendrite (= pdvec) and spine

    - defined as +ve when spvec y component is -ve,
    when dvec is oriented along y axis (i.e. pdvec is oriented along x axis)
    and spine is in +x direction
    '''
    printif(verbose, "Spine vector", spvec)
    printif(verbose, "Dendrite normal0", pdvec)
    printif(verbose, "Angle to y (for round z):", orient_den)
    spvecr = rotate_about_centre(spvec[0:2], (0,0), orient_den)
    printif(verbose, "Rotated spine vector", spvecr)
    spvec_xy  = n.array( spvec, copy=True )
    spvec_xy[2] = 0.
    if n.sign(pdvec[0]) <> n.sign(spvec_xy[0]):
        pdvec *= -1
    printif(verbose, "Dendrite normal", pdvec)
    arz = angle_between( spvec_xy, pdvec )
    if n.isnan(arz):
        # arbitrarily set it to zero - shouldn't make any difference
        arz = 0

    if n.sign(spvecr[1]) <> n.sign(spvecr[0]):
        arz *= -1
    printif(verbose, "angle round z", arz, "rads")
    return arz

def sp_extract_angle_toz( spvec, verbose=False ):
    '''calculate angle_toz
    angle between spvec and z
    (need to multiply by res since model
       will be done with xy = z res)'''
    printif( verbose, "To z:", spvec )
    atz = pi/2 - angle_between( spvec, (0.,0.,-1.) )
    printif( verbose, "angle to z", atz, "rads" )
    return atz

def sp_extract_dend_angle( dvecr, spvecr, verbose=False ):
    '''calculate angle of dendrite to z:
    angle = tan^{-1}(z/sqrt(x^2 + y^2))

    whether the angle is +ve or -ve depends on up/down of dendrite in +y,
    when spine is pointing in -x direction

    1) rotate spvec in x-y so that d_vec lies along y axis (one way or other)
    2) if spvec x is > 0, rotate d_vec by 180 deg in x-y
    3) angle_den is angle between x-y plane and dendrite, in +y direction
    4) also have to take into account anisotropy in pixel spacing x-y vs z
    '''
    x, y, z = list(dvecr)
    printif( verbose, "x,y,z:", x, y, z )
    printif( verbose, "spine vector", spvecr )
    #if y == 0:
    #    dtoy = 0
    #else:
    #    dtoy = n.arctan(x/y) # angle to +y axis
    #printif(verbose, "Angle to y (dend. z):", dtoy)
    #spvecr = rotate_about_centre( spvec[0:2], (0,0), -dtoy )
    #printif(verbose, "Rotated spine vec", spvecr)
    #dvecr = rotate_about_centre( (x,y), (0,0), -dtoy )
    #printif(verbose, "Rotated dendrite vec", dvecr)
    angle_den = n.arctan(z/n.sqrt(x**2 + y**2))
    if dvecr[1] < 0:
        angle_den *= -1
    if spvecr[0] > 0:
        angle_den *= -1
    printif( verbose, "Z angle of dendrite:", angle_den )
    return angle_den

def correct_diam_for_blurring(blurred_diam, verbose=False):
    filedir = path.split(__file__)[0]
    filename = 'binarydiam_vs_blurreddiam.nc'
    relfile = filedir + '/' + filename
    if not os.path.isfile(relfile):
        raise RuntimeError('Binary-to-blurred diameter' + \
                           'correction data file is not present.')
    bin_blr = ncdf.r(relfile)

    bins = bin_blr[0]
    blrs = bin_blr[1]

    k = 5
    s = 1.
    nest = -1
    interpvals = linspace( 0, 1, 251 )
    # 251 simply gives enough points to give one point close to 0.5
    tckp, u = splprep( [bins, blrs], s=s, k=k, nest=nest )
    bins_f, blrs_f = splev( interpvals, tckp )

    dif_to_val = abs(blurred_diam - blrs_f) 
    closest_pt = where( dif_to_val == dif_to_val.min() )

    err = array(dif_to_val[closest_pt] / blurred_diam * 100.)[0]
    printif( verbose, "fitting error = %0.2f%s" % (err, '%') )
    bin_val = float( bins_f[closest_pt] )

    if graphic and verbose:
        p.plot(bins, blrs, 'bo', label='model')
        p.plot(bins_f, blrs_f, 'r-', label='fit')
        p.axvline(bins_f[closest_pt], color='k')
        p.axhline(blrs_f[closest_pt], color='k')
        p.xlim(0,2)
        p.ylim(0,2)
    
    return bin_val

def line_intersect(pa1, pa2, pb1, pb2):
    '''return the intersection of two 3-D vectors

    (int, lamz) = line_intersect(pa1, pa2, pb2, pb2)

    inputs:

      px1, px2 = two points (3-D vectors) on line x

      NB: the z-component of line a direction vector
          is ignored and the appropriate value
          returned to allow intersection,

    outputs:

      int = intersection co-ordinates

      lamz = z-component of line a'''
    
    Da = pa1 - pa2
    Db = pb1 - pb2
    a = matrix( (Da[0:2], -Db[0:2]) )
    b = matrix( (pb1[0:2] - pa1[0:2]) )
    (m,n) = tuple(array(b * a.I).flatten())
    lamz = (pa1[2] - pb1[2]) / (Da[2] * m + Db[2] * n)
    intscn = pa1 + Da * m
    
    return intscn, lamz

#def build_binarydiam_vs_blurreddiam():

# CODE GRAVEYARD ---------------------------------------------------------------

#     scale_factor = array( (viewscl, viewscl, 1) )
#     # useful scaling for viewing in pyvis  
#     zoomsz = boxsz * scale_factor
#     zoomed = rebin.rebin_neighbour( tozoom, tuple( zoomsz ) )
#     pv.SetImg( zoomed, 'spine' )
#     pv.ClearClickCoords()
#     # get coords of spine head, neck and dendrite (2 pts)
#     print "click (in this order):\n" + \
#           " 1- spine head\n" + \
#           " 2- spine neck centre\n" + \
#           " 3- spine origin\n" + \
#           " 4- another point on dendrite\n"
#     return tozoom, sofs

# def pull_spine_coords(pv, s, res, inname, sofs):
#     '''pulls spine co-ordinates from pv,
#     converts into parameters for spine_model,
#     and writes into a file

#     needs pv, small volume, res,
#     filename of original file and spine co-ordinates'''

#     # extract click_coords and convert
#     # to match scaling in display
#     coords = pv.GetClickCoords()
#     sclfac = array( (10,10,1) )
#             # defined by scaling for spine display in 
#             # extract_spine and extract_random_spine
#     if len(coords) < 4:
#         print "[process_spine] Not enough co-ords entered."
#         return None
#     else:
#         fname = raw_input('Select file to store data: ')
#         spname = raw_input('Give this spine a name: ')
        
#         p = { 'orig_file' : inname, 'coords' : str(sofs) }
#         # extract points from data
#         den = array( coords[-1] ) / sclfac
#         da = array( coords[-2] ) / sclfac
#         spb = array( coords[-3] ) / sclfac
#         spa = array( coords[-4] ) / sclfac
#         (ori, zslope) = line_intersect(spa, spb, da, den)
#         denv = s[den[0], den[1], den[2]]
#         print "dendrite 2", tuple( den ), "value", denv
#         oriv = s[ori[0], ori[1], ori[2]]
#         print "origin", tuple( ori ), "value", oriv

#         sp_vec = spb - spa
#         dvec = den - da
#         (spprof, spprofpts) = extract_line( s, spa, sp_vec )
#         xs = arange( len(spprof) ) * res[0]
#         clf()
#         plot( xs, spprof )

#         neck = copy(ori)
#         head = copy(ori)
#         neckv = s[neck[0], neck[1], neck[2]]
#         print "neck", tuple( neck ), "value", neckv
#         headv = s[head[0], head[1], head[2]]
#         print "head", tuple( head ), "value", headv

#         # measure dend_diam
#         # extract in-plane line perpendicular to dendrite (=profile)
#         #dvec = den - ori
#         pd_vec = perpz( dvec )
#         print "s shape", s.shape
#         print "den", den
#         print "pd_vec", pd_vec
#         (crossprof, profpts) = extract_line( s, den, pd_vec )
#         xs = arange( len(crossprof) ) * res[0]
#         plot(xs, crossprof)
#         (A, dend_diam, c) = fitgauss1d( crossprof, r=xs )
#         dend_diam = abs(dend_diam)
#           # sometimes fitting gets confused with sign
#         profit = gauss1d( xs, A, dend_diam, c )
#         plot(xs, profit)
#         print "dendrite width =", dend_diam, "microns"
#         p['dend_diam'] = dend_diam
             
#         # calculate neck_length
#         # (defined as from centre of dendrite)
#         sp_vec = head - ori
#         nl = norm( sp_vec * res )
#         print "neck length =", nl, "microns"
#         p['neck_length'] = nl

#         # calculate angle_roundz
#         # define cylinder based on head and origen points
#         vsz = s.shape
#         mask_cyl = cylinder(vsz, head, (head - ori), 0.75/res.max())
#         print "mask_cyl shape", mask_cyl.shape
#         masked = s * mask_cyl
#         print "max masked", masked.max()
#           # try 1 micron width for now
#           # arbitrary threshold of 10 (should cut out noise)
#         (s_mask_coords, s_mask_weights) = vol2coords( s * mask_cyl, 10 )
#         print "s mask coords", s_mask_coords
#         (a, sp_vec) = fitlineNdw(s_mask_coords, s_mask_weights)
#         print "sp_vec", sp_vec, "as", type(sp_vec)
#         sp_vec *= res # correct for anisotropic pixel spacing
#         print "fitted vector angle", sp_vec
        
#         # angle between perp-to-dendrite (= pd_vec) and spine
#         sp_vec_xy  = array( sp_vec )
#         sp_vec_xy[2] = 0.
#         arz = angle_between( sp_vec_xy, pd_vec )
#         # convert to positive angle between 0 and pi/2
#         if arz < 0:
#             arz *= -1
#         if arz > pi / 2.:
#             arz = pi - arz
#         print "angle round z", arz, "rads"
#         p['angle_roundz'] = arz

#         # calculate angle_toz
#         # angle between sp_vec and z
#         # (need to multiply by res since model
#         #   will be done with xy = z res)
#         atz = pi/2 - angle_between( sp_vec * res, (0,0,1) )
#         if atz < 0:
#             atz *= -1
#         print "angle to z", atz, "rads"
#         p['angle_toz'] = atz

#         # calculate normalized intensity values
#         neck_int = neckv / oriv.astype(float)
#         p['neck_int'] = neck_int
#         head_int = headv / oriv.astype(float)
#         p['head_int'] = head_int
        
#         # write results to file
#         write_params( fname, spname, complete_params( p ) )




# def extract_random_spine2(pv, fname=None, dims=None, s=None):
#     '''reads a randomly generated spine from a file created by file
#     (or from a pre-existing array), displays the spine in pyvis and
#     presents co-ordinate gathering instructions for user

#     Usage: (s, pxres, fname, sofs) =
#        extract_random_spine( pv, fname=fname, dims )

#     or: (s, pxres, fname, sofs) = extract_random_spine( pv, s=s )'''
    
#     if s == None:
#         s = fromfile(fname, dtype=float64).reshape(dims)
#     sz = array( s.shape )
#     newsz = sz * array( (10, 10, 1) )
#     news = rebin.rebin_neighbour(s, tuple( newsz ) )
#     pv.SetImg(news, 'dummy spine')
    
#     pv.ClearClickCoords()
#     # get coords of spine head, neck and dendrite (2 pts)
#     print "click (in this order):\n" + \
#           " 1- spine line - end\n" + \
#           " 2- spine line - other end\n" + \
#           " 3- dendrite line - end\n" + \
#           " 4- another point on dendrite\n" + \
#           "and then run \'t = pf.process_spine_coords(pv, s, res, fname, sofs)\'"
#     if fname == None:
#         fname = 'randsp'

#     sofs = array( (0,0,0) )
#     pxres = array( (0.156, 0.156, 0.4) )
#     return [s, pxres, fname, sofs]


# def get_spine_points(pv=None, sp=None, fname=None, dims=None, \
#                        pxres=(0.156, 0.156, 0.4)):
#     '''reads a randomly generated spine from a file created by file
#     (or from a pre-existing array), displays the spine in pyvis and
#     presents co-ordinate gathering instructions for user

#     Usage:
#       lpars = get_spine_points( pv, fname=<datafilename>, \
#                                   dims=<dims of datafile> )
#       calc_spine_params( lpars, store_fname, spname )'''

#     if pv == None:
# 		pv = pyvis()
    
#     if sp == None:
#         fformat = fname.split('.')[1]
#         if fformat == 'nc':
#             from ncdf import rncdf
#             sp = rncdf(fname)
#         elif fformat == 'pybin':
#             sp = fromfile(fname, dtype=float64).reshape(dims)
#         else:
#             raise ValueError("File format unknown.")
    
#     sz = array( sp.shape )
#     newsz = sz * array( (viewscl, viewscl, 1) )
#     news = rebin.rebin_neighbour(sp, tuple( newsz ) )
#     pv.SetImg(news, fname)

#     # get coords of spine head, neck and 2 pts of dendrite
#     pv.ClearClickCoords()
#     print "click (in this order):\n" + \
#           " 1  head\n" + \
#           " 2  neck centre\n" + \
#           " 3  spine origin\n" + \
#           " 4  dendrite\n"
#     if fname == None:
#         fname = 'sp'
#     sofs = array( (0,0,0) )
#     return [pv, sp, pxres, fname, sofs]


#viewscl = 10   # scaling factor for spine display
#viewbigscl = 2 # scaling factor for big dataset display
