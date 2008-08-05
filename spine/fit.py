import numpy as n, pylab as pl
import spine, spine.model
from spine.utils import model_pxres
from fit3d import fitlineNdw, vol2coords
from binary import biggest
from scipy.ndimage import median_filter
import time
from gaussian import fwhm2k

try:
    from enthought.mayavi.tools import mlab
    from enthought.mayavi.sources.array_source import ArraySource
    from enthought.mayavi.modules.iso_surface import IsoSurface
    from enthought.mayavi.modules.outline import Outline
    display_mode = '3-D'
except ImportError:
    display_mode = '2-D'

def fit_dendrite_vector(real, threshold=0.75):
    real=n.asarray(real, dtype=float)
    # get threshold from histogram
    footprint = n.ones((3,3,1),dtype=float)
    
    #footprint[1] = True
    #footprint[:,1] = True

    filt = median_filter(real, footprint=footprint)
    filt /= filt.max()
    
    # threshold dendrite
    mask = biggest(filt > threshold)
    
    cds, wi = vol2coords(mask * real)
    a, D = fitlineNdw(cds, wi)

    return a, D

def create_fit_weights(sp, fwhm=2.5):
    '''calculate an array of co-efficients by which to weight source data when
    fitting dendrite. Should favour pixels located close to user_pt origin.

    fwhm in the full-width at half-maximum, in micron, of the resulting kernel
    '''
    ks = n.array([fwhm2k(i) for i in list(fwhm/sp.pxspc)])
    centre = sp.user_pts[1] # get origin pt
    dimvals = n.abs(n.indices( sp.source_data.shape ) - \
              centre[...,n.newaxis, n.newaxis, n.newaxis]) \
          / ks[...,n.newaxis, n.newaxis, n.newaxis]
    Is = n.exp(-(dimvals**2).sum(0))
    return Is / Is.max()
    

def residual(p, d_a, d_D, vol_size, \
             blur1, blur2, sp, graphic=False, scene=None, iteration=None):
    '''Returns sum of squares of difference between model data and source
    data, aligned at origin. Takes smallest common region for comparison.
    
    d_a - in micron
    d_D - unitlength
    '''
    d_a = d_a.copy()
    d_a /= model_pxres
    head_diam, neck_diam, neck_length, dend_diam, \
               angle_roundz, angle_toz, A, t = p

    # apply non-negativity penalty (will prevent more than one test with
    # negative values)
    penalty = 0
    if (head_diam < 0) or (neck_diam < 0) or \
       (neck_length < 0) or (dend_diam < 0) or (A < 0):
        penalty = 1000

    # and by making everything positive, this means negative
    # values won't prevent continuation of the routine
    head_diam   = n.abs(head_diam)
    neck_diam   = n.abs(neck_diam)
    if (head_diam < neck_diam):
        penalty = 1000
    neck_length = n.abs(neck_length)
    dend_diam   = n.abs(dend_diam)    
    A = n.abs(A)
    
    print ""
    o = d_a + t * d_D # must be in model pixel co-ordinates
    #angle_den = n.arctan( d_D[2] / n.sqrt( d_D[0]**2 + d_D[1]**2 ) )
    angle_den = sp.params['angle_den']
    #print "Residual angle_den: %0.3f" % angle_den
    print "%s\n hd=%5.2f nd=%5.2f nl=%5.2f arz=%5.2f atz=%5.2f" % \
          (time.asctime(), head_diam, neck_diam, neck_length, \
           angle_roundz, angle_toz)

    model = spine.model.build( \
        px_res = spine.utils.model_pxres, \
        pars = {'origin_x'    : o[0], \
                'origin_y'    : o[1], \
                'origin_z'    : o[2], \
                'angle_den'   : angle_den, \
                'dend_diam'   : dend_diam, \
                'head_diam'   : head_diam, \
                'neck_diam'   : neck_diam, \
                'neck_length' : neck_length, \
                'angle_roundz': angle_roundz, \
                'angle_toz'   : angle_toz
                }, \
        vol_size = vol_size)

    mid_pxspc = n.zeros(3) + model_pxres * 4
    m_scl = spine.model.image_steps( model, model_pxres, mid_pxspc, sp.pxspc, \
                               sp.fwhm, blur1 = blur1, blur2 = blur2 )
    del model
    print "A", A, "m_scl s&m", m_scl.sum(), m_scl.max()
    print "Real: sum=%0.2f & max=%0.2f" % \
          (sp.rotated_data.sum(), sp.rotated_data.max())
    m_scl *= A / m_scl.sum()
    print "Model: sum=%0.2f & max=%0.2f" % (m_scl.sum(), m_scl.max())

    if type(iteration) == list:
        iteration[0] += 1
    
    if graphic:
        if display_mode == '2-D':
            overall_max = n.array(( m_scl.max(), sp.rotated_data.max() )).max()
            f0 = pl.figure()
            ax0 = f0.add_subplot(221)
            ax0.imshow(m_scl.max(2), vmin=0, vmax=overall_max)
            ax0.set_title("Model X-Y")
            ax1 = f0.add_subplot(222)
            ax1.imshow(m_scl.max(1), vmin=0, vmax=overall_max)
            ax1.set_title("Model X-Z")
            ax2 = f0.add_subplot(223)
            ax2.imshow(sp.rotated_data.max(2), vmin=0, vmax=overall_max)
            ax2.set_title("Real X-Y")
            ax3 = f0.add_subplot(224)
            ax3.imshow(sp.rotated_data.max(1), vmin=0, vmax=overall_max)
            ax1.set_title("Real X-Z")
        else:
            if scene == None:
                scene = mlab.figure()
                msrc = ArraySource()
                rsrc = ArraySource()
            else:
                if len(scene.children) >= 2:
                    msrc = scene.children[0]
                    isom = msrc.children[0].children[1]
                    isom.contour.contours = [0.]
                    #isor.contour.contours = [0.]
                    msrc.scalar_data = m_scl
                    #rsrc = scene.children[1]
                    #rsrc.scalar_data = sp.rotated_data
                    #isor = msrc.children[0].children[1]
                else:
                    scene.scene.background = (1.,1.,1.)
                    scene.scene.foreground = (0.,0.,0.)

                    msrc = ArraySource()
                    msrc.scalar_data = m_scl
                    msrc.spacing = list(sp.pxspc)
                    msrc._update_image_data_fired()
                    rsrc = ArraySource()
                    rsrc.scalar_data = sp.rotated_data
                    rsrc.spacing = list(sp.pxspc)
                    rsrc._update_image_data_fired()
                    scene.add_child(msrc)
                    msrc.add_module(Outline())
                    isom = IsoSurface()
                    isom.actor.actor.property.representation = 'wireframe'
                    isom.actor.actor.property.line_width = 1.0
                    isom.actor.mapper.scalar_visibility = False
                    msrc.add_module(isom)
                    isom.actor.actor.property.color = (0.,0.,1.)
                    scene.add_child(rsrc)
                    rsrc.add_module(Outline())
                    isor = IsoSurface()
                    isor.actor.mapper.scalar_visibility = False
                    rsrc.add_module(isor)
                    isor.actor.actor.property.color = (0.,0.5,1.)
                    isor.actor.actor.property.opacity = 0.5
                    isor.contour.contours = [rsrc.scalar_data.max() * 0.08]
                    scene.scene.z_plus_view()
                isom.contour.contours = [msrc.scalar_data.max() * 0.08]
                #msrc.scalar_data = m_scl
                #rsrc.scalar_data = sp.rotated_data
                #misoval = n.minimum(0.08, msrc.scalar_data.max())
                #risoval = n.minimum(0.08, rsrc.scalar_data.max()*0.08)
                scene.scene.z_plus_view()
                if type(iteration) == list:
                    scene.scene.save_png('%s_fit%03d.png' % \
                                         (sp.name, iteration[0]))
                    
    # calculate residual
    resid = n.abs(m_scl - sp.rotated_data).sum()
    print "Difference:", resid
    return resid + penalty

# CODE GRAVEYARD -------------------------------------------------------------

# # ---------------------------------------------------------
# # second take on the model fitting thing - uses parameter
# # space halving to find best fit point
# # ---------------------------------------------------------

#from spine.utils import read_params, write_params, complete_params, model_pxres_ar
#from pylab import plot, axhline, axvline, xlabel, ylabel
#from pylab import setp, subplot, clf
#from copy import copy
#from rebin import congrid
#from scipy.ndimage import rotate, median_filter
#from vectors import rotate_about_centre
#from coordhandling import centre_of_mass, align_stacks_big_by_offs

# def midtup(tup):
#     return (tup[0] + tup[1]) / 2.


# def flteq(a,b):
#     '''used to check if two floats are the same (within pixel res)'''
#     # dodgy
#     return n.abs(a - b) < model_pxres / 2.


# def tupeq(a,b):
#     '''checks if both tuple''s values are the same (within pixel res)'''
#     return flteq(a[0], b[0]) and flteq(a[1], b[1])


# def test_model(hd, nd, spine):
#     '''creates a model based on the supplied parameters and returns
#     the neck and head intensities'''
#     (sp, ori, neck, head, skl) = \
#          spine_model( px_res = model_pxres, \
#                       neck_length = spine.params['neck_length'], \
#                       dend_diam = spine.params['dend_diam'], \
#                       angle_roundz = spine.params['angle_roundz'], \
#                       angle_toz = spine.params['angle_toz'], \
#                       angle_den = spine.params['angle_den'], \
#                       head_diam = hd, \
#                       neck_diam = nd )
#     del skl
#     (spi, psf) = spine_image( sp_mdl = sp, pxr = model_pxres, \
#                               fwhms = spine.fwhm )
#     del sp, psf
#     # - correct point co-ords for convolve's padding
#     # -- no longer needed as using processing.fftconvolve2
#     #psfsize = array( psf.shape )
#     #ori += (psfsize - 1)/2
#     #neck += (psfsize - 1)/2
#     #head += (psfsize - 1)/2
    
#     mdiabs = spi[ori[0], ori[1], ori[2]]
#     mni = spi[neck[0], neck[1], neck[2]] / mdiabs.astype(float)
#     mhi = spi[head[0], head[1], head[2]] / mdiabs.astype(float)
#     plot( (mhi,), (mni,), 'bo' )
#     # return (mhi, mni, spi, head, neck, sp, psf)
#     return (mhi, mni)
# #    print "    -> model intensities - head %0.3f neck %0.3f" % (mhi, mni)


# def fit_head2( spine, nd, estrange, error=0.01 ):
#     '''fits spine model neck to given match given intensity
#     (using head and neck diams if supplied) '''

#     print "*** fitting HEAD with neck diam", nd
# #    print "    starting at range (%0.3f, %0.3f)" % (estrange[0], estrange[1])

#     hd = 0
#     done = False
#     print "    Btm | Top  |  Val  |    I   |   Ref   "
#     while not done:
#         oldhd = hd

#         hd = midtup(estrange)
#         hd -= hd % model_pxres
#         print "H (%5.3f, %5.3f) %5.3f ->" % (estrange[0], estrange[1], hd),
#         if flteq(hd, oldhd):
#             print "Terminated because repeating."
#             done = True
#         else:
#             #print "    checking middle (%0.3f) of range (%0.3f, %0.3f)" % \
#             #      (hd, estrange[0], estrange[1])
#             (mhi,mni) = test_model( hd, nd, spine )
#             print "%5.3f" % (mhi),

#             if mhi > spine.params['head_int']:
#                 print "> %5.3f" % spine.params['head_int']
#                 estrange[1] = hd
#                 #print "       head int %0.3f too hi (> %f), " \
#                 #      "adjusting top of range" % \
#                 #      (mhi, spine.params['head_int'])
#             else:
#                 print "< %5.3f" % spine.params['head_int']
#                 estrange[0] = hd
#                 #print "       head int %0.3f too low (< %f), " \
#                 #      "adjusting bottom of range" % \
#                 #      (mhi, spine.params['head_int'])

#             # problem is that we get repetition before closeness
#             # need to change finishing criterion
#             done = n.abs( mhi - spine.params['head_int'] ) < error
#             if done:
#                 print "Terminated because got close."
#             # notdone = abs( mhi - spine.params['head_int'] ) > \
#             #          error
#     return (hd, mhi, mni)


# def fit_neck2( spine, hd, estrange, error=0.01 ):
#     '''fits spine model neck to given match given intensity
#     (using head and neck diams if supplied) '''
#     print "*** fitting NECK with head diam", hd
# #    print "    starting at range (%0.3f, %0.3f)" % (estrange[0], estrange[1])

#     nd = 0
#     done = False
#     print "    Btm | Top  |  Val  |    I   |   Ref   "
#     while not done:
#         oldnd = nd

#         nd = midtup(estrange)
#         nd -= nd % model_pxres
#         print "N (%5.3f, %5.3f) %5.3f ->" % (estrange[0], estrange[1], nd),
#         if flteq(nd, oldnd):
#             done = True
#             print "Terminate because repeating."
#         else:
#             #print "    checking middle (%0.3f) of range (%0.3f, %0.3f)" % \
#             #        (nd, estrange[0], estrange[1])
#             (mhi,mni) = test_model(hd, nd, spine)
#             print "%5.3f" % (mni),

#             if mni > spine.params['neck_int']:
#                 estrange[1] = nd
#                 print "> %5.3f" % spine.params['neck_int']
#                 #print "       neck int %0.3f too hi (> %f), " \
#                 #      "adjusting top of range" % \
#                 #      (mni, spine.params['neck_int'])
#             else:
#                 print "< %5.3f" % spine.params['neck_int']
#                 estrange[0] = nd
#                 #print "       neck int %0.3f too low (< %f), " \
#                 #      "adjusting bottom of range" % \
#                 #      (mni, spine.params['neck_int'])
            
#             # problem is that we get repetition before closeness
#             # need to change finishing criterion
#             done = n.abs( mni - spine.params['neck_int'] ) < error
#             if done:
#                 print "Terminated because got close."

#             # notdone = abs( mni - spine.params['neck_int'] ) > \
#             #          error

#     return (nd, mhi, mni)


# def one_fit_loop( spine, estrange, ishead, otherparm, goal):
#     '''Performs one loop of fitting - improving estimate of range within
#     which param, ishead tells whether head (= True) or neck (= False), lies.
#     Tests mid-point of estrange and checks if >|< desired intensity and
#     updates estrange accordingly'''
#     parm = midtup(estrange)
#     parm -= parm % model_pxres
#     if ishead:
#         print "H",
#     else:
#         print "N",
#     print "(%5.3f, %5.3f) %5.3f -> " % (estrange[0], estrange[1], parm),
#     if ishead:
#         (mhi,mni) = test_model(parm, otherparm, spine)
#         mi = mhi
#     else: # is neck
#         (mhi,mni) = test_model(otherparm, parm, spine)
#         mi = mni
#     print "%5.3f" % (mi),
#     newrange = copy(estrange)
#     if mi > goal:
#         newrange[1] = parm
#         print "> %5.3f" % goal
#     else: # mi < goal
#         print "< %5.3f" % goal
#         newrange[0] = parm
#     return newrange, (mhi, mni)


# # ----------------------------------------------------------------------------


# def batch_fit( filename, psf_fwhms ):

#     pars = read_params( filename, '*' )
#     fit_marker = '_fit'


#     # extract list of spnames from pars
#     spnames = []
#     for par in pars:
#         spnames.append( par['spname'] )

#     # parse through parameter sets
#     for par in pars:
#         spn = par['spname']

#         # check if done
#         if fit_marker not in spn:
#             if spn + fit_marker in spnames:
#                 print spn + ' is already fitted.'
#             else:
#                 print spn + ' needs fitting.'
#                 sp = spine.Spine( spn, model_pxres_ar, psf_fwhms )
#                 sp.params = par
#                 sp.fit_model()
#                 write_params( filename, spn + '_fit', \
#                               complete_params( sp.params ))
#                 #fit_model2( filename, spn )


# # ----------------------------------------------------------------------------


# def downsample_model(model, orig_pxres, new_pxres):
#     cal_factor = orig_pxres / new_pxres
#     newdims = n.asarray(model.shape) * cal_factor
#     return congrid(model, newdims)


# def align_stacks_minimum(a, oa, b, ob):
#     '''take common region between a and b, centred so
#     that oa and ob overlap'''
#     oa = n.asarray(oa)
#     ob = n.asarray(ob)
    
#     min_lc = n.vstack((oa, ob)).min(0)
#     min_uc = n.vstack((a.shape - oa, \
#                        b.shape - ob)).min(0)
#     min_shape = min_lc + min_uc
#     a_slice = [ slice(i,j) for (i,j) in \
#                 zip(list(oa - min_lc), \
#                     list(oa + min_uc)) ]
#     min_a = a[a_slice]
#     b_slice = [ slice(i,j) for (i,j) in \
#                 zip(list(ob - min_lc), \
#                     list(ob + min_uc)) ]
#     min_b = b[b_slice]
#     return min_a, min_b

    
# #-----------------------------------------------------------------------------

# def compare_models(sp):
#     '''just for use in debugging test_accuracy spines (where orient_den == 0)
#     '''
#     # check params are valid
#     if n.isnan(sp.params['head_diam']):
#         sp.params['head_diam'] = 0.33
#     if n.isnan(sp.params['neck_diam']):
#         sp.params['neck_diam'] = 0.1
            
#     # fit dendrite vector
#     dvec_a, dvec = fit_dendrite_vector(sp.source_data)
#     # define dvec as pointing towards +ve z
#     if n.sign(dvec[2]) == -1:
#         dvec *= -1

#     sp.params['orient_den'] = 0. # just when debugging fake spines

#     head_diam = sp.params['head_diam']
#     neck_diam = sp.params['neck_diam']
#     neck_length = sp.params['neck_length']
#     dend_diam = sp.params['dend_diam']
#     angle_roundz = sp.params['angle_roundz']
#     angle_toz = sp.params['angle_toz']

#     dvec *= sp.pxspc / model_pxres
#     dvec_a *= sp.pxspc / model_pxres

#     # pre-calculate some constant stuff for use in the loops
#     model_size = (n.asarray(sp.source_data.shape) * \
#                   sp.pxspc / model_pxres).round()
#     model_inds = n.indices(model_size)

#     t = 0.
#     o = dvec_a + t * dvec
#     angle_den = n.arctan( dvec[2] / n.sqrt( dvec[0]**2 + dvec[1]**2 ) )
#     print "Angle den", angle_den
    
#     print "  hd=", head_diam, "nd=", neck_diam, "\n", \
#           "  nl=", neck_length, "\n", \
#           "  origin", o

#     mbits = spine_model_new( \
#         px_res = model_pxres, \
#         pars = {'origin_x'    : o[0], \
#                 'origin_y'    : o[1], \
#                 'origin_z'    : o[2], \
#                 'angle_den'   : angle_den, \
#                 'dend_diam'   :  sp.params['dend_diam'], \
#                 'head_diam'   : head_diam, \
#                 'neck_diam'   : neck_diam, \
#                 'neck_length' : neck_length, \
#                 'angle_roundz':  sp.params['angle_roundz'], \
#                 'angle_toz'   :  sp.params['angle_toz']
#                 }, \
#         vol_size = model_size,
#         vol_inds = model_inds)
#     return mbits[0], sp.source_data

# def model_merit(p, vol_size, sp, flip=False, graphic=False):
#     '''Returns sum of squares of difference between model data and source
#     data, aligned at origin. Takes smallest common region for comparison.'''
#     # alternative method might have been to zero pad as required to
#     # make all of both volumes common
    
#     #    head_diam, neck_diam, neck_length, A = p
#     head_diam, neck_diam, neck_length, A, ox, oy, oz = p
#     head_diam = n.abs(head_diam)
#     neck_diam = n.abs(neck_diam)
#     neck_length = n.abs(neck_length)
    
#     #print "Time:", time.asctime()
#     ts = [time.time()] # t0
#     print "%s hd=%6.4f nd=%6.4f nl=%6.4f" % (time.asctime(), head_diam, neck_diam, neck_length)
#           #"  A=", A, "\n", \
#           #"  origin", ox, oy, oz

#     # rotate real data in x-y so that vector points up (if it has z extent)
#     # along y axis i.e. (0,1,0)
#     deg = 180. / n.pi

#     # x-y rotation of real data to align to models
#     #real_c = n.round( n.asarray( sp.source_data.shape )/2. )
#     real = rotate( sp.source_data, \
#                    sp.params['orient_den'] * deg, \
#                    axes=(0,1), reshape=True)
#     ts.append( time.time() ) # t1
    
#     mbits = spine_model_new( \
#         px_res = model_pxres, \
#         pars = {'origin_x'    :ox, \
#                 'origin_y'    :oy, \
#                 'origin_z'    :oz, \
#                 'angle_den'   :sp.params['angle_den'], \
#                 'dend_diam'   :sp.params['dend_diam'], \
#                 'head_diam'   :head_diam, \
#                 'neck_diam'   :neck_diam, \
#                 'neck_length' :neck_length, \
#                 'angle_roundz':sp.params['angle_roundz'], \
#                 'angle_toz'   :sp.params['angle_toz']
#                 }, \
#         vol_size = vol_size)

#     ts.append(time.time()) # t2
#     m0 = mbits[0]
#     del mbits

#     mid_pxspc = n.zeros(3) + model_pxres * 2
#     m_scl = spine_image_steps( m0, model_pxres, mid_pxspc, sp.pxspc, sp.fwhm )

#     ts.append(time.time()) # t3

#     m_scl_sum = m_scl.sum()
#     m_scl *= A / m_scl_sum

#     if flip:
#         m_scl = m_scl[:,::-1]
#         title = "(flipped)"
#     else:
#         title = ""
        
#     if graphic:
#         overall_max = n.array(( m_scl.max(), real.max() )).max()
#         f0 = pl.matshow(m_scl.max(2), vmin=0, vmax=overall_max)
#         ax0 = f0.axes[0]
#         ax0.set_title("Model" + title)
        
#         f1 = pl.matshow(real.max(2), vmin=0, vmax=overall_max)
#         ax1 = f1.axes[0]
#         ax1.set_title("Real")

#     # calculate error
#     print "Sum-of-squares:", ((m_scl - real)**2).sum()
#     ts.append(time.time()) # t4

#     print "Rotate real", ts[1] - ts[0]
#     print "Make model", ts[2] - ts[1]
#     print "Image model + downsample", ts[3] - ts[2]
#     print "Scale [+ flip] + sum of sqs", ts[4] - ts[3]
    
#     return (m_scl - real).flatten()




# def model_merit_sq(p, sp):
#     err = model_merit(p, sp)
#     return (err**2).sum()


# def fit_spine_angle_toz(p, d_a, d_D, vol_size, vol_inds, sp, \
#                     flip=False, graphic=False):
#     '''Fits unblurred cylinders to spine/dendrite, to get
#     spine angles - solves for t, angle_roundz and angle_toz'''
#     print ""
#     print "Time:", time.asctime()

#     angle_toz, t = p
#     o = d_a + t * d_D
#     angle_den = n.arctan( d_D[2] / n.sqrt( d_D[0]**2 + d_D[1]**2 ) )
#     print " origin:", o, "\n angle to z:", angle_toz

#     real = congrid(sp.rotated_data, vol_size, method='spline')

#     model = 




# def fit_spine_angle(p, d_a, d_D, vol_size, vol_inds, sp, \
#                     flip=False, graphic=False):
#     '''Fits unblurred cylinders to spine/dendrite, to get
#     spine angles - solves for t, angle_roundz and angle_toz'''
#     print ""
#     print "Time:", time.asctime()

#     angle_roundz, angle_toz, t = p
#     o = d_a + t * d_D
#     angle_den = n.arctan( d_D[2] / n.sqrt( d_D[0]**2 + d_D[1]**2 ) )
#     print " origin:", o, "\n angle round z:", angle_roundz, \
#           "\n angle to z:", angle_toz

#     real = congrid(sp.rotated_data, vol_size, method='spline')
#     mbits = spine_model_new( \
#         px_res = model_pxres, \
#         pars = {'origin_x'    : o[0], \
#                 'origin_y'    : o[1], \
#                 'origin_z'    : o[2], \
#                 'angle_den'   : angle_den, \
#                 'dend_diam'   : 0.5, \
#                 # just needs to be smaller
#                 # than real dendrite
#                 'head_diam'   : 0.09, \
#                 'neck_diam'   : 0.09, \
#                 'neck_length' : 3.0, \
#                 # just has to be longer than spine
#                 'angle_roundz': angle_roundz, \
#                 'angle_toz'   : angle_toz \
#                 }, \
#         vol_size = vol_size,
#         vol_inds = vol_inds)
    
#     m0 = mbits[0]
#     del mbits
    
#     m0 *= real.max() / m0.max()
    
#     if flip:
#         m0 = m0[:,::-1]
#         title = "(flipped)"
#     else:
#         title = ""

#     if graphic:
#         overall_max = n.array(( m0.max(), real.max() )).max()
#         f0 = pl.matshow(m0.max(2), vmin=0, vmax=overall_max)
#         ax0 = f0.axes[0]
#         ax0.set_title("Model" + title)
        
#         f1 = pl.matshow(real.max(2), vmin=0, vmax=overall_max)
#         ax1 = f1.axes[0]
#         ax1.set_title("Real")

#     err = (m0 - real).flatten()
#     print "Sum-of-squares:", (err**2).sum()
#     return err



# def acceptable():
#     '''returns acceptable tolerance for fitting procedures
#     (distance in microns)'''
#     return 0.002




# def fit_head(hi, nl, dd, arz, atz, \
#              hd = 0.25, nd = 0.1 ):
#     '''fits spine head to match given head intensity'''
#     print "Fitting HEAD with neck diam", nd

#     step_size = 0.195
#     step_scale = 0.5
#     wasbig = None
#     notdone = True
#     pxres = 0.015

#     subplot(121)
#     axhline(hi, color='r')

#     while notdone:

#         (sp, neck_ori, nec_cen, head, skl) = \
#              spine_model( px_res = pxres, \
#                           neck_length = nl, \
#                           dend_diam = dd, \
#                           angle_roundz = arz, \
#                           angle_toz = atz, \
#                           head_diam = hd, \
#                           neck_diam = nd )
#         (spi, psf) = spine_image( sp, pxr = pxres() )

#         # correct point co-ords for convolve's padding
#         psfsize = array( psf.shape )
#         neck_ori += (psfsize - 1)/2
#         nec_cen += (psfsize - 1)/2
#         head += (psfsize - 1)/2

#         mdiabs = spi[neck_ori[0], neck_ori[1], neck_ori[2]]
#         mhi = spi[head[0], head[1], head[2]] / mdiabs.astype(float)

#         plot( (hd,), (mhi,), 'ro' )
#         print "  guess: d=", hd, "i=", mhi

#         # correct guess value (for head) based on direction of error
#         # can I be smarter and get some idea of magnitude of the error?
#         notdone = (abs( mhi - hi ) > acceptable()) and (step_size > pxres)

#         if notdone:
#             if mhi > hi:
#                 # model is too big
#                 dr = -1
#                 if not wasbig:
#                     step_size *= step_scale
#                 wasbig = True
#             else:
#                 dr = 1
#                 if wasbig:
#                     step_size *= step_scale
#                 wasbig = False

#             hd += step_size * dr
#             print "  stepping by", step_size * dr, "to d=", hd


#     print "fitted head diameter", hd
#     return hd




# def fit_neck(ni, nl, dd, arz, atz, \
#              hd = 0.25, nd = 0.1 ):
#     '''fits spine head to match given head intensity'''
#     print "Fitting NECK with head diam", hd

#     step_size = 0.195
#     step_scale = 0.5
#     wasbig = None
#     notdone = True
#     pxres = 0.015
    
#     subplot(122)
#     axhline(ni, color='b')

#     while notdone:
#         (sp, neck_ori, nec_cen, head, skl) = \
#              spine_model( px_res = pxres, \
#                           neck_length = nl, \
#                           dend_diam = dd, \
#                           angle_roundz = arz, \
#                           angle_toz = atz, \
#                           head_diam = hd, \
#                           neck_diam = nd )
#         (spi, psf) = spine_image( sp, pxr = pxres() )

#         # correct point co-ords for convolve's padding
#         psfsize = array( psf.shape )
#         neck_ori += (psfsize - 1)/2
#         nec_cen += (psfsize - 1)/2
#         head += (psfsize - 1)/2

#         mdi_abs = spi[neck_ori[0], neck_ori[1], neck_ori[2]]
#         mni = spi[nec_cen[0], nec_cen[1], nec_cen[2]] / mdi_abs.astype(float)

#         print "  guess: d=", nd, "i=", mni
#         plot( (hd,), (mni,), 'bo' )

#         # correct guess value (for head) based on direction of error
#         # can I be smarter and get some idea of magnitude of the error?
#         notdone = (abs( mni - ni ) > acceptable()) and (step_size > pxres)

#         if notdone:
#             if mni > ni:
#                 # model is too big
#                 dr = -1
#                 if not wasbig:
#                     step_size *= step_scale
#                 wasbig = True
#             else:
#                 dr = 1
#                 if wasbig:
#                     step_size *= step_scale
#                 wasbig = False

#             nd += step_size * dr
#             print "  stepping by", step_size * dr, "to d=", nd
            

#     print "fitted neck diameter", nd
#     return nd




# def fit_model(fname, spname):
#     '''fits head and neck diams to data supplied in fname'''
    
#     # read data file
#     vals = read_params( fname, spname )
#     nl = vals['neck_length']
#     dd = vals['dend_diam']
#     arz = vals['angle_roundz']
#     if arz > pi/2:
#         arz = arz - pi
#     atz = vals['angle_toz']
#     ni = vals['neck_int']
#     hi = vals['head_int']

#     print "DESIRED VALUES: neck int", ni, "head int", hi

#     hplot = subplot(121)
#     setp(hplot, xlim=(0.1,0.4))
#     nplot = subplot(122)
#     setp(nplot, xlim=(0.1,0.4))

#     # pick starting values
#     ndfit0 = 0.1
#     hdfit0 = 0.25
#     error = 0.005

#     notdone = True

#     while notdone:
#         hdfit1 = fit_head( hi, nl, dd, arz, atz, nd = ndfit0, hd = hdfit0 )
#         ndfit1 = fit_neck( ni, nl, dd, arz, atz, hd = hdfit1, nd = ndfit0 )
#         notdone = (abs(ndfit0 - ndfit1) > error) \
#                   or (abs(hdfit0 - hdfit1) > error)
#         (hdfit0, ndfit0) = (hdfit1, ndfit1)

#     print "Fitted params:", hdfit1, ndfit1

#     vals['head_diam'] = hdfit1
#     vals['neck_diam'] = ndfit1
#     write_params( fname, spname + '_fitted', complete_params( vals ) )


# from model_merit

#     frpr = pl.matshow(sp.source_data.max(2))
#     axrpr = frpr.axes[0]
#     axrpr.set_title("Whole Unrotated Real")
#     real_or = sp.user_pts[2]
#     axrpr.plot( (real_or[1] + 0.5, ), (real_or[0] + 0.5, ), 'bo')

#     # calculate rotated origin point
#     real_unroto = sp.user_pts[2]
#     real_rotc = n.round( n.asarray( real_rot.shape )/2. )
    
#     # not sure about this next line
#     real_roto_xy = rotate_about_centre( \
#         real_unroto[0:2], real_c[0:2], sp.params['orient_den'] - n.pi/2. ) \
#         - real_c[0:2] + real_rotc[0:2]
#     real_o = n.asarray( (list(real_roto_xy) + [real_unroto[2]]) )

# #     fr = pl.matshow(real.max(2))
# #     axr = fr.axes[0]
# #     axr.set_title("Whole Rotated Real")
# #     axr.plot( (real_o[1] + 0.5, ), (real_o[0] + 0.5, ), 'bo')

#     fm = pl.matshow(m_scl.max(2))
#     axm = fm.axes[0]
#     axm.set_title("Whole Model")
#     axm.plot( (m_scl_orig[1] + 0.5, ), (m_scl_orig[0] + 0.5, ), 'bo')

#     ffl = pl.matshow(mfl_scl.max(2))
#     axfl = ffl.axes[0]
#     axfl.set_title("Whole Model Flipped")
#     axfl.plot( (mfl_scl_orig[1] + 0.5, ), (mfl_scl_orig[0] + 0.5, ), 'bo')

#     # now handle model data
#     m_origin = n.asarray(m_orig)
#    mfl_orig = n.array((m_orig[0], m.shape[1] - m_orig[1], m_orig[2]))
