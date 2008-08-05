import numpy as n
from vectors import unitvec
from primitives import cylinder, sphere, clip_plane
# from processing import fftconvolve2 as fftconvolve
from scipy.signal import fftconvolve
from gaussian import gauss3d
from myutils import dparse, pdbg
from rebin import congrid, rebin_average
import spine.utils
from scipy import weave

'''angles behave as follows:
angle_toz +ve = down
          -ve = up
angle_roundz +ve = towards up end (+ve angle_den)
             -ve = towards down end (-ve angle_den)
'''

def build( px_res = spine.utils.model_pxres, \
           pars = {'origin_x' : 150, \
                   'origin_y' : 100, \
                   'origin_z' : 80}, \
           vol_size = (200,200,160), \
           need_skl = False):
    '''an implementation of spine_model_new using weave.inline'''
    pars = dparse({ 'origin_x'     : 0., \
                    'origin_y'     : 0., \
                    'origin_z'     : 0., \
                    'angle_den'    : 0., \
                    'dend_diam'    : 1., \
                    'neck_diam'    : 0.5, \
                    'neck_length'  : 0.5, \
                    'head_diam'    : 1. , \
                    'angle_roundz' : 0., \
                    'angle_toz'    : 0.
                    }, pars)

    #print "Size: %d, %d, %d" % (vol_size[0],vol_size[1],vol_size[2])
    ox = pars['origin_x']
    oy = pars['origin_y']
    oz = pars['origin_z']
    #print "Origin: %d, %d, %d" % (ox,oy,oz)
    angle_den = pars['angle_den']
    dend_diam = pars['dend_diam']
    neck_diam = pars['neck_diam']
    neck_length = pars['neck_length']
    head_diam = pars['head_diam']
    angle_roundz = pars['angle_roundz']
    angle_toz = pars['angle_toz']

    if px_res == None:
        raise ValueError('Pixel resolution must be defined.')
    assert( type( px_res ) == float )

    assert( len(vol_size) == 3 )
    xsz,ysz,zsz = [int(x) for x in list(vol_size)]

    dpx = dend_diam / px_res
    nwpx = neck_diam / px_res
    hpx = head_diam / px_res
    nlpx = neck_length / px_res
    dend_vec = (0., n.cos(angle_den),n.sin(angle_den))
    dDx, dDy, dDz = [float(x) for x in list(dend_vec)]
    nec_vec_dir = n.array( (n.cos( angle_roundz ) * n.cos( angle_toz ), \
                            # x component
                            n.sin( angle_roundz ) * n.cos( angle_toz ), \
                            # y component
                            n.sin( angle_toz )
                            # z component
                            ) )
    nDx, nDy, nDz = [float(x) for x in list(nec_vec_dir)]
    nec_in_den = (dpx/2.) / n.sqrt( 1 - n.dot(nec_vec_dir, dend_vec)**2 )
    nec_vec_len = (nlpx + nec_in_den + hpx/2.)
    nec_vec = nec_vec_len * nec_vec_dir
    neck_ori = n.array( ( ox, oy, oz ) )
    head_pos = n.asarray(neck_ori - nec_vec)
    hx, hy, hz = [float(x) for x in list(head_pos)]
    dr = dpx / 2.
    nr = nwpx / 2.
    hr = hpx / 2.
    ar = n.empty((xsz,ysz,zsz), dtype=n.uint8)
    code = '''
    #line 14
    bool inden, inneck, inhead, pl0, pl1;
    double Yi, Yj, Yk;
    double tmp, ddi, ndi;

    // iterate through volume - should need to do only once
    for (int i = 0; i<xsz; i++)
        for (int j = 0; j<ysz; j++)
            for (int k = 0; k<zsz; k++)
            {
                // make dendrite
                Yi = (double)i - ox; Yj = (double)j - oy; Yk = (double)k - oz;
                ddi = dDx * Yi + dDy * Yj + dDz * Yk;
                tmp = sqrt( pow(Yi - ddi * dDx, 2) +
                  pow(Yj - ddi * dDy, 2) + pow(Yk - ddi * dDz, 2) );
                inden = (tmp < (double)dr);

                // make neck - whole length
                ndi = nDx * Yi + nDy * Yj + nDz * Yk;
                tmp = sqrt( pow(Yi - ndi * nDx, 2) +
                  pow(Yj - ndi * nDy, 2) + pow(Yk - ndi * nDz, 2) );
                // clip neck - dendrite side
                pl0 = (nDx * ((double)i - ox) + nDy * ((double)j - oy) +
                  nDz * ((double)k - oz)) <= 0;
                //pl1 = (1==1);
                // clip - head side
                pl1 = (-nDx * ((double)i - hx) - nDy * ((double)j - hy) -
                  nDz * ((double)k - hz)) <= 0;
                inneck = (tmp < nr) & pl0 & pl1;

                // make head
                tmp = sqrt(pow((double)i - hx, 2) + pow((double)j - hy, 2)
                  + pow((double)k - hz, 2));
                inhead = (tmp < hr);
                
                // assign point
                ar(i,j,k) = (inden | inneck | inhead);
            }
    '''
    weave.inline( code, ['ar', \
                         'xsz','ysz','zsz', \
                         'ox', 'oy', 'oz', \
                         'dDx', 'dDy', 'dDz', \
                         'nDx', 'nDy', 'nDz', \
                         'hx', 'hy', 'hz', \
                         'dr', 'nr', 'hr'], \
                  type_converters=weave.converters.blitz )
    return ar

def image_steps(model, model_pxspc, mid_pxspc, \
                final_pxspc, final_fwhms, \
                blur1 = None, blur2 = None):
    
    model = n.cast[float](model)
    model_pxspc = n.asarray(model_pxspc)
    mid_pxspc = n.asarray(mid_pxspc)
    final_pxspc = n.asarray(final_pxspc)
    final_fwhms = n.asarray(final_fwhms)    

    mid_fwhms = mid_pxspc * 2
    mid_fwhms_px = mid_fwhms / model_pxspc
    if blur1 == None:
        blur1 = gauss3d(mid_fwhms_px[0], mid_fwhms_px[1], \
                        mid_fwhms_px[2], 1e-3)
    # pre-filter
    mid_big = fftconvolve( model, blur1, mode='same')
    # downsample to mid_pxspc
    scale_factor = model_pxspc / mid_pxspc
    middims = (n.asarray(model.shape) * scale_factor).round()
    mid = congrid(mid_big, middims, method='linear')

    second_fwhms = n.sqrt(final_fwhms**2 - mid_fwhms**2) / mid_pxspc
    if blur2 == None:
        blur2 = gauss3d(second_fwhms[0], second_fwhms[1], \
                        second_fwhms[2], 1e-3)
    final_big = fftconvolve( mid, blur2, mode='same')
    scale_factor2 = mid_pxspc / final_pxspc
    finaldims = (middims * scale_factor2).round()
    final = congrid(final_big, finaldims, method='linear')
    return final

def build_with_values( px_res = None, \
                             dend_diam = 1.5, \
                             head_diam = 1.5, \
                             neck_diam = 0.5, \
                             neck_length = 4., \
                             angle_roundz = 0.,
                             angle_toz = 0., \
                             space = 25, \
                             head_val = 1., \
                             neck_val = 1., \
                             dend_val = 1.,
                             dend_len = 0.):
    '''used for generating colour-coded fit model figure - should be updated
    to match spine_model_weave implementation'''
    if px_res == None:
        raise ValueError('pxr must be defined')

    if angle_roundz > (n.pi/2):
        print "Parameter error: angle_roundz <= pi/2 please"
        return None
    if angle_toz > (n.pi/2):
        print "Parameter error: angle_toz <= pi/2 please"
        return None

    dpx = dend_diam / px_res
    nwpx = neck_diam / px_res
    hpx = head_diam / px_res
    nlpx = neck_length / px_res
    dlpx = dend_len / px_res
    
    # added dpx/2., hpx/2. to nlpx <= changing definition of neck length
    print "Error here and below"
    nec_vec = (nlpx + dpx/2. + hpx/2.) * \
              n.array( (n.cos( angle_roundz ) * n.cos( angle_toz ), \
                      # x component
                      n.sin( angle_roundz ) * n.cos( angle_toz ), \
                      # y component
                      n.sin( angle_toz )
                      # z component
                      ))

    vol_width = dpx / 2. + abs( nec_vec[0] ) + hpx / 2. + space * 2
    vol_height = n.array( (dpx, nwpx, hpx) ).max() \
                 + abs( nec_vec[2] ) + space * 2
    vol_depth = n.array( (dpx, nwpx, hpx, dlpx) ).max() \
               + abs( nec_vec[1] ) + space * 2
    sz = n.array( (vol_width, vol_depth, vol_height) ).round().astype(int)

    neck_ori = n.array( (vol_width - (space + dpx / 2.), \
                        vol_depth / 2. + nec_vec[1] / 2., \
                        vol_height / 2. + nec_vec[2] / 2.) \
                      ).round(0).astype(int)

    head_pos = n.array(neck_ori - nec_vec).astype(int)

    # make dendrite
    vol = cylinder( sz, neck_ori, (0,1,0), n.array(dpx / 2.).round() \
                    ).astype(float) * dend_val
    skvol = cylinder( sz, neck_ori, (0,1,0), 1 ).astype(float) \
            * dend_diam

    # make neck
    nvol = cylinder( sz, neck_ori, unitvec( nec_vec ), \
                     n.array(nwpx / 2.).round() )
    nvol = nvol & clip_plane( sz, neck_ori, -unitvec( nec_vec ) )
    nvol = nvol & clip_plane( sz, head_pos, unitvec( nec_vec ) )
    # vol[n.where(nvol == True)] = nvol[n.where(nvol == True)].astype(float) \
    vol[nvol == True] = nvol[nvol == True].astype(float) \
                               * neck_val
    del nvol
    
    sknvol = cylinder( sz, neck_ori, unitvec( nec_vec ), 1);
    sknvol = sknvol & clip_plane( sz, neck_ori, -unitvec( nec_vec ) )
    sknvol = (sknvol & clip_plane( sz, head_pos, unitvec( nec_vec ) ) \
              ).astype(float)
    sknvol *= neck_diam
    # sknvol[n.where(skvol > 0)] = skvol[n.where(skvol > 0)]
    sknvol[skvol > 0] = skvol[skvol > 0]
    del skvol

    # make head
    hvol = sphere( sz, head_pos, n.array( hpx / 2.).round() )
    # vol[n.where(hvol == True)] = hvol[n.where(hvol == True)].astype(float) \
    vol[hvol == True] = hvol[hvol == True].astype(float) \
                               * head_val
    del hvol

    # nlpx... removed - 0.5 * dpx - 0.5 * hpx in-line with new definition nlpx
    dist_to_nec_cen = nlpx / 2 + 0.5 * dpx
    nec_cen = n.array(neck_ori - unitvec( neck_ori - head_pos ) \
                    * dist_to_nec_cen).astype(int)

#     vol[nec_cen[0], nec_cen[1], nec_cen[2]] = 2
#     vol[head_pos[0], head_pos[1], head_pos[2]] = 2
#     vol[neck_ori[0], neck_ori[1], neck_ori[2]] = 2

    sknvol[head_pos[0], head_pos[1], head_pos[2]] = head_diam

    return vol, neck_ori, nec_cen, head_pos, sknvol

# ----------------------------------------------------------------------------

# def spine_model( px_res = None, \
#                  dend_diam = 1.5, \
#                  head_diam = 1., \
#                  neck_diam = 0.5, \
#                  neck_length = 2., \
#                  angle_roundz = 0.,
#                  angle_toz = 0., \
#                  angle_den = 0., \
#                  space = 25,
#                  dend_len = 0.,
#                  need_skl = False):

#     if px_res == None:
#         raise ValueError('pxr must be defined')

#     if angle_roundz > (n.pi/2.):
#         print "Parameter error: angle_roundz <= pi/2 please"
#         return None
#     if angle_toz > (n.pi/2.):
#         print "Parameter error: angle_toz <= pi/2 please"
#         return None

#     dpx = dend_diam / px_res
#     #    print "dpx", dpx, dpx/2.
#     nwpx = neck_diam / px_res
#     hpx = head_diam / px_res
#     #    print "hpx", hpx, hpx/2.
#     nlpx = neck_length / px_res
#     #    print "nlpx", nlpx
#     dlpx = dend_len / px_res

#     dend_vec = (0, n.cos(angle_den),n.sin(angle_den))
#     nec_vec_dir = n.array( (n.cos( angle_roundz ) * n.cos( angle_toz ), \
#                             # x component
#                             n.sin( angle_roundz ) * n.cos( angle_toz ), \
#                             # y component
#                             n.sin( angle_toz )
#                             # z component
#                             ) )

#     # added dpx/2., hpx/2. to nlpx <= changing definition of neck length
#     nec_in_den = (dpx/2.) / n.sqrt( 1 - n.dot(nec_vec_dir, dend_vec)**2 )
#     nec_vec_len = nlpx + nec_in_den + hpx/2.

#     # print "nec_vec", nec_vec_len, n.linalg.norm(nec_vec_dir)
#     nec_vec = nec_vec_len * nec_vec_dir

#     # print "nec_vec", nec_vec

#     vol_width = dpx / 2. + abs( nec_vec[0] ) + hpx / 2. + space * 2
#     vol_height = n.array( (dpx, nwpx, hpx) ).max() \
#                  + abs( nec_vec[2] ) + space * 2
#     vol_depth = n.array( (dpx, nwpx, hpx, dlpx) ).max() \
#                + abs( nec_vec[1] ) + space * 2
#     sz = n.array( (vol_width, vol_depth, vol_height) ).round().astype(int)

#     neck_ori = n.array( (vol_width - (space + dpx / 2.), \
#                         vol_depth / 2. + nec_vec[1] / 2., \
#                         vol_height / 2. + nec_vec[2] / 2.) \
#                       ).round(0).astype(n.int)
#     # print "neck ori:", neck_ori

#     head_pos = n.asarray(neck_ori - nec_vec).round(0).astype(int)

#     # print "head_pos:", head_pos

#     # make dendrite
#     vol = cylinder( sz, neck_ori, dend_vec, n.array(dpx / 2.).round() )
    
#     if need_skl:
#         skvol = cylinder( sz, neck_ori, dend_vec, 1 ).astype(float) \
#                 * dend_diam
        
#     # make head
#     vol = n.asarray(vol | sphere( sz, head_pos, n.asarray( hpx / 2.).round() ))

#     # make neck
#     nvol = cylinder( sz, neck_ori, unitvec( nec_vec ), \
#                      n.array(nwpx / 2.).round() )
#     nvol = nvol & clip_plane( sz, neck_ori, -unitvec( nec_vec ) )
#     nvol = nvol & clip_plane( sz, head_pos, unitvec( nec_vec ) )

#     if need_skl:
#         sknvol = cylinder( sz, neck_ori, unitvec( nec_vec ), 1);
#         sknvol = sknvol & clip_plane( sz, neck_ori, -unitvec( nec_vec ) )
#         sknvol = (sknvol & clip_plane( sz, head_pos, unitvec( nec_vec ) ) \
#                   ).astype(float)
#         sknvol *= neck_diam
#         # sknvol[n.where(skvol > 0)]  = skvol[n.where(skvol > 0)]
#         sknvol[skvol > 0]  = skvol[skvol > 0]
    
#     # nlpx... removed - 0.5 * dpx - 0.5 * hpx in-line with new definition nlpx
#     dist_to_nec_cen = nlpx / 2 + nec_in_den
#     nec_cen = n.array(neck_ori - unitvec( neck_ori - head_pos ) \
#                     * dist_to_nec_cen).astype(int)
    
#     vol = n.asarray(vol | nvol).astype(int)
    
#     if need_skl:
#         sknvol[head_pos[0], head_pos[1], head_pos[2]] = head_diam
#     else:
#         sknvol = None

#     return vol, neck_ori, nec_cen, head_pos, sknvol

# def spine_image(sp_mdl = None, pxr = None, fwhms = None):

#     if pxr == None:
#         raise ValueError('pxr must be defined')

#     if fwhms == None:
#         raise ValueError('fwhms must be defined')
    
#     if sp_mdl == None:
#         sp_mdl = spine_model( px_res = pxr )

#     fwhms = n.asarray(fwhms)

#     (fx,fy,fz) = tuple(fwhms / pxr)
#     psf = gauss3d(fx, fy, fz, 1e-3)
#     #print "fftconvolve shapes:", sp_mdl.shape, psf.shape
    
#     spine_image = fftconvolve( n.cast[float](sp_mdl), \
#                                n.cast[float](psf), mode='same')
#     return spine_image, psf

# def spine_model_new( px_res = None, \
#                      pars = {}, \
#                      vol_size = (200,200,160), \
#                      vol_inds = None, \
#                      need_skl = False):

#     # by definition dendrite runs along y axis
#     # (and if applicable, low end -y, high end +y)

#     pars = dparse({ 'origin_x'     : 150, \
#                     'origin_y'     : 100, \
#                     'origin_z'     : 80, \
#                     'angle_den'    : 0, \
#                     'dend_diam'    : 1, \
#                     'neck_diam'    : 0.5, \
#                     'neck_length'  : 0.5, \
#                     'head_diam'    : 1. , \
#                     'angle_roundz' : 0, \
#                     'angle_toz'    : 0
#                     }, pars)

#     origin_x = pars['origin_x']
#     origin_y = pars['origin_y']
#     origin_z = pars['origin_z']
#     angle_den = pars['angle_den']
#     dend_diam = pars['dend_diam']
#     neck_diam = pars['neck_diam']
#     neck_length = pars['neck_length']
#     head_diam = pars['head_diam']
#     angle_roundz = pars['angle_roundz']
#     angle_toz = pars['angle_toz']

#     #pdbg("debug verbose", "Parameters", pars)
    
#     if px_res == None:
#         raise ValueError('Pixel resolution must be defined.')

#     if vol_inds == None:
#         print 'Calculating indices in model:'
#         vol_inds = n.indices(vol_size)

#     dpx = dend_diam / px_res
#     #print "dpx", dpx, dpx/2.
#     nwpx = neck_diam / px_res
#     hpx = head_diam / px_res
#     #print "hpx", hpx, hpx/2.

#     nlpx = neck_length / px_res
#     #print "nlpx", nlpx
#     # dlpx = dend_len / px_res

#     dend_vec = (0, n.cos(angle_den),n.sin(angle_den))
#     # print "dend_vec", dend_vec

#     nec_vec_dir = n.array( (n.cos( angle_roundz ) * n.cos( angle_toz ), \
#                             # x component
#                             n.sin( angle_roundz ) * n.cos( angle_toz ), \
#                             # y component
#                             n.sin( angle_toz )
#                             # z component
#                             ) )
    
#     #print "Nec vec dir length:", n.linalg.norm(nec_vec_dir)
#     #print "Dend vec length:", n.linalg.norm(dend_vec)
#     #nec_in_den = (dpx/2.) / n.sin( n.arccos( n.dot(nec_vec_dir, dend_vec) ) )
#     nec_in_den = (dpx/2.) / n.sqrt( 1 - n.dot(nec_vec_dir, dend_vec)**2 )
#     #print "Dendrite radius:", dpx/2.
#     #print "Nec length in dendrite:", nec_in_den

#     nec_vec_len = (nlpx + nec_in_den + hpx/2.)
#     nec_vec = nec_vec_len * nec_vec_dir
#     #print "nec_vec", nec_vec_len, nec_vec_dir
    
#     sz = n.asarray( vol_size ).round().astype(int)

#     neck_ori = n.array( ( origin_x, origin_y, origin_z ) )
#     #print "neck ori:", neck_ori

#     head_pos = n.asarray(neck_ori - nec_vec) #.round(0).astype(int)
#     #print "head_pos:", head_pos

#     # make dendrite - by definition runs along y axis (=> x = 0)
#     #print sz, neck_ori, dend_vec, n.array(dpx / 2.).round()
#     vol = cylinder( sz, neck_ori, dend_vec, n.array(dpx / 2.).round(), \
#                     inds = vol_inds)
    
#     if need_skl:
#         skvol = cylinder( sz, neck_ori, dend_vec, 1, \
#                           inds = vol_inds).astype(float32) \
#                 * dend_diam
        
#     # make head
#     vol |= sphere( sz, head_pos, n.asarray( hpx / 2.).round(), \
#                    inds = vol_inds)
# #    vol = sphere( sz, head_pos, n.asarray( hpx / 2.).round() )

#     # make neck
#     nvol = cylinder( sz, neck_ori, unitvec( nec_vec ), \
#                      n.asarray(nwpx / 2.).round(), inds = vol_inds )
#     nvol &= clip_plane( sz, neck_ori, -unitvec( nec_vec ), \
#                         inds = vol_inds)
# #    nvol = clip_plane( sz, neck_ori, -unitvec( nec_vec ) )
#     nvol &= clip_plane( sz, head_pos, unitvec( nec_vec ), \
#                         inds = vol_inds)

#     if need_skl:
#         sknvol = cylinder( sz, neck_ori, nec_vec_dir, 1, \
#                            inds = vol_inds);
#         sknvol &= clip_plane( sz, neck_ori, -nec_vec_dir, \
#                               inds = vol_inds)
#         sknvol &= clip_plane( sz, head_pos, nec_vec_dir, \
#                               inds = vol_inds)
#         sknvol = n.cast[float32](sknvol)
#         sknvol *= neck_diam
#         # sknvol[n.where(skvol > 0)]  = skvol[n.where(skvol > 0)]
#         sknvol[skvol > 0]  = skvol[skvol > 0]
#         del skvol
    
#     dist_to_nec_cen = nlpx / 2 + nec_in_den

#     #print 'case 1:', unitvec( neck_ori - head_pos )
#     #print 'case 2:', unitvec( nec_vec )
    
#     nec_cen = n.asarray(neck_ori - unitvec( neck_ori - head_pos ) \
#                         * dist_to_nec_cen).astype(int)

#     vol |= nvol
#     del nvol
#     vol = n.cast[int](vol)
    
#     if need_skl:
#         sknvol[head_pos[0], head_pos[1], head_pos[2]] = head_diam
#     else:
#         sknvol = None

#     return vol, neck_ori, nec_cen, head_pos, sknvol

# ----------------------------------------------------------------------------
