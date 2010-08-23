import numpy as np
from matplotlib.projections import register_projection, PolarAxes
from matplotlib.transforms import Affine2D, Affine2DBase, Bbox, \
    BboxTransformTo, IdentityTransform, Transform, TransformWrapper
from matplotlib.path import Path

class SplitLambertTransform(Transform):
    """
    The base Lambert transform.
    
    Converts theta, phi 3-d direction data co-ordinates to x, y into
    (2,1) rectangle axis co-ordinates.
    """
    input_dims = 2
    output_dims = 2
    is_separable = False

    def __init__(self, resolution):
        '''Create a new Split Lambert transform.
        
        Resolution is the number of steps to interpolate between each
        input line segment to approximate its path in curved Lambert space.
        '''
        Transform.__init__(self)
        self._resolution = resolution
    
    def _fix_theta(self, theta):
        return np.abs(np.pi - np.abs(np.pi - theta % (2 * np.pi)))
    
    def _fix_phi(self, phi):
        return phi % (2 * np.pi)
    
    def transform(self, thph):
        """
        Override the transform method to implement the custom transform.
        
        The input and output are Nx2 numpy arrays.
        
        Corrects all theta values to 0 < pi, and phi inputs to 0 < phi < 2 * pi.
        
        Theta is angle from (0,0,1) i.e. up axis,
        and phi is angle from (0,1,0) i.e. right.
        """
        # get values
        theta = self._fix_theta(thph[:, 0:1])
        phi = self._fix_phi(thph[:,1:2])
        
        # sort into flipped and not
        top_hemi = np.atleast_1d((theta < (np.pi / 2.)).squeeze())
        
        # flip theta and phi values for top hemisphere
        theta[top_hemi] = np.pi - theta[top_hemi]
        phi[top_hemi] = (np.pi - phi[top_hemi]) % (2 * np.pi)
        
        # 3d -> 2d polar co-ordinates
        rho = 2 * np.sin((np.pi - theta)/2.)
        rho /= 4 * np.sqrt(2.)
        psi = phi
        tr = np.concatenate((psi, rho), 1)

        # 2d polar -> 2d cartesian co-ordinates
        polarTrans = PolarAxes.PolarTransform()
        xy = polarTrans.transform(tr)

        xy[top_hemi, 0] -= 0.25
        xy[~top_hemi, 0] += 0.25
        return xy

        # This is where things get interesting.  With this projection,
        # straight lines in data space become curves in display space.
        # This is done by interpolating new values between the input
        # values of the data.  Since ``transform`` must not return a
        # differently-sized array, any transform that requires
        # changing the length of the data array must happen within
        # ``transform_path``.
    def transform_path(self, path):
        ipath = path.interpolated(self._resolution)
        nvert = self.transform(ipath.vertices)        
        ncode = np.array([y for x,y in ipath.iter_segments(simplify=False)])
        #npath = Path(nvert, ncode)
        # must change all codes in that segment to Path.MOVETO
        #print '\n'.join([str(x) + str(y) for 
        #                 x,y in npath.iter_segments(simplify=False)])
        # could try making line segment where x switches sign a MOVETO
        # currently will depend on there being no intermediate points
        # in that line segment
        nx = nvert[:,0]
        flip = [x + 1 for x in np.nonzero(np.diff(nx > 0))]
        ncode[flip] = Path.MOVETO
        npath = Path(nvert, ncode)
        
        # print [y for x,y in npath.iter_segments(simplify=False)]
        return npath

    transform_path_non_affine = transform_path

    def inverted(self):
        return InvertedSplitLambertTransform()
    inverted.__doc__ = Transform.inverted.__doc__

class InvertedSplitLambertTransform(Transform):
    input_dims = 2
    output_dims = 2
    is_separable = False
    
    def transform(self, xy):
        x = xy[:, 0:1]
        y = xy[:, 1:2]
        
        # work out which are in the top hemisphere
        top_hemi = np.atleast_1d((x < 0.).squeeze())
        # then collapse them all into one circle
        xy[top_hemi, 0] += 0.25
        xy[~top_hemi, 0] -= 0.25

        # convert to polar co-ordinates
        polarInvTrans = PolarAxes.InvertedPolarTransform()
        
        tr = polarInvTrans.transform(xy)
        # replace psi nans
        tr[:,0:1][np.isnan(tr[:,0:1])] = 0.

        # convert 2d polar -> 3d polar
        rho = tr[:,1:2]
        psi = tr[:,0:1]
        theta = np.pi - 2 * np.arcsin(2 * np.sqrt(2.) * rho)
        phi = psi
        
        # flip theta and phi values for top hemisphere
        theta[top_hemi] = np.pi - theta[top_hemi]
        phi[top_hemi] = (np.pi - phi[top_hemi]) % (2 * np.pi)
            
        return np.concatenate((theta, phi), 1)
    transform.__doc__ = Transform.transform.__doc__

    def inverted(self):
        # The inverse of the inverse is the original transform... ;)
        return SplitLambertTransform()
    inverted.__doc__ = Transform.inverted.__doc__

    
def test_split_lambert_transform():
    slt = SplitLambertTransform()
    tp = np.array([[0, 0],
                   [np.pi, 0],
                   [np.pi/2., 0.],
                   [np.pi/2., np.pi]])
    answers = np.array([[-0.25, 0],
                        [0.25, 0],
                        [0.5, 0],
                        [0, 0]])
    assert(np.allclose(slt.transform(tp), answers))
    islt = slt.inverted()
    print "********"
    assert(np.allclose(islt.transform(answers), tp))
    
# def mprint(a, a_, b, b_):
#     assert a_.shape[0] == b_.shape[0]
#     afield = '%s' % ('%' + '%d' % (len(a) + 1) + 's ')
#     bfield = '%s' % ('%' + '%d' % (len(b) + 1) + 's ')
#     n_lines = a_.shape[0]
#     for i in xrange(n_lines):
#         if i == 0:
#             a__ = a
#             b__ = b
#         else:
#             a__ = ''
#             b__ = ''
#         print afield % a__ + ' '.join(['%5.2f' % x for x in a_[i]]) + '  ' + \
#             bfield % b__ + ' '.join(['%5.2f' % y for y in b_[i]]) 

        
