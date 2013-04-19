import numpy as np
import matplotlib
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
    
    def transform_non_affine(self, thph):
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
        rho = (2 * np.sin((np.pi - theta)/2.))
        rho /= 4 * np.sqrt(2.)
        psi = phi % (2 * np.pi)
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
    def transform_path_non_affine(self, path):
        ipath = self.interpolate_path(path) # interpolate in data-space
        
        nvert = self.transform(ipath.vertices) # transform interpolated pts
        ncode = np.array([y for x,y in ipath.iter_segments(simplify=False)])
        nx = nvert[:,0]
        flip = [x + 1 for x in np.nonzero(np.diff(nx > 0))]
        ncode[flip] = Path.MOVETO
        npath = Path(nvert, ncode)
        return npath
    transform_path_non_affine.__doc__ = \
        Transform.transform_path_non_affine.__doc__

    def interpolate_path(self, path):
        """
        Returns a new path resampled to length N x steps.  Does not
        currently handle interpolating curves.
        """
        steps = self._resolution
        if steps == 1:
            return self
        
        vertices = circular_interpolation(path.vertices, steps)
        codes = path.codes
        if codes is not None:
            new_codes = Path.LINETO * np.ones(((len(codes) - 1) * steps + 1, ))
            new_codes[0::steps] = codes
        else:
            new_codes = None
        return Path(vertices, new_codes)

    if matplotlib.__version__ < '1.2':
        transform = transform_non_affine
        transform_path = transform_path_non_affine
        transform_path.__doc__ = Transform.transform_path.__doc__

    def inverted(self):
        return InvertedSplitLambertTransform()
    inverted.__doc__ = Transform.inverted.__doc__

def circular_interpolation(a, steps):
    '''
    special case for interpolation of circular co-ordinates
    i.e. phi has to be done `mod 2 * pi`
    '''
    if steps == 1:
        return a

    steps = np.floor(steps).astype(int)
    
    new_length = ((len(a) - 1) * steps) + 1
    new_shape = list(a.shape)
    new_shape[0] = new_length
    result = np.zeros(new_shape, a.dtype)

    result[0] = a[0]
    a0 = a[0:-1]
    a1 = a[1:  ]
    # problem is:
    # * only with phi dimension
    # * when shortest distance between a0 and a1 goes through
    #   2n * np.pi (for n E N), rather than just difference
       
    delta = (a1 - a0)
    ow = np.abs(delta[:,1]) > np.pi # go other way on these ones
    delta[ow,1] = (2 * np.pi - np.abs(delta[ow,1])) * np.sign(delta[ow,1]) * -1
    delta /= steps
    
    for i in range(1, int(steps)):
        result[i::steps] = delta * i + a0
    result[steps::steps] = a1

    return result
    
class InvertedSplitLambertTransform(Transform):
    input_dims = 2
    output_dims = 2
    is_separable = False
    
    def transform_non_affine(self, xy):
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
    transform_non_affine.__doc__ = Transform.transform_non_affine.__doc__

    if matplotlib.__version__ < '1.2':
        transform = transform_non_affine

    def inverted(self):
        # The inverse of the inverse is the original transform... ;)
        return SplitLambertTransform()
    inverted.__doc__ = Transform.inverted.__doc__
    
def test_split_lambert_transform():
    slt = SplitLambertTransform(100.)
    tp = np.array([[np.pi, 0],
                   [np.pi/2., 0.],
                   [np.pi/2., np.pi]])
    answers = np.array([[0.25, 0],
                        [0.5, 0],
                        [0, 0]])
    np.testing.assert_array_almost_equal(slt.transform(tp), answers)
    islt = slt.inverted()
    np.testing.assert_array_almost_equal(islt.transform(answers), tp)
    
