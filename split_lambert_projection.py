from matplotlib.axes import Axes
from matplotlib import cbook, docstring
from matplotlib.patches import Patch, cbook, transforms, artist
from matplotlib.path import Path
from matplotlib.ticker import Formatter, Locator, NullLocator, \
    FixedLocator, NullFormatter
from matplotlib.transforms import Affine2D, Affine2DBase, Bbox, \
    BboxTransformTo, IdentityTransform, Transform, TransformWrapper
from matplotlib.projections import register_projection, PolarAxes
from matplotlib.projections.geo import GeoAxes, LambertAxes
import matplotlib.spines as mspines
import matplotlib.axis as maxis
from copy import copy
from split_lambert_transforms import SplitLambertTransform

import numpy as np

class TwoCircle(Patch):
    """
    A scale-free ellipse.
    """
    def __str__(self):
        return "TwoCircle(%s, %s; %s)" % (self.center[0], self.center[1],
                                          self.radius)

    @docstring.dedent_interpd
    def __init__(self, xy, radius, **kwargs):
        """
        xy : array_like
          center of two circles
        radius : scalar
          size of each circle
          
        Valid kwargs are:
        %(Patch)s
        """
        Patch.__init__(self, **kwargs)
        self.center = xy
        self.radius = radius
        self.width = 4. # two x unit circle (i.e. from +1 to -1)
        self.height = 2. # one x unit circle
        print "h", self.height
        path = copy(Path.unit_circle())
        n_pts = path.vertices.shape[0]
        path.vertices = np.tile(path.vertices, [2,1])
        path.vertices[:n_pts,0] -= 1
        path.vertices[n_pts:,0] += 1
        path.codes = np.tile(path.codes, [2])
        self._path = path
        # Note: This cannot be calculated until this is added to an Axes
        self._patch_transform = transforms.IdentityTransform()

    def _recompute_transform(self):
        """NOTE: This cannot be called until after this has been added
                 to an Axes, otherwise unit conversion will fail. This
                 makes it very important to call the accessor method and
                 not directly access the transformation member variable.
        """
        center = (self.convert_xunits(self.center[0]),
                  self.convert_yunits(self.center[1]))
        width = self.convert_xunits(self.width)
        height = self.convert_yunits(self.height)
        self._patch_transform = transforms.Affine2D() \
            .scale(1. / width, 1. / height) \
            .translate(*center)

    def get_path(self):
        """
        Return the vertices of the rectangle
        """
        return self._path

    def get_patch_transform(self):
        self._recompute_transform()
        return self._patch_transform

    def contains(self,ev):
        if ev.x is None or ev.y is None: return False,{}
        x, y = self.get_transform().inverted().transform_point((ev.x, ev.y))
        return (((x - 1) * (x - 1) + y * y) <= 1) or \
            (((x + 1) * (x + 1) + y * y) <= 1), {}
    
def elliptical_spine(cls, axes, center, xrad, yrad, **kwargs):
    '''
    (staticmethod) Returns an elliptical :class: `Spine`.
    '''
    path = Path.unit_circle()
    spine_type = 'circle' # since aspect = 2 and axes rect is 0-1, 0-1
    result = cls(axes, spine_type, path, **kwargs)
    set_patch_ellipse(result, center, xrad, yrad)
    return result

def set_patch_ellipse(cls, center, xrad, yrad):
    cls._patch_type = 'circle'
    cls._center = center
    cls._width = xrad * 2
    cls._height = yrad * 2
    cls._angle = 0
    cls.set_transform(cls.axes.transAxes)

class SplitLambertAxes(Axes):
    """
    A custom class for the split Lambert projection, an equal-area map
    projection using one circle for each hemisphere.
    """
    name = 'split_lambert'
    resolution = 75

    def __init__(self, *args, **kwargs):
        Axes.__init__(self, *args, **kwargs)
        self.set_aspect(0.5, adjustable='box', anchor='C')
        self.cla()

    def _init_axis(self):
        self.xaxis = maxis.XAxis(self) # xaxis == theta == latitude
        self.yaxis = maxis.YAxis(self) # yaxis == phi == longitude
        
        # Do not register xaxis or yaxis with spines -- as done in
        # Axes._init_axis() -- until SplitLambertAxes.xaxis.cla() works.
        # self.spines['split_lambert'].register_axis(self.yaxis)
        self._update_transScale()

    def cla(self):
        """
        Override to set up some reasonable defaults.
        """
        # Don't forget to call the base class
        Axes.cla(self)

        # Set up a default grid spacing
        self.set_theta_grid(np.pi/4.)
        self.set_phi_grid(np.pi/2.)
        self.set_phi_grid_ends(np.pi/8.)

        # Turn off minor ticking altogether
        self.xaxis.set_minor_locator(NullLocator())
        self.yaxis.set_minor_locator(NullLocator())

        # Do not display ticks
        self.xaxis.set_ticks_position('none')
        self.yaxis.set_ticks_position('none')

        # The limits on this projection are fixed -- they are not to
        # be changed by the user.  This makes the math in the
        # transformation itself easier, and since this is a toy
        # example, the easier, the better.
        Axes.set_xlim(self, 0, np.pi)
        Axes.set_ylim(self, 0, 2 * np.pi)

    def _set_lim_and_transforms(self):
        """
        This is called once when the plot is created to set up all the
        transforms for the data, text and grids.
        """
        # There are three important coordinate spaces going on here:
        #
        #    1. Data space: The space of the data itself
        #
        #    2. Axes space: The unit rectangle (0, 0) to (1, 1)
        #       covering the entire plot area.
        #
        #    3. Display space: The coordinates of the resulting image,
        #       often in pixels or dpi/inch.

        # This function makes heavy use of the Transform classes in
        # ``lib/matplotlib/transforms.py.`` For more information, see
        # the inline documentation there.

        # The goal of the first two transformations is to get from the
        # data space (in this case longitude and latitude) to axes
        # space.  It is separated into a non-affine and affine part so
        # that the non-affine part does not have to be recomputed when
        # a simple affine change to the figure has been made (such as
        # resizing the window or changing the dpi).

        # 1) The core transformation from data space into
        # rectilinear space defined in the SplitLambertTransform class.
        # Needs to transform theta, phi into x,y on a unit rectangle
        # (i.e. into axis co-ordinates).
        self.transProjection = SplitLambertTransform(self.resolution)

        # 2) The above has an output range that is not in the unit
        # rectangle, so scale and translate it so it fits correctly
        # within the axes. Need to scale height up by 2 and move origin.
        self.transAffine = Affine2D().scale(1,2).translate(0.5, 0.5)

        # 3) This is the transformation from axes space to display
        # space.
        self.transAxes = BboxTransformTo(self.bbox)

        # Now put these 3 transforms together -- from data all the way
        # to display coordinates.  Using the '+' operator, these
        # transforms will be applied "in order".  The transforms are
        # automatically simplified, if possible, by the underlying
        # transformation framework.
        self.transData = \
            self.transProjection + \
            self.transAffine + \
            self.transAxes
        
        self._xaxis_pretransform = \
            Affine2D() \
            .scale(1.0, 2 * np.pi) \
            .translate(0.0, np.pi/4.)
        self._xaxis_transform = \
            self._xaxis_pretransform + \
            self.transData

        self._yaxis_pretransform = \
            Affine2D().scale(np.pi/2. - 1e-8, 1.0).translate(0.0, 0.)           
        self._yaxis_transform = \
            self._yaxis_pretransform + \
            self.transData

    def get_xaxis_transform(self, which='grid'):
        """
        Override this method to provide a transformation for the
        x-axis grid and ticks.
        """
        assert which in ['tick1','tick2','grid']
        return self._xaxis_transform

    def get_yaxis_transform(self,which='grid'):
        """
        Override this method to provide a transformation for the
        y-axis grid and ticks.
        """
        assert which in ['tick1','tick2','grid']
        return self._yaxis_transform
    
    def _gen_axes_patch(self):
        """
        Override this method to define the shape that is used for the
        background of the plot.  It should be a subclass of Patch.
        """
        return TwoCircle((0.5, 0.5), 0.25)

    def _gen_axes_spines(self):
        # circular_spine(axes, center, radius)
        return {'top' : elliptical_spine(mspines.Spine, \
                self, (0.25, 0.5), 0.25, 0.5),
                'bottom' : elliptical_spine(mspines.Spine, \
                self, (0.75, 0.5), 0.25, 0.5)}
    
    # Prevent the user from applying scales to one or both of the
    # axes.  In this particular case, scaling the axes wouldn't make
    # sense, so we don't allow it.
    def set_xscale(self, *args, **kwargs):
        if args[0] != 'linear':
            raise NotImplementedError
        Axes.set_xscale(self, *args, **kwargs)

    def set_yscale(self, *args, **kwargs):
        if args[0] != 'linear':
            raise NotImplementedError
        Axes.set_yscale(self, *args, **kwargs)

    # Prevent the user from changing the axes limits.  In our case, we
    # want to display the whole sphere all the time, so we override
    # set_xlim and set_ylim to ignore any input.  This also applies to
    # interactive panning and zooming in the GUI interfaces.
    def set_xlim(self, *args, **kwargs):
        Axes.set_xlim(self, 0, np.pi)
        Axes.set_ylim(self, 0, 2.0 * np.pi)
    set_ylim = set_xlim

    def format_coord(self, theta, phi):
        """
        Override this method to change how the values are displayed in
        the status bar.
        """
        # U03D1 = theta, U03D5 = phi
        return u'\u03D1 %f, \u03D5 %f' % (theta, phi)

    class RadianFormatter(Formatter):
        """
        This is a custom formatter. Returns values as factors of pi.
        """
        def __init__(self, round_to=0.01):
            self._round_to = round_to

        def __call__(self, x, pos=None):
            pix = x / np.pi
            n_dp = len(str(self._round_to).split('.')[-1])
            rads = round(pix / self._round_to) * self._round_to
            # U03C0 = lowercase pi
            format_str = "%0." + "%0d" % n_dp + "f" + u'\u03C0'
            return format_str % rads

    def set_theta_grid(self, rads):
        """
        Set the number of radians between each theta grid,
        i.e. circle of latitude

        This is an example method that is specific to this projection
        class -- it provides a more convenient interface to set the
        ticking than set_xticks would.
        """
        # Set up a FixedLocator at each of the points, evenly spaced
        # by radians.
        number = 0 #(np.pi / rads) + 1
        self.xaxis.set_major_locator(
            FixedLocator(
                np.linspace(0, np.pi, number, True)[1:-1]))
        self.xaxis.set_major_formatter(self.RadianFormatter(0.01))

    def set_phi_grid(self, rads):
        """
        Set the number of radians between each phi grid.

        This is an example method that is specific to this projection
        class -- it provides a more convenient interface than
        set_yticks would.
        """
        # Set up a FixedLocator at each of the points, evenly spaced
        # by radians.
        number = 0 # (2 * np.pi / rads) + 1
        self.yaxis.set_major_locator(
            FixedLocator(
                np.linspace(0, 2 * np.pi, number, True)[1:-1]))
        self.yaxis.set_major_formatter(self.RadianFormatter(0.01))

    def set_phi_grid_ends(self, rads):
        """
        Set the theta(s) at which to stop drawing the phi grids.

        Often, in geographic projections, you wouldn't want to draw
        longitude gridlines near the poles.  This allows the user to
        specify the location at which to stop drawing longitude grids.

        This is an example method that is specific to this projection
        class -- it provides an interface to something that has no
        analogy in the base Axes class.
        """
        phi_cap = rads #* np.pi
        self._yaxis_pretransform \
            .clear() \
            .scale(np.pi/2. - rads - 1e-8, 1.0) \
            .translate(rads, 0.)
                    #phi_cap * 2.0) \


    def get_data_ratio(self):
        """
        Return the aspect ratio of the data itself.

        This method should be overridden by any Axes that have a
        fixed data ratio.
        """
        return 1.0

    # Interactive panning and zooming is not supported with this projection,
    # so we override all of the following methods to disable it.
    def can_zoom(self):
        """
        Return True if this axes support the zoom box
        """
        return False
    def start_pan(self, x, y, button):
        pass
    def end_pan(self):
        pass
    def drag_pan(self, button, key, x, y):
        pass    
    
# Now register the projection with matplotlib so the user can select
# it.
register_projection(SplitLambertAxes)

# Now make a simple example using the custom projection.

# import matplotlib.pyplot as plt
# ax = plt.subplot(111, projection="split_lambert")
# H = np.pi
# #p = plt.plot([0, 3/8.*H, 3/8.*H, 11/8.*H],
# #             [0,      0,   H/6.,    H/6.], "o--")
# p = plt.plot([0, H], [0, 0], 'o-')

# plt.show()
