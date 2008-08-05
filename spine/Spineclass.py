import numpy as n, ncdf
import scipy.optimize as opt
from scipy.ndimage import rotate
from gaussian import gauss3d
from myutils import printif
import myutils
from vectors import angle_between, rotate_about_centre
from spine.extract import *
import spine.fit, spine.model
from spine.utils import model_pxres, read_params
from copy import copy
import os.path

posmods = ['pynetcdf', 'Scientific.IO.NetCDF', 'Scientific_NetCDF']
tried = 0
for mod in posmods:
    try:
        exec("from %s import NetCDFFile" % mod)
        break
    except ImportError:
        tried += 1
if tried == len(posmods):
    raise ImportError("No module containing NetCDFFile found.")

graphic = True
try:
    import pylab
except:
    graphic = False

try:
    from enthought.mayavi.tools import mlab
    from enthought.mayavi.sources.array_source import ArraySource
    from enthought.mayavi.modules.api import IsoSurface, Outline, Glyph
    from sources_scatter_source import ScatterSource
    enth = True
except ImportError:
    enth = False

stage_initd    = 1001
stage_userptd  = 1002
stage_measured = 1003
stage_fitd     = 1004

fit_thresh = 0.7

class Spine:
    def __init__(self, *args, **kwargs):
        '''sp = Spine(spine_name,
           pixel_spacing,
           psf_fwhms,
           source_data=None,
           user_pts=None)
or
sp = Spine(filename)

Creates a spine object incorporating image data, and methods to reconstruct
sub-resolution information.

:Parameters:
    pixel_spacing : (3,) tuple
        pixel spacing of volume data in microns (x,y,z)
    psf_fwhms     : (3,) tuple 
        fwhms of psf in microns (x,y,z)
    source_data   : (i,j,k) ndarray
        volume data containing spine
    user_pts      : list of 3x(3,) tuples
        co-ordinates of three key points in spine--
        (1) spine head,
        (2) neck-dendrite intersection,
        (3) 2nd pt on dendrite

    filename      : string
        path to spine ncdf-saved file

:Returns:
    spine         : Spine instance

Examples
--------
After loading a spine, the usual progression of use is:
  sp.pick_by_click(filename, pv, varname="var0")
or
  sp.pick_by_coords(filename, coord, varname="var0")
then
  sp.click_user_pts(pv)
or
  sp.load_user_pts(coords_list)
then
  sp.optimize()

At any time, sp.write_ncdf() can be called to save the current state into a
netcdf file.
'''
        if len(args) == 3:
            # first usage form
            if kwargs.has_key('source_data'):
                source_data = kwargs['source_data']
            else:
                source_data = n.zeros((1,1,1), dtype=n.integer)
            if kwargs.has_key('user_pts'):
                self.load_user_pts(kwargs['user_pts'])
            self.picked = False
            self.stage = stage_initd
            self.user_pts = []
            self.params = {}
            self.big_offset = (0,0,0)
            self.big_name = ""
            self.source_data = n.asarray(source_data)
            self.pxspc = n.array( args[1] )
            self.fwhm = n.array( args[2] )
            self.name = args[0]
            print("Pre-check name: %s" % self.name)
            unique = myutils.find_unique_filename(self.name + '.txt')
            self.name = os.path.splitext(unique)[0]
            print("Post-check name: %s" % self.name)
            self.rotated_data = None
        elif len(args) == 1:
            # second usage form - read from file
            fname = args[0]
            ncfile = NetCDFFile( fname, 'r' )
            var = ncfile.variables['source_data']
            data = var[:]
            if NetCDFFile.__module__ == 'Scientific.IO.NetCDF':
                if data.typecode == 'd':
                    data = data.astype('uint8')
            else:
                if data.dtype.kind == 'i':
                    data = data.astype('uint8')
            self.source_data = n.asarray(data)
            self.name = ncfile.name
            self.pxspc = n.asarray(ncfile.pxspc)
            self.fwhm = n.asarray(ncfile.fwhm)
            self.stage = int(ncfile.stage)
            self.big_name = ncfile.big_name
            self.big_offset = n.asarray(ncfile.big_offset)
            self.user_pts = []
            self.params = {}
            if self.stage >= stage_userptd:
                self.user_pts.append( ncfile.headpt )
                #self.user_pts.append( ncfile.neckpt )
                self.user_pts.append( ncfile.origpt )
                self.user_pts.append( ncfile.dendpt )
            if self.stage >= stage_measured:
                self.params['angle_toz'] = float( ncfile.angle_toz )
                self.params['angle_roundz'] = float( ncfile.angle_roundz )
                self.params['neck_length'] = float( ncfile.neck_length )
                self.params['dend_diam'] = float( ncfile.dend_diam )
                self.params['head_int'] = float( ncfile.head_int )
                self.params['neck_int'] = float( ncfile.neck_int )
                try:
                    self.params['angle_den'] = float( ncfile.angle_den )
                except AttributeError:
                    #print "No dendrite z angle found. Defaulting to zero."
                    self.params['angle_den'] = 0.
                try:
                    self.params['orient_den'] = float( ncfile.angle_den )
                except AttributeError:
                    #print "No dendrite x-y angle found. Defaulting to zero."
                    self.params['orient_den'] = 0.
            if self.stage >= stage_fitd:
                self.params['head_diam'] = float( ncfile.head_diam )
                self.params['neck_diam'] = float( ncfile.neck_diam )
            ncfile.close()

    def write_ncdf(self):
        '''Writes a ncdf file containing the volume.
        '''
        fname = self.name + '.nc'
        ncfile = NetCDFFile( fname, 'w' )
        dimslist = ['datax','datay','dataz']
        sd_shape = self.source_data.shape
        dimssize = [sd_shape[0], sd_shape[1], sd_shape[2]]
        for i in xrange( len(dimslist) ):
            ncfile.createDimension( dimslist[i], self.source_data.shape[i] )
        # store source_data - may be empty - but okay
        var = ncfile.createVariable( 'source_data', 'd', tuple( dimslist ) )
        var[:] = self.source_data
        # may also need to store rotated_data sometime??
        
        # store attributes
        # defined by __init__ for sure
        ncfile.name = self.name
        ncfile.pxspc = self.pxspc
        ncfile.fwhm = self.fwhm
        ncfile.stage = self.stage
        ncfile.big_name = self.big_name
        ncfile.big_offset = self.big_offset
        # defined by collect_user_pts
        if self.stage >= stage_userptd:
            ncfile.headpt = self.user_pts[0]
            ncfile.origpt = self.user_pts[1]
            ncfile.dendpt = self.user_pts[2]
        # defined by calc_params
        if self.stage >= stage_measured:
            ncfile.angle_toz = self.params["angle_toz"]
            ncfile.angle_roundz = self.params["angle_roundz"]
            ncfile.neck_length = self.params["neck_length"]
            ncfile.dend_diam = self.params["dend_diam"]
            ncfile.head_int = self.params["head_int"]
            ncfile.neck_int = self.params["neck_int"]
            ncfile.angle_den = self.params["angle_den"]
            ncfile.orient_den = self.params["orient_den"]
        # defined by optimize
        if self.stage >= stage_fitd:
            ncfile.head_diam = self.params["head_diam"]
            ncfile.neck_diam = self.params["neck_diam"]
        ncfile.close()

    def click_user_pts(self, pv):
        '''Prompts user to select spine key points.
        Calls collect_user_pts_cb.

        Usage:
          click_user_pts(pv)'''
        self.user_pts = []
        pv.AddImg(self.source_data, self.name)
        pv.DisconnectClickCallback()
        pv.ConnectClickCallback(self.click_user_pts_cb, pv)
        print "Click at spine head...",

    def save_params(self, setname="Construction"):
        fname = self.name + '.txt'
        if os.path.exists(fname):
            f = open(fname, 'a')
        else:
            f = open(fname, 'w')
        f.write("%s\n" % setname)
        for k,v in self.params.iteritems():
            f.write("%16s : %f\n" % (k,v))
        f.write("\n")
        f.close()

    def save_extra(self, name='', val=n.nan):
        fname = self.name + '.txt'
        if os.path.exists(fname):
            f = open(fname, 'a')
        else:
            f = open(fname, 'w')
        f.write("%16s : %f\n" % (name, val))
        f.write("\n")
        f.close()

    def click_user_pts_cb(self, x, y, z, *args):
        pv = args[0]
        if not x is None and not y is None:
            x = int( x )
            y = int( y )
            print "adding point:", x, y, z, "."
            self.user_pts.append( (x, y, z) )
            if len(self.user_pts) == 1:
                print "Click at spine origin...",
            elif len(self.user_pts) == 2:
                print "Click at second dendrite point...",
            elif len(self.user_pts) >= 3:
                print "All points collected."
                pv.DisconnectClickCallback()
                self.stage = stage_userptd
                self.calc_params()

    def load_user_pts(self, coordslist):
        '''Enter pre-determined coordinates of spine user pts.

        Usage:
          spine.load_user_pts(coordslist)
        where coords_list is a list of 4x (3,) coord tuples
        in the order head, neck, origin, dendrite'''
        if len(coordslist) == 3:
            self.user_pts = coordslist
        else:
            raise DeprecationWarning("Neck point is no longer required.")
            self.user_pts = [coordslist[0], coordslist[2], coordslist[3]]
        self.stage = stage_userptd
        self.calc_params()
        
    def calc_params(self, verbose=False):
        '''Calculates the initial parameters
          angle to z,
          angle round z,
          dendrite diameter,
          neck length
        from the user selected pts.
        '''
        # save original params
        self.save_params()
        # extract click_coords and correct for display scaling
        if len(self.user_pts) < 3:
            raise StandardError( "Not enough co-ordinates entered." )
        
        p = {}
        # extract points from data
        den  = n.asarray( self.user_pts[2] )
        ori  = n.asarray( self.user_pts[1] )
        head = n.asarray( self.user_pts[0] )

        # ---- handle dvec and rotation ----
        
        # should weight dendrite fitting by proximity to origin
        #fit_weights = spine.fit.create_fit_weights(self)
        dvec_a, dvec = spine.fit.fit_dendrite_vector(\
            self.source_data, threshold=fit_thresh)
        
        dvec = unitvec(dvec * self.pxspc)
        print "dvec", dvec    
        dvec_a *= self.pxspc
        pdvec = perpz( dvec )
        # arbitrarily select x axis if dendrite extends perfectly in z
        if n.any(n.isnan(pdvec)):
            pdvec = n.array([1.,0.,0.])

        spvec = (head - ori) * self.pxspc
        printif( verbose, "Spine vector:", spvec )
        # ^ vector from origin to head, corrected for anisotropy
        # by definition, spvec must extend in +ve x direction

        # calculates angle relative to +ve y axis
        p['orient_den'] = extract_orient_den(dvec, spvec, self.pxspc)
        #n.arctan( dvec[0] / float(dvec[1]) )

        print 'orient_den (corrected)', p['orient_den']
        deg = 180. / n.pi
        # x-y rotation of real data to align to models
        self.rotated_data = rotate( self.source_data, \
                                    p['orient_den'] * deg, \
                                    axes=(0,1), reshape=False)
        self.rotated_data /= self.rotated_data.max()
        # NB rotated_data is NOT currently stored in the ncdf files
        # rotation is about volume centre, making dvec_x,y = (0,1)
        
        self.dvecr = dvec.copy()
        self.dvecr[0:2] = rotate_about_centre( self.dvecr[0:2], (0,0), \
                                               -p['orient_den'] )
        self.dvecr_a = dvec_a.copy()
        #print "dvec_a", dvec_a
        centre = n.asarray(self.source_data.shape) * self.pxspc / 2.
        #print "centre", centre
        self.dvecr_a[0:2] = rotate_about_centre( self.dvecr_a[0:2], \
                                                 centre[0:2], -p['orient_den'] )

        spvecr = spvec.copy()
        spvecr[0:2] = rotate_about_centre( spvecr[0:2], (0,0), -p['orient_den'] )
        # ---- calc other params ----
        p['angle_den'] = sp_extract_dend_angle( self.dvecr, spvecr )
        print "Calc angle_den: %0.3f" % p['angle_den']
        p['dend_diam'] = sp_extract_dend_diam( den, pdvec, \
                                               self.source_data, self.pxspc,\
                                               self.fwhm[0] )
        p['head_diam'] = 0.5 # can't know it so take a wild guess
        
        p['neck_length'] = sp_extract_neck_length( \
            spvec, dvec, p['head_diam'], p['dend_diam'] )
        p['angle_roundz'] = sp_extract_angle_roundz( spvec, pdvec, \
                                                     p['orient_den'] )
        p['angle_toz'] = sp_extract_angle_toz( spvec )

        self.params = complete_params( p )
        self.stage = stage_measured
        # write results to file
        self.save_params('Calculated')
        #write_params( outfname, self.name, complete_params( p ) )

    def rotate_data(self):
        '''run in place of calc_params when params are specified rather than
        calculated'''
        # need to have orient_den defined before running
        deg = 180. / n.pi
        # x-y rotation of real data to align to models
        self.rotated_data = rotate( self.source_data, \
                                    -self.params['orient_den'] * deg, \
                                    axes=(0,1), reshape=True)
        self.rotated_data /= self.rotated_data.max()
        # NB this is NOT currently stored in the ncdf files

    def pprint(self):
        '''Print spine details - not complete'''
        if self.npts == 4:
            upts = "defined"
        else:
            upts = "undefined"
        if self.params.keys() != []:
            params = "defined"
        else:
            params = "undefined"
        stufftoprint = {"spine"         : self.name,              \
                        "from data"     : self.big_name,          \
                        "at location"   : self.big_offset,        \
                        "data size"     : self.source_data.shape, \
                        "pixel spacing" : self.pxspc,             \
                        "user points"   : upts,                   \
                        "parameters"    : params}
        for k,v in stufftoprint.iteritems():
            print "%15s %15s" % (k, v)

    def load_params_from_file(self, fname):
        '''Loads the spine parameters into spine instance from a text
        file as read/written by spine.utils.read/write_params'''
        self.params = read_params( fname, self.name )

    def optimize(self, single=True, method='Nelder-Mead', graphic=True, \
                 t0=None):
        ''': Parameters :
    t0 : float
        optional override of initial t parameter, given in model_pxres unit
        '''
        # check params are valid
        if n.isnan(self.params['head_diam']):
            self.params['head_diam'] = 0.5
        if n.isnan(self.params['neck_diam']):
            self.params['neck_diam'] = 0.25
        
        # need to get dvec, dvec_a again...
        #dvec_a, dvec = spine.fit.fit_dendrite_vector(self.rotated_data,
        #                                             threshold=fit_thresh)
        # not sure about the 0.5 - but makes it match with the threshold
        # in calc_params
        dvec_a, dvec = self.dvecr_a, self.dvecr
        print "Dendrite vector:", dvec, "\n & point:", dvec_a
        #dvec *= self.pxspc
        #dvec = unitvec(dvec)
        #dvec_a *= self.pxspc# / model_pxres

        # pre-calculate some constant stuff for use in the loops
        model_size = (n.asarray(self.rotated_data.shape) * \
                      self.pxspc / model_pxres).round().astype(int)

        arz0 = self.params['angle_roundz']
        atz0 = self.params['angle_toz']

        # fit diameters and neck length
        mid_pxspc = n.zeros(3) + model_pxres * 4
        mid_fwhms = mid_pxspc * 2
        mid_fwhms_px = mid_fwhms / model_pxres
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

        if graphic and enth:
            sc = mlab.figure()
        else:
            sc = None

        if t0 == None:
            # fit origin
            # manual iteration through possible values to find minimum resid
            # - need to identify limits to t because 0 may be off to one side
            backlim = n.where((dvec > 0), n.zeros((3,)), model_size * model_pxres)
            print "backlim", backlim
            forwardlim = n.where((dvec < 0), n.zeros((3,)), model_size * \
                                 model_pxres)
            print "forwardlim", forwardlim
            
            tstep = 0.25  # in multiples of dvec - i.e. microns
            tsteps = [0.] # start at origin, at least 
            nbacksteps = ((backlim - dvec_a) / \
                          (-dvec * tstep)).min().round().astype(int)
            backsteps = range(-1, -nbacksteps, -1)
            print "back steps", backsteps
            backsteps.reverse()
            tsteps = backsteps + tsteps
            nforwardsteps = ((forwardlim - dvec_a) / \
                             (dvec * tstep)).min().round().astype(int)
            print "back steps", nbacksteps, "forward steps", nforwardsteps
            forwardsteps = range(1, nforwardsteps)
            print "forward steps", forwardsteps
            tsteps += forwardsteps
            tsteps = n.asarray(tsteps) * tstep
            # need to convert tsteps, dvec and dvec_a from micron to model_pxres
            tsteps /= model_pxres
            print "t steps", tsteps
            print "Finding optimal t..."
            t0 = tsteps[0]
            p0 = (hd0, nd0, nl0, dd0, arz0, atz0, A0, tsteps[0])
            ermin = spine.fit.residual( p0, dvec_a, dvec, \
                                        model_size, blur1, blur2, \
                                        self, graphic, sc )
            for tx in tsteps:
                print "t ==", tx
                p0 = (hd0, nd0, nl0, dd0, arz0, atz0, A0, tx)
                err = spine.fit.residual( p0, dvec_a, dvec, \
                                          model_size, blur1, blur2, \
                                          self, graphic, sc )
                if err < ermin:
                    t0 = tx
                    ermin = err
                
        # actually do optimisation below here
        p0 = (hd0, nd0, nl0, dd0, arz0, atz0, A0, t0)
        if single:
            err = spine.fit.residual( p0, dvec_a, dvec, \
                                      model_size, blur1, blur2, \
                                      self, graphic, sc )
            return err
        else:
            if method =='Nelder-Mead':
                iteration = [0,]
                pf = opt.fmin(spine.fit.residual, p0, \
                              args=(dvec_a, dvec, \
                                    model_size, blur1, blur2, \
                                    self, graphic, sc, iteration), \
                              full_output=1, xtol=0.01)
                self.params['head_diam']    = pf[0][0]
                self.params['neck_diam']    = pf[0][1]
                self.params['neck_length']  = pf[0][2]
                self.params['dend_diam']    = pf[0][3]
                self.params['angle_roundz'] = pf[0][4]
                self.params['angle_toz']    = pf[0][5]
                self.stage = stage_fitd
                self.save_params('Fitted')
                self.save_extra('Residual', pf[1])
                return pf

    def pickbycoords(self, fname, coords, varname='var0', compat=True):
        big = ncdf.r(fname, compat=compat, varname=varname)
        self.big_name = fname
        self.big_offset = tuple( coords )
        boxsz = array( (30, 30, 10) )
        sp = extract_box( big, tuple( coords ), tuple( boxsz ) )
        self.source_data = sp
        self.picked = True

    def pickbyclick(self, fname, pv=None, varname='var0', compat=True):
        '''Opens netcdf file and displays in pyvis so that
        user can select a spine for further analysis.
        
        Usage:
        sp.pick_data(file_name, pv[, varname='var0'][, compat=True])
        '''
        if pv is None:
            raise NameError("pv needs to be assigned to an instance of pyvis.")

        self.big_name = fname
        big = ncdf.r(fname, compat=compat, varname=varname)
        pv.AddImg(big, fname)
        pv.ConnectClickCallback(self.pickbyclick_cb, big)

    def pick_data(self, source, location, **kwargs):
        '''Picks source data for spine

        Usage:
          pick_data( source, location, **kwargs )
        where:
          source = filename or ndarray
          location = (3,) co-ords (in list, tuple or ndarray)
                     OR pyvis instance for getting a click at location

        **kwargs - varname= for netCDF file
                 - compat= mode for netCDF read
                 - boxsize= size of box to extract
        '''
        if type(source) == type(''):
            # read data from file
            self.bigname = source
            if kwargs.has_key('compat'):
                compat = kwargs['compat']
            else:
                compat = True
            if kwargs.has_key('varname'):
                varname = kwargs['varname']
            else:
                varname = 'var0'
            data = ncdf.r(source, compat=compat, varname=varname)
        elif type(source) == n.ndarray:
            self.bigname = 'array'
            # data is provided
            data = source

        if type(location) in [list, tuple, n.ndarray]:
            # co-ords supplied - FILL IN
            self.big_offset = tuple(location)
            if kwargs.has_key('boxsize'):
                boxsz = kwargs['boxsize']
            else:
                boxsz = (30, 30, 10)
            sp = extract_box( data, tuple( location ), tuple( boxsz ) )
            self.source_data = sp
            self.picked = True
        else:
            # click for co-ords - assume args[1] is pv - FILL IN
            location.AddImg(data, self.bigname)
            location.DisconnectClickCallback() # make sure no others registered
            location.ConnectClickCallback(self._pick_data_cb, data)

    def _pick_data_cb(self, x, y, z, *args):
        '''callback routine for pick_data - shouldn''t need user calling
        '''
        print "Got click at ", x, y ,z
        self.big_offset = int(x), int(y), z
        big = args[0]
        boxsz = array( (30, 30, 10) )
        coords = (x, y, z)
        sp = extract_box( big, tuple( coords ), tuple( boxsz ) )
        self.source_data = sp
        self.picked = True
        
    def generate_model(self):
        # check all params are entered
        needed = ( 'dend_diam', 'head_diam', 'neck_diam', 'neck_length', \
                   'angle_roundz', 'angle_toz', 'angle_den' )
        for i in needed:
            if not i in self.params.keys():
                if i == 'head_diam':
                    self.params['head_diam'] = 0.5
                elif i == 'neck_diam':
                    self.params['neck_diam'] = 0.25
                else:
                    print i, "is needed to generate a model."
                    return
            elif n.isnan(self.params[i]):
                if i == 'head_diam':
                    self.params['head_diam'] = 0.5
                elif i == 'neck_diam':
                    self.params['neck_diam'] = 0.25
                else:
                    print i, "is needed to generate a model."
                    return
        # construct model
        # need to get dvec, dvec_a again...
        dvec_a, dvec = spine.fit.fit_dendrite_vector(self.rotated_data, \
                                                     threshold=fit_thresh)
        print "Dendrite vector:", dvec, "\n & point:", dvec_a
        dvec   *= self.pxspc / model_pxres
        dvec_a *= self.pxspc / model_pxres

        # pre-calculate some constant stuff for use in the loops
        model_size = (n.asarray(self.rotated_data.shape) * \
                      self.pxspc / model_pxres).round().astype(int)

        # fit angles and origin
        t0 = 0.
        arz0 = self.params['angle_roundz']
        atz0 = self.params['angle_toz']

        # fit diameters and neck length
        mid_pxspc    = n.zeros(3) + model_pxres * 4
        mid_fwhms    = mid_pxspc * 2
        mid_fwhms_px = mid_fwhms / model_pxres
        blur1        = gauss3d(mid_fwhms_px[0], mid_fwhms_px[1], \
                               mid_fwhms_px[2], 1e-3)
        second_fwhms = n.sqrt(self.fwhm**2 - mid_fwhms**2) / mid_pxspc
        blur2        = gauss3d(second_fwhms[0], second_fwhms[1], \
                               second_fwhms[2], 1e-3)

        A0  = self.rotated_data.sum()
        nd0 = self.params['neck_diam']
        hd0 = self.params['head_diam']
        nl0 = self.params['neck_length']
        dd0 = self.params['dend_diam']
        o = dvec_a + t0 * dvec
        angle_den = n.arctan( dvec[2] / n.sqrt( dvec[0]**2 + dvec[1]**2 ) )
        model = spine.model.build( \
            px_res = spine.utils.model_pxres, \
            pars = {'origin_x'    : o[0], \
                    'origin_y'    : o[1], \
                    'origin_z'    : o[2], \
                    'angle_den'   : angle_den, \
                    'dend_diam'   : self.params['dend_diam'], \
                    'head_diam'   : self.params['head_diam'], \
                    'neck_diam'   : self.params['neck_diam'], \
                    'neck_length' : self.params['neck_length'], \
                    'angle_roundz': self.params['angle_roundz'], \
                    'angle_toz'   : self.params['angle_toz']
                    }, \
            vol_size = model_size)
        return model

    def display3d(self, what='data', isoval=0.08, show_dvec=False):
        if enth:
            sc = mlab.figure()
            sc.scene.background = (1.,1.,1.)
            sc.scene.foreground = (0.,0.,0.)
            src = ArraySource()
            if what == 'data' or what == 'rotated':
                src.spacing = list(self.pxspc)
            elif what == 'model':
                src.spacing = n.array(spine.utils.model_pxres).repeat(3)
            src._update_image_data_fired()
            if what == 'data':
                src.scalar_data = self.source_data
            elif what == 'rotated':
                if self.rotated_data <> None:
                    src.scalar_data = self.rotated_data
                else:
                    print "Rotated data does not exist. Run calc_params()."
            elif what == 'model':
                m = self.generate_model()
                src.scalar_data = m
            if show_dvec:
                # need to get dvec, dvec_a
                dvec_a, dvec = spine.fit.fit_dendrite_vector( \
                    src.scalar_data, threshold=fit_thresh)
                print "Dendrite vector:", dvec, "\n & point:", dvec_a
                dvec *= self.pxspc# / model_pxres
                dvec = unitvec(dvec)
                dvec_a *= self.pxspc# / model_pxres
                dvecsrc = ScatterSource()
                dvecsrc.points = dvec_a[n.newaxis,...]
                dvecsrc.vector_data = dvec[n.newaxis,...]
                gly = Glyph()
                dvecsrc.add_child(gly)
                sc.add_child(dvecsrc)
            iso = IsoSurface()
            out = Outline()
            src.add_child(iso)
            src.add_child(out)
            sc.add_child(src)
            iso.contour.contours=[isoval]
            return sc

    def test_origins(self, single=True):
        sc = self.display3d('rotated', 0.05)
        
        # need to get dvec, dvec_a again...
        #dvec_a, dvec = spine.fit.fit_dendrite_vector(self.rotated_data,
        #                                             threshold=fit_thresh)
        #dvec *= self.pxspc
        #dvec = unitvec(dvec)
        #dvec_a *= self.pxspc
        dvec_a, dvec = self.dvecr_a, self.dvecr
        print "Dendrite vector:", dvec, "\n & point:", dvec_a

        orig = ScatterSource()
        ogly = Glyph()
        orig.points = dvec_a[n.newaxis]
        orig.add_child(ogly)
        ogly.glyph.glyph.scale_factor = 0.2
        ogly.glyph.glyph_source = ogly.glyph.glyph_list[1]
        orig.vector_data = dvec[n.newaxis,...]
        ogly.glyph.color_mode = 0 #"no_coloring"
        ogly.glyph.scale_mode = 3 #"no_coloring"
        sc.add_child(orig)
        
        # pre-calculate some constant stuff for use in the loops
        model_size = (n.asarray(self.rotated_data.shape) * \
                      self.pxspc / model_pxres).round().astype(int)
        print "model_size", model_size

        # fit angles and origin
        #t0 = 0. # need a better guess at initial position for this
        # manual iteration through possible values to find minimum resid
        # - need to identify limits to t because 0 may be off to one side
        backlim = n.where((dvec > 0), n.zeros((3,)), \
                          model_size * model_pxres)
        print "backlim", backlim
        forwardlim = n.where((dvec < 0), n.zeros((3,)), \
                             model_size * model_pxres)
        print "forwardlim", forwardlim
        
        tstep = 0.25  # in multiples of dvec - i.e. microns
        tsteps = [0.] # start at origin, at least 
        nbacksteps = ((backlim - dvec_a) / \
                      (-dvec * tstep)).min().round().astype(int)
        backsteps = range(-1, -nbacksteps, -1)
        print "back steps", backsteps
        backsteps.reverse()
        tsteps = backsteps + tsteps
        nforwardsteps = ((forwardlim - dvec_a) / \
                         (dvec * tstep)).min().round().astype(int)
        print "back steps", nbacksteps, "forward steps", nforwardsteps
        forwardsteps = range(1, nforwardsteps)
        print "forward steps", forwardsteps
        tsteps += forwardsteps
        tsteps = n.asarray(tsteps) * tstep
        print "t steps", tsteps
        print "Finding optimal t..."
        points = (dvec_a[n.newaxis,...] + dvec[n.newaxis,...] * \
                 tsteps[...,n.newaxis])
        # calculates points in _micron_ through volume
        
        print "first and last points =", points[[0,-1]]
        opts = ScatterSource()
        opts.points = points
        opts.scalar_data = n.ones(points.shape[0])
        
        gly = Glyph()
        opts.add_child(gly)
        gly.glyph.glyph_source = gly.glyph.glyph_list[4] # works
        gly.glyph.glyph.scale_factor = 0.11              # works
        gly.glyph.color_mode = 0 # "no_coloring"
        gly.glyph.scale_mode = 3 # "no_scaling"
        sc.add_child(opts)
        ogly.actor.actor.property.color = (0.0, 1.0, 0.0)
        gly.actor.actor.property.color = (1.0, 0.0, 0.0) # works
