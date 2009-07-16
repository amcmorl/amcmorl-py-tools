import numpy
import coordhandling
import scipy.ndimage
import pyvis
import ncdf
from glob import glob
from rebin import congrid
import numpy as n
posmods = ['pupynere', 'pynetcdf',
           'Scientific.IO.NetCDF', 'Scientific_NetCDF']
tried = 0
for mod in posmods:
    try:
        exec("from %s import NetCDFFile" % mod)
        break
    except ImportError:
        tried += 1
if tried == len(posmods):
    raise ImportError("No module containing NetCDFFile found.")
import os.path
import wx
import numpil
import matplotlib
from copy import copy
matplotlib.interactive(False)
matplotlib.use("WXAgg")
matplotlib.rcParams['text.usetex'] = False
import Image

from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg, \
     FigureManager, NavigationToolbar2WxAgg
from matplotlib.figure import Figure
from matplotlib.axes import Subplot

def bytescl(name):
    maxz = 96
    ncfile = NetCDFFile(name, 'r')
    var = ncfile.variables['reconstruction']
    testxy = n.asarray(var[0,:,:])
    testz = n.asarray(var[:,0,0])
    sh = (testxy.shape[1], testxy.shape[0], testz.shape[0])
    print "stack shape", sh
    arsize = sh[0] * sh[1] * sh[2]
    if arsize > (1024 * 1024 * maxz):
        tedges = n.arange(int(sh[2] / maxz) + 1) * maxz
        bedges = (n.arange(int(sh[2] / maxz) + 1) + 1) * maxz - 1
        if bedges[-1] > sh[2] - 1:
            bedges[-1] = sh[2] - 1
        #hz = sh[2] / 2
        #print "half way point", hz
        stackmax = 0
        for i in xrange(len(list(tedges))):
            print "Loading", name, "from", tedges[i], "to", bedges[i], "..."
            stack = n.asarray(var[tedges[i]:bedges[i],:,:])
            thismax = stack.max()
            print "max of sub-volume is", thismax
            if thismax > stackmax:
                stackmax = thismax
        print "max of volume is", stackmax
        conv = 255. / stackmax
        
        print "writing ncdf"
        nc_out = NetCDFFile(name.split('.')[0] + '_byte.nc', 'w')
        print "creating dimensions"
        dimslist = ['dim0', 'dim1', 'dim2']
        for i, dim in enumerate(dimslist):
            nc_out.createDimension( dim, sh[i] )
        print "determining typecode",
        #        if stack.dtype.type == numpy.uint8:
        typecode = 'i'
        print "is", typecode
        #else:
        #    typecode = 'd'
        #    print " is d"
        print "creating variable"
        varout = nc_out.createVariable( 'var0', typecode, tuple(dimslist) )

        for i in xrange(len(list(tedges))):
            print "Loading", name, "from", tedges[i], "to", bedges[i], "..."
            stack = n.array(var[tedges[i]:bedges[i],:,:])
            print "Converting sub-volume..."
            stack *= conv
            print "assigning data to volume of z size", \
                  tedges[i], "to", bedges[i]
            varout[:,:,tedges[i]:bedges[i]] = stack.astype( \
                numpy.int32).transpose(2,1,0)
        print "closing file"
        nc_out.close()

    else:
        stack = n.array(var[:])
        conv = 255. / stack.astype(float).max()
        stack *= conv
        ncdf.w(name.split('.')[0] + '_byte.nc', \
               stack.astype(numpy.uint8).transpose((2,1,0)))

def zeromin(x):
    '''makes sure x is lower-bounded by zero'''
    if x < 0:
        return 0
    else:
        return x

vzeromin = numpy.vectorize(zeromin)

def big_max_projection(name):
    '''calculates a maximum projection of a very large volume
    by loading the volume in two separate runs'''
    maxz = 96
    ncfile = NetCDFFile(name, 'r' )
    var = ncfile.variables['var0']

    testxy = var[:,:,0]
    testz = var[0,0,:]
    sh = (testxy.shape[0], testxy.shape[1], testz.shape[0])
    maxes = n.zeros((sh[0],sh[1],2),dtype = n.uint8)
    print "stack shape", sh   
    tedges = n.arange(int(sh[2] / maxz) + 1) * maxz
    bedges = (n.arange(int(sh[2] / maxz) + 1) + 1) * maxz - 1
    if bedges[-1] > sh[2] - 1:
        bedges[-1] = sh[2] - 1
    for i in xrange(len(list(tedges))):
        print "Loading", name, "from", tedges[i], "to", bedges[i], "..."
        stack = var[...,tedges[i]:bedges[i]]
        maxes[...,1] = stack.max(2)
        maxes[...,0] = maxes.max(2)
    return maxes[...,0], testz.shape[0]

    

def join(a, pta, b, ptb=None, method='avg'):
    '''creates a new projection in which projections of vola and volb
    are overlapped such that pta in vola is mapped on to ptb in volb'''

    if ptb == None:
        ptb = (0,0)
    if not type(pta) == numpy.ndarray:
        pta = numpy.array(pta)
    if not type(ptb) == numpy.ndarray:
        ptb = numpy.array(ptb)
    
    if not method in ['avg', 'a', 'b']:
        raise ValueError('Method must be one of \'avg\', \'a\' or \'b\'')

    asz = numpy.array(a.shape)
    bsz = numpy.array(b.shape)

    maxofs = numpy.vstack((pta, ptb)).max(0)
    bigsz = maxofs + numpy.vstack((asz - pta, bsz - ptb)).max(0)

    # global co-ordinates
    lca = maxofs - pta
    uca = lca + asz
    lcb = maxofs - ptb
    ucb = lcb + bsz

    big = numpy.zeros( tuple(bigsz) )

    if method == 'avg':
        print "Averaging common region..."
        mask = numpy.zeros( tuple(bigsz), dtype=bool)
        mask[lca[0]:uca[0], \
             lca[1]:uca[1]] = (a > 0)
        mask[lcb[0]:ucb[0], \
             lcb[1]:ucb[1]] = mask[lcb[0]:ucb[0], \
                                   lcb[1]:ucb[1]] & (b > 0)
        big[lca[0]:uca[0], \
            lca[1]:uca[1]] = a
        big[lcb[0]:ucb[0], \
            lcb[1]:ucb[1]] += b
        big[numpy.where(mask == True)] /= 2
        
    elif method == 'a':
        print "Placing first dataset in common region..."
        big[lcb[0]:ucb[0], \
            lcb[1]:ucb[1]] = b
        big[lca[0]:uca[0], \
            lca[1]:uca[1]] = a

    elif method == 'b':
        print "Placing second dataset in common region..."
        big[lca[0]:uca[0], \
            lca[1]:uca[1]] = a
        big[lcb[0]:ucb[0], \
            lcb[1]:ucb[1]] = b
    return big

def montage(vola, vac, volb, vbc, subsz=30, disp=False, method='avg'):
    '''takes 2 stacks, with points at a common feature, calculates
    the correlation offset between the two and returns a montaged
    image containing information from both stacks

    Parameters:
      vola = first volume
      vac = point in vola
      volb = second volume
      vbc = point in volb
      subsz = size of boxes to correlate (=30)
      disp = whether to display result (=False)
    '''
    vas = coordhandling.extract_box( vola, vac, subsz )
    vbs = coordhandling.extract_box( volb, vbc, subsz )
    
    cor = scipy.ndimage.correlate( vas, vbs )
    coroff = numpy.where(cor == cor.max())

    voffs = vac + coroff - vbc
    print "Offset = ", voffs

    if disp:
        print "Creating montage volume..."
        nv = join(vola, voffs, volb, method=method)
        return nv, voffs
    else:
        return voffs

def extract_vol(name, vol_pt, subsz, zsz):
    vol_start = numpy.hstack((vol_pt - subsz/2, 0))
    vol_start = vzeromin(vol_start)
    vol_fullz_sz = numpy.hstack((subsz, subsz, zsz))
    vol_bit_fullz = ncdf.r(name, start=vol_start, size=vol_fullz_sz)
    print "loading volume of shape", vol_bit_fullz.shape
    vol_ptz = coordhandling.zmax(vol_bit_fullz, (subsz/2, subsz/2))
    vol_range = numpy.array((vol_ptz - subsz / 2, vol_ptz + subsz / 2))
    vol_range = vol_range.clip(min=0, max=zsz)
    vol_bit = vol_bit_fullz[:,:,vol_range[0]:vol_range[1]]
    print "correlating volume of shape", vol_bit.shape
    return vol_bit, vol_ptz


class PlotFigure(wx.Frame):
    def __init__(self, parent):
        self.parent = parent
        wx.Frame.__init__(self, None, -1, "Neuron stack tiler")

        self.fig = Figure((9,8), 75)
        self.canvas = FigureCanvasWxAgg(self, -1, self.fig)
        self.toolbar = NavigationToolbar2WxAgg(self.canvas)
        self.toolbar.Realize()

        self.figmgr = FigureManager(self.canvas, 1, self)
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.canvas, 1, wx.LEFT|wx.TOP|wx.GROW)
        sizer.Add(self.toolbar, 0, wx.GROW)
        self.SetSizer(sizer)
        self.Fit()
        self.canvas.mpl_connect('button_press_event', self.ClickHandler)

    def plot_data(self, img1, name1, img2, name2):
        self.aname = name1
        self.bname = name2

        a = self.fig.add_subplot(121)
        a.imshow(numpy.flipud(img1.transpose()))
        a.set_title(name1)
        b = self.fig.add_subplot(122)
        b.imshow(numpy.flipud(img2.transpose()))
        b.set_title(name2)
        self.toolbar.update()

    def GetToolBar(self):
        return self.toolbar

    def ClickHandler(self, evt):
        #if not (self.mainframe.toolbar.panmode or \
        #        self.mainframe.toolbar.zoomode):
        if self.toolbar.mode == '':
            if evt.inaxes:
                x = evt.xdata
                y = evt.ydata
                if evt.inaxes.title.get_text() == self.aname:
                    self.parent.vola_pt = x,y
                    print "Click was in vol a at pt", x, y
                else:
                    self.parent.volb_pt = x,y
                    print "Click was in vol b at pt", x, y
            else:
                print "Click was outside axes"            
  
class Tiler:
    def __init__(self, aname, bname, method='avg', subsz=30):
        '''gets in two stacks, takes user input to select a common
        feature in both images and returns the optimal offset between
        them i.e. the co-ordinates of the (0,0,0) point of volb in
        the vola image space'''
        if os.path.isfile(aname):
            self.aname = aname
        else:
            raise ValueError(aname + " is not a valid file.")
        if os.path.isfile(bname):
            self.bname = bname
        else:
            raise ValueError(bname + " is not a valid file.")
        self.subsz = subsz
        self.method = method
        print "opening", aname, "..."
        self.vola_proj, self.vola_zsz = big_max_projection(aname)

        print "opening", bname, "..."
        self.volb_proj, self.volb_zsz = big_max_projection(bname)

        self.frame = PlotFigure(self)
        self.frame.Show(True)
        self.frame.plot_data(self.vola_proj, aname, self.volb_proj, bname)
        print "Select points in each volume"
        print "then run \"run_montage()\""
       
    def run_montage(self):
        if (not self.vola_pt is None) and (not self.volb_pt is None):
            subsz = self.subsz
            
            vola_pt = numpy.floor(self.vola_pt).astype(int)
            vola_bit, vola_ptz = extract_vol(self.aname, vola_pt, \
                                  subsz, self.vola_zsz)
            volb_pt = numpy.floor(self.volb_pt).astype(int)
            volb_bit, volb_ptz = extract_vol(self.bname, volb_pt, \
                                   subsz, self.volb_zsz)

            print "correlating subvolumes..."
            cor = scipy.ndimage.correlate( vola_bit.astype(float), \
                                           volb_bit.astype(float) )
            corov = numpy.hstack(numpy.where(cor == cor.max()))
            corcorov = corov - subsz/2

            apt = numpy.hstack((vola_pt, vola_ptz))
            bpt = numpy.hstack((volb_pt, volb_ptz))

            print "corrected correlation offset", corcorov

            self.ptb_in_a = apt + corcorov - bpt

            return self.ptb_in_a
            #result = montage(self.vola, pta, volb, ptb, disp=False, \
            #                 method=method)
        
    def get_result(self,method='avg'):
        return join(self.vola_proj, self.ptb_in_a[0:2], self.volb_proj, \
                    method=method)

def calculate_global_offsets(fname):
    '''calculates global offsets for each stack from relative offsets
    to neighbouring stacks

    Parameters:
      fname : filename of offsets file
    '''
    
    f = open(fname)
    lines = f.readlines()
    unglobd = []
    globd = {}
    for line in lines:
        if not line[0] == '#' or line.strip() == '':
            bits = line.split()
            entry = { 'rel':bits[1],
                     'ref':bits[0] }
            unglobd.append(entry)
            if bits[2].isalpha():
                entry['rel'] += '|' + bits[2]
                entry['slice'] = bits[3]
                entry['offset'] = n.asarray([int(x) for x in bits[4:7]])
                # need to handle split stacks
            else:
                # next three entries are x y z
                entry['offset'] = n.asarray([int(x) for x in bits[2:5]])
    while len(unglobd) > 0:
        entry = unglobd.pop(0) # take first item in list
        if globd.has_key(entry['ref']):
            globd[entry['rel']] = globd[entry['ref']] + entry['offset']
        elif entry['ref'] == 'glob':
            globd[entry['rel']] = n.asarray((0,0,0))
        else:
            unglobd.append(entry)
    return globd 

def get_stacks(gdir, daydirs=True, ending='_byte.nc'):
    '''returns paths to all stacks (assumed to have _byte.nc at end
    of filename in dated (6 digit named) directories'''
    if daydirs:
        dirs = glob(gdir+'[0-9][0-9][0-9][0-9][0-9][0-9]/')
    else:
        dirs = [gdir]
    stacks = []
    for wdir in dirs:
        stacks += glob(wdir + '*' + ending)
    return stacks

def generate_montage(gdir, scale=3):
    '''constructs a maximum projection montage from stack data and
    offsets file'''
    offsets_file = gdir + 'stack_offsets.txt'
    global_offsets = calculate_global_offsets(offsets_file)
    stacks = get_stacks(gdir)
    
    # match stack name with offset
    for short_code in global_offsets.keys():
        short_name = short_code.split('|')[0]
        i = 0
        found = False
        while not found:
            #if stacks[i].find(short_name) <> -1:
            if os.path.basename(stacks[i]).split('_')[0] == short_name:
                found = True
                global_offsets[short_code] = \
                                           [global_offsets[short_code], \
                                            stacks[i]]
            i += 1

    # calculate biggest size of array
    # assume 1024x1024 + biggest offset

    # find biggest x and biggest y offsets
    maxx, maxy = 0, 0
    minx = miny = 0
    for info in global_offsets.itervalues():
        maxx = n.maximum(maxx, info[0][0])
        maxy = n.maximum(maxy, info[0][1])
        minx = n.minimum(minx, info[0][0])
        miny = n.minimum(miny, info[0][1])

    bigsz = n.asarray((maxx + 1024 + n.abs(minx), \
                       n.abs(miny) + 1024 + maxy)) / float( scale )
    
    big = n.zeros(n.hstack((bigsz,2)))
    for info in global_offsets.itervalues():
    #sts = ['upa', 'som', 'rta', 'upcd', 'rtb', 'rtc', 'dnr', 'dnba', \
    #       'dnc', 'dnl', 'dnga', 'dngb']
    #for info in [global_offsets[x] for x in sts]:
        dat = ncdf.r(info[1]).max(2)
        datsz = n.asarray(dat.shape)
        datnewshape = n.round(datsz / float( scale )).astype(int)
        dat = congrid(dat, datnewshape)
        startx = (-minx + info[0][0]) / float( scale )
        starty = (-miny + info[0][1]) / float( scale )
        
        big[startx:startx + datnewshape[0], \
            starty:starty + datnewshape[1],1] = dat
        big[:,:,0] = big.max(2)
    return big.max(2)

def make_max_projections(wdir):
    files = glob(wdir + '/*/*_byte.nc')
    for f in files:
        ar = ncdf.r(f)
        armax = ar.max(2)
        newname = f.split('.')[0] + '_max.png'
        armax *= 255 / armax.max()
        im = numpil.numpy2PIL(armax.astype(n.uint8))
        print "Writing %s..." % newname
        im.save(newname)

def make_max_projection(f):
    ar = ncdf.r(f)
    armax = ar.max(2)
    newname = f.split('.')[0] + '_max.png'
    armax *= 255 / armax.max()
    im = numpil.numpy2PIL(armax.astype(n.uint8))
    print "Writing %s..." % newname
    im.save(newname)

def get_max_projections(gdir):
    '''returns paths to all max projections - files two directories deep,
    with name ending in _max.png'''
    stacks = glob(gdir+'*/*_max.png')
    return stacks

def generate_max_projections_montage(gdir, scale=3):
    '''constructs a maximum projection montage from maximum projections
    and offsets file'''
    offsets_file = gdir + 'stack_offsets.txt'
    global_offsets = calculate_global_offsets(offsets_file)
    stacks = get_max_projections(gdir)
    
    # match stack name with offset
    for short_code in global_offsets.keys():
        short_name = short_code.split('|')[0]
        i = 0
        found = False
        while not found:
            if stacks[i].find(short_name) <> -1:
                found = True
                global_offsets[short_code] = \
                                           [global_offsets[short_code], \
                                            stacks[i]]
            i += 1

    # calculate biggest size of array
    # assume 1024x1024 + biggest offset

    # find biggest x and biggest y offsets
    maxx, maxy = 0, 0
    minx = miny = 0
    for info in global_offsets.itervalues():
        maxx = n.maximum(maxx, info[0][0])
        maxy = n.maximum(maxy, info[0][1])
        minx = n.minimum(minx, info[0][0])
        miny = n.minimum(miny, info[0][1])

    bigsz = n.asarray((maxx + 1024 + n.abs(minx), \
                       n.abs(miny) + 1024 + maxy)) / float( scale )
    
    big = n.zeros(n.hstack((bigsz,2)))
    for info in global_offsets.itervalues():
        #sts = ['upa', 'som', 'rta', 'upcd', 'rtb', 'rtc', 'dnr', 'dnba', \
        #       'dnc', 'dnl', 'dnga', 'dngb']
        #for info in [global_offsets[x] for x in sts]:
        dat = numpil.PIL2numpy(Image.open(info[1]))
        datsz = n.asarray(dat.shape)
        datnewshape = n.round(datsz / float( scale )).astype(int)
        dat = congrid(dat, datnewshape)
        startx = (-minx + info[0][0]) / float( scale )
        starty = (-miny + info[0][1]) / float( scale )
        
        big[startx:startx + datnewshape[0], \
            starty:starty + datnewshape[1],1] = dat
        big[:,:,0] = big.max(2)
    return big.max(2)
