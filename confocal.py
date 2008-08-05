import numpy as n
import string, os, glob, Image, numpil
# required for imglob and load_stack
import pylab, gaussian, scipy.stats
# additionally required for meas_fwhms
import processing
# required for align_stack
import unittest
import pyvis

memlimit = 100663296 # maximum number of elements to load at once

def imglob(exf, ndig=3):
    ''' returns list of files matching pattern of exf
    with remaining digits of filename stripped and replaced with
    ndig * [0-9]'''
    
    f_bits = os.path.splitext(exf)
    f_noext = f_bits[0] # extract base file without extension
    #f_nodig = f_noext.rstrip(string.digits) # remove all digits at the end
    f_nodig = f_noext[:-ndig] # remove last three characters
    f_glob = f_nodig + ''.join(['[0-9]' for x in range(ndig)]) + \
             f_bits[1]
    #print "Glob string", f_glob
    fs = glob.glob(f_glob)
    fs.sort()
    return fs


def load_stack(exf, ndig=3, verbose=False, start=None, stop=None):
    ''' loads a stack of confocal images including the file exf
    and with ndig * number of digits in filename

    start = first image of stack
    stop  = last image of stack + 1'''
    
    # find stack files
    files = imglob(exf, ndig)[start:stop]
    z = len(files)
    
    # load first image to get shape and type
    if len(files) > 0:
        if verbose:
            print "Opening file", files[0]
        im = Image.open(files[0])
        x,y = im.size

        # load stack files
        sl = numpil.PIL2numpy( im )
        ar = n.zeros((y,x,z), dtype=sl.dtype)
        ar[...,0] = sl
        files.pop(0)
        for i, fn in enumerate( files ):
            if verbose:
                print "Opening file", fn
            im = Image.open( fn )
            ar[...,i + 1] = numpil.PIL2numpy( im )
        return ar
    else:
        print "No files found"
        return None

def stack_size(exf, ndig=3, verbose=False, start=None, stop=None):
    # find stack files
    files = imglob(exf, ndig)[start:stop]
    z = len(files)
    
    # load first image to get shape and type
    if len(files) > 0:
        if verbose:
            print "Opening file", files[0]
        im = Image.open(files[0])
        x,y = im.size
    return n.array((x,y,z))

def load_image(file):
    im = n.array(Image.open(file))
    return im


def write_stack(stack, namebase):
    '''writes a series of tif images, one for each plane from stack,
    with a name in the format namebase%%%.tif'''

    znum = stack.shape[2]
    for i in xrange(znum):
        fname = "%s%03d.tif" % (namebase, i)
        print "Writing %s..." % fname
        im = numpil.numpy2PIL(stack[:,:,i])
        im.save(fname)


def make_overview(firstfile, ndig=3, verbose=False):
    files = imglob(firstfile, ndig)
    z = len(files)
    maxz = 96
    
    # load first image to get shape and type
    if len(files) > 0:
        if verbose:
            print "Opening file", files[0]
        im = Image.open(files[0])
        y,x = im.size
        if x*y*z < memlimit:
            stack = load_stack(firstfile, ndig)
            over = stack.max(2)
        else:
            maxes = n.zeros((x,y,2), dtype=n.uint8)
            tedges = n.arange(int(z / maxz) + 1) * maxz
            bedges = (n.arange(int(z / maxz) + 1) + 1) * maxz - 1
            if bedges[-1] > z - 1:
                bedges[-1] = z - 1
            for i in xrange(len(list(tedges))):
                stack = load_stack(firstfile, start=tedges[i], stop=bedges[i])
                maxes[...,1] = stack.max(2)
                maxes[...,0] = maxes.max(2)
            over = maxes[...,0]
        return over
    else:
        return None

def pick_cent(coords):
    '''Picks middle (index-wise) of the co-ordinates returned by numpy.where'''
    xar, yar, zar = coords
    nc = xar.shape[0]
    hnc = nc / 2
    return n.array([xar[hnc]]), \
           n.array([yar[hnc]]), \
           n.array([zar[hnc]])
    

def meas_fwhms(psf, xyspc, zspc, offset, graphic=False):
    '''Measures the FWHMs of a psf, using the brightest point as the centre.
    Seems to underestimate the peak slightly, which will overestimate the
    FWHM. The PSF should be a 3-D volume (e.g. an averaged psf), probably
    loaded with confocal.load_stack(''first_img_name.tif'')  '''

    cent = pick_cent( n.where( psf == psf.max() ) )
    xprof = psf[:, cent[1], cent[2]].squeeze()
    xpts = n.arange(xprof.shape[0]) * xyspc
    zprof = psf[cent[0], cent[1], :].squeeze()
    zpts = n.arange(zprof.shape[0]) * zspc

    if graphic:
        pylab.clf()
        pylab.plot(xpts, xprof, label='x data')
        pylab.plot(zpts, zprof, label='z data')

    Ax, fwhmx, cx = gaussian.fitgauss1d( xprof - offset, r=xpts, \
                                         p0=(xprof.max(), \
                                             0.5 * xpts.max(), \
                                             int(cent[0]) * xyspc ) )
    Az, fwhmz, cz = gaussian.fitgauss1d( zprof - offset, r=zpts, \
                                         p0=(zprof.max(), \
                                             0.5 * zpts.max(), \
                                             int(cent[2]) * zspc) )
    if graphic:
        print "   %5s  %5s %5s" % ("Peak", "FWHM", "Centre")
        print "x: %5.2f %5.2f %5.2f" % (Ax, fwhmx, cx)
        print "z: %5.2f %5.2f %5.2f" % (Az, fwhmz, cz)
    
        xmodel = gaussian.gauss1d(xpts, Ax, fwhmx, cx) + offset
        zmodel = gaussian.gauss1d(zpts, Az, fwhmz, cz) + offset
        
        pylab.plot(xpts, xmodel, label='x fit')
        pylab.plot(zpts, zmodel, label='z fit')
        
        pylab.legend()

    return (fwhmx, fwhmz)


# ----------------- Helper routines for stack watcher ------------------------

def pixel_stats(pixels):
    '''calculates signal-to-noise ratio for pixels (could be an image)
    or subscripted region of an image'''
    mean = pixels.mean()
    vari = scipy.stats.var( pixels.flatten() )
    print "n = %d, mean = %0.1f var = %0.2f" % \
          (pixels.size, mean, vari )
    print " => photon number = %0.1f" % (mean ** 2 / vari)
    print " => SNR = %f" % (mean / n.sqrt(vari))


def shift_image(image, offset, pad_val):
    '''Shifts an image by offset, padding empty spaces with pad_val.

    Description:
    
       Returns a shifted version of an image, viewed through the same
       ''window'' as the original. Empty spaces created by the shift
       are filled with pad_val

    Inputs:

       image   -- a 2-D array
       offset  -- a (2,) array
       pad_val -- (default == 0) the value to pad the slices by when shifting

    Outputs:

       out -- a 2-D array of the same shape and dtype as image, with
              values of image shifted by offset amount.
    '''
    res = n.zeros_like(image) + pad_val

    # handle x
    if offset[0] < 0:
        destxlo, destxhi = None, offset[0]
        src_xlo, src_xhi = -offset[0], None
    elif offset[0] == 0:
        destxlo, destxhi = None, None
        src_xlo, src_xhi = None, None
    else: # offset[0] > 0
        destxlo, destxhi = offset[0], None
        src_xlo, src_xhi = None,-offset[0]

    # handle y
    if offset[1] < 0:
        destylo, destyhi = None, offset[1]
        src_ylo, src_yhi = -offset[1], None
    elif offset[1] == 0:
        destylo, destyhi = None, None
        src_ylo, src_yhi = None, None
    else: # offset[1] > 0
        destylo, destyhi = offset[1], None
        src_ylo, src_yhi = None,-offset[1]

#     print "res[%s:%s, %s:%s] = image[%s:%s, %s:%s]" % \
#           (destxlo, destxhi, destylo, destyhi, \
#            src_xlo, src_xhi, src_ylo, src_yhi)
    res[destxlo:destxhi, destylo:destyhi] = image[src_xlo:src_xhi, \
                                                  src_ylo:src_yhi]
    return res


def align_stack(stack, pad_val = 0):
    '''Aligns adjacent images of a stack to correct for slight shift.

    Description:
    
       Calculates the correlation between adjacent images of a stack
       and shifts the lower image to minimize the correlated offset
       along the z (3rd) axis.

    Inputs:

       stack   -- a 3-dimensional array
       pad_val -- (default == 0) the value to pad the slices by when shifting

    Outputs:

       out -- a 3-dimensional array of the same shape and dtype as stack
              with aligned z-levels.
    '''
    fig = pylab.figure(1)
    ax = fig.add_subplot(111)
    res = n.zeros_like(stack)
    offsets = n.zeros((2,stack.shape[2]))
    centre = n.array(stack[:,:,0].shape) / 2.
    res[:,:,0] = stack[:,:,0]
    for i in xrange(stack.shape[2] - 1):
        cor = processing.fft_correlate2(stack[:,:,i], stack[:,:,i+1])
        ax.imshow(cor)
        pylab.show()
        offsets[:,i] = n.array(n.where(cor == cor.max())).squeeze() \
                       - centre
        print "z-level: %2d, offset %2d, %2d; shift: %2d, %2d" % \
              ( i, offsets[0,i], offsets[1,i], \
                offsets.sum(1)[0], offsets.sum(1)[1] )
        
        res[:,:,i+1] = shift_image(stack[:,:,i+1], offsets.sum(1), \
                                   pad_val = pad_val)
        # sum over the 1st dim of offsets to account for previous shifts

    return res, offsets


class stack_watcher:
    '''observes an incoming stack and reports
    various image statistics to the user upon request

    Inputs:
      name - string name of first file in stack
      bg   - dict, lower and upper corners of background square
             specified as two tuples'''

    def __init__(self, name, bg={'start':None, 'stop':None}):
        self.count = 0
        self.name = name
        self.bg = bg
        self.bgs = []
        self.maxproj = None

    def update(self):
        all_files = imglob(self.name)
        new_files = all_files[self.count:]
        for f in new_files:
            print "Opening image %s" % (f)
            im = Image.open( f )
            self.latest = numpil.PIL2numpy( im )
            
            if self.maxproj is None:
                self.maxproj = self.latest
            else:
                self.maxproj = n.array((self.maxproj, self.latest)).max(0)

            if (self.bg['start'] is not None) and \
                   (self.bg['stop'] is not None):
                sl = [ slice(i,j) for (i,j) in \
                      zip(list(self.bg['start']), list(self.bg['stop'])) ]
                self.bgs.append(self.latest[sl].mean())
        self.count = len(all_files)
                
    def show_latest(self):
        pylab.imshow(self.latest)

    def show_max(self):
        pylab.imshow(self.maxproj)

    def plot_bg(self):
        pylab.plot(self.bgs)


class TestConfocalFunctions(unittest.TestCase):

    def test_pick_cent(self):
        a = n.zeros((10,10,10))
        a[4:7,5,5] = 1
        coords = n.where(a == 1)
        res = n.asarray(pick_cent(coords)).squeeze()
        self.assertTrue( n.all( res == n.asarray((5,5,5)) ) )

    def test_meas_fwhms(self):
        x, z = (3, 7)
        pxscl = 0.1
        tpsf = gaussian.gauss3d( x, x, z )
        measd = n.asarray( meas_fwhms(tpsf, pxscl, pxscl, 0) )
        self.assertTrue( n.allclose(measd, n.asarray((x, z)) * pxscl ))


def test():
    suite = unittest.TestLoader().loadTestsFromTestCase( \
        TestConfocalFunctions)
    unittest.TextTestRunner(verbosity=2).run(suite)
