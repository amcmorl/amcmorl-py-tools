import numpy
#from numpy.fft import fftshift
from scipy.signal import fftn, ifftn
from scipy.signal import conjugate, correlate
from scipy.fftpack import fftshift, ifftshift
from pylab import figure, close, show
import gaussian

def _centrek(k, sh):
    res = numpy.zeros(sh,dtype=float)
    ndims = len(sh)
    
    pta = numpy.round((numpy.array(sh) - numpy.array((k.shape))) / 2.)
    ptb = (pta + numpy.array((k.shape))).astype(int)
    execstr = 'res[' + \
              ','.join(['pta[%d]:ptb[%d]' % (x,x) for x in xrange(ndims)]) + \
              '] = k'
    exec( execstr )
    return res


def fftconvolve2(in1, in2=None, IN2=None):
    """Convolve two N-dimensional arrays using FFT. Saves on memory
    compared to scipy.signal.fftconvolve by not giving 'full' option
    and only calculating for the largest size in each dimension.

    Alternatively takes precomputed FFT of fftshifted input 2
    """
    s1 = numpy.array(in1.shape)
    s2 = numpy.array(in2.shape)
    size = tuple(numpy.maximum(s1, s2))
    #print "max size is", size
                 
    if in1.shape != size:
        #print "Centering image 1"
        in1 = _centrek(in1, size)
    if in2.shape != size:
        #print "Centering image 2"
        in2 = _centrek(in2, size)

    if (s1.dtype.char in ['D','F']) or (s2.dtype.char in ['D', 'F']):
        cmplx=True
    else: cmplx=False

    IN1 = fftn( in1, size )
    # IN2 = fftn( fftshift(in2), size )
    IN1 *= fftn( fftshift(in2), size )

    # IN1 *= IN2
    # del IN2
    ret = ifftn(IN1)
    del IN1
    if not cmplx:
        ret = numpy.real(ret)
    return ret


def test_fftconvolve2():
    pylab.clf()
    length = 200
    
    func = numpy.zeros(length)
    func[10:80] = 1
    func[120:150] = 2

    kern_len = 30
    kern = gaussian.gauss1d(numpy.arange(kern_len), 1, 11, kern_len/2.)
    kern /= kern.sum()

    theirs_direct = numpy.convolve(func, kern, mode='same')
    theirs_fft = scipy.signal.fftconvolve(func, kern, mode='same')
    mine = fftconvolve2(func, kern)

    diff = abs(theirs_fft - mine)
    slope = numpy.hstack( (abs(numpy.diff(theirs_fft)), (1e-6)) )
    ddiff = abs(slope - diff)

    # plot function
    #pylab.plot(func, 'g-', label='raw')

    # plot convolutions
    #pylab.plot(theirs_direct, 'r-', label='direct')
    #pylab.plot(theirs_fft, 'ro', label='fft', markersize=3)
    #pylab.plot(mine, 'bo-', label='mine', markersize=3)

    # plot difference
    #pylab.plot(diff, 'c-', label='diff')
    #pylab.plot(slope, 'co', label='slope', markersize=3)
    #pylab.plot(ddiff, 'm-', label='normalised diff')

    #pylab.legend()
    
    # plot kernels
    #a = pylab.axes([0.2,0.6,0.2,0.2])
    #pylab.title('kernels')
    #pylab.plot(kern, 'k-', label='kcen')
    # pylab.xlim((0,kern.shape[0]-1))

    if all( ddiff < 0.02 ):
        return True
    else:
        return False


def fft_correlate2(in1, in2):
    '''Calculates, using FFTs, the n-D correlation between in1 and in2.
    Saves on memory compared to scipy.signal.fftconvolve by not giving
    ''full'' option and only calculating for the largest size in each
    dimension. Only works for shape(in1) == shape(in2).'''

    this_shape = numpy.array(in2.shape)
    if in1.shape != in2.shape:
        raise ValueError("Input arrays must be of same shape.")
                 
    if (in1.dtype.char in ['D','F']) or (in2.dtype.char in ['D', 'F']):
        cmplx=1
    else: cmplx=0

    IN1 = fftn( in1, shape = this_shape )
    IN2 = conjugate(fftn( fftshift(in2), shape = this_shape ))

    IN1 *= IN2
    del IN2
    ret = ifftn(IN1)
    del IN1
    if not cmplx:
        ret = numpy.real(ret)
    return ret


def test_fft_correlate2():
    lg = 100.
#     fig = figure(1)
#     fig.clf()
#     ax1 = fig.add_subplot(211)
#     ax2 = fig.add_subplot(212)
    x = numpy.arange(lg)
    y1 = gaussian.gauss1d(x,1.,10.,lg/2.)
    y2 = gaussian.gauss1d(x,1.,10.,lg/2.+8.)
#     ax1.plot(x,y1,'b-')
#     ax1.plot(x,y2,'r-')
    
    theirs = correlate(y1, y2, mode='same')
    mine = fft_correlate2(y1, y2)

#     ax2.plot(x,theirs/theirs.max(),'b--o')
#     ax2.plot(x,mine/mine.max(),'r--x')

    if numpy.allclose( theirs, mine ):
        return True
    else:
        return False
