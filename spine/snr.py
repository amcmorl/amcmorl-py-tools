import numpy as _n, pylab as _p, Image, matplotlib as mpl
from scipy.stats import poisson, var, sem
from primitives import draw_circle
from gaussian import gauss2d
from scipy.ndimage import convolve, median_filter
from pyvis import pyvis

# -- helper functions --

gpsfsize=0.5

def add_noise(img, maxphotons, normalise=True, pc=0.95):
    sh = img.shape
    photons = img.flatten() / img.max() * maxphotons
    noisy = poisson.rvs(photons, \
                        size=photons.size).reshape(sh).astype(float)
    if normalise:
        #blowout = noisy.max() / photons.max()
        #print "blowout", blowout
        uplim = img.size * pc
        robustmax = noisy.flatten()[noisy.argsort(None)][uplim:].mean()
        #print "Real max =", noisy.max(), "robust max =", robustmax
        return noisy / robustmax
    else:
        return noisy

def get_slope(X, Y, n):
    x = X - X.mean()
    y = Y - Y.mean()
    return ((x * y).sum() + n * X.mean() * Y.mean()) \
           / ((x**2).sum() + n * X.mean()**2)

def build_data(imagesize = (512,512), beadsize=0.5, eps=1e-12, \
               pixelsize=0.05, psfsize=0.3, npmax=100, bg=5):
    radius = beadsize / pixelsize / 2.
    bdlocs = [(406, 163), (175, 270), (202, 216), (356, 475), (100, 434), \
              (236, 359), (433, 444), (295, 272), (93, 394), (61, 137)]
    img = _n.zeros(imagesize, dtype=bool)
    for loc in bdlocs:
        img = img | draw_circle(radius, size=imagesize, centre=loc)
    psf = gauss2d(psfsize / pixelsize, psfsize / pixelsize)
    img = convolve(img.astype(float), psf)
    img = img / img.max() + bg
    #print img.max()
    #img = poisson.rvs(img.flatten() + eps, \
    #                  size=img.size).reshape(imagesize).astype(float)
    return img

# -- unscaled tests --

# def test_unscaled_single(ax, img=None, size=(512,512), npmax=200., \
#                          nbins=40, low_thresh=0., hi_thresh=0.7):
#     if img == None:
#         img = _n.random.uniform(size=size)
#         img*=npmax
#         noisy = add_noise(img, npmax, normalise=False)
#         dif = img - noisy
#     means = _n.zeros(nbins)
#     varis = _n.zeros(nbins)
#     bins = _n.linspace(0., img.max(), num=nbins+1, endpoint=True)
#     for i in range(nbins):
#         bin = (img > bins[i]) & (img <= bins[i+1])
#         means[i] = img[bin].mean()
#         varis[i] = var(dif[bin])
#     scale = means.max()
#     #print "Scale factor", scale
#     ax.plot(means, varis, color='k')

#     valids = ~_n.isnan(varis) & ~_n.isnan(means) & \
#              (means < hi_thresh * img.max()) & (means > low_thresh * img.max())
#     slope = get_slope(means[valids], varis[valids], nbins)
#     np = 1./slope
#     print "    np = %.1f, should be %d" % (np, npmax)
#     fxs = _n.array((0.,img.max()))
#     fys = _n.array((0.,slope * img.max()))
#     ax.plot(fxs, fys, ':', color='k')

# def test_unscaled_double(ax, seed=None, size=(512,512), npmax=200., \
#                          nbins=40, low_thresh=0., hi_thresh=0.7, pc=0.95, \
#                          colour='b'):  
#     if seed == None:
#         seed = _n.random.uniform(size=size)
#     seed *= npmax
#     imga = add_noise(seed, npmax, normalise=False, pc=pc)
#     imgb = add_noise(seed, npmax, normalise=False, pc=pc)
#     dif = imga - imgb    

#     means = _n.zeros(nbins)
#     varis = _n.zeros(nbins)
#     bins = _n.linspace(0., seed.max(), num=nbins+1, endpoint=True)
#     for i in range(nbins):
#         bin = (seed > bins[i]) & (seed <= bins[i+1])
#         means[i] = seed[bin].mean()
#         varis[i] = var(dif[bin])/2.
#     scale = means.max()
#     #print "Scale factor", scale
#     ax.plot(means, varis, color=colour)

#     valids = ~_n.isnan(varis) & ~_n.isnan(means) & \
#              (means < hi_thresh * seed.max()) & \
#              (means > low_thresh * seed.max())
#     slope = get_slope(means[valids], varis[valids], nbins)
#     np = 1./slope
#     print "    np = %.1f, should be %d" % (np, npmax)
#     fxs = _n.array((0.,seed.max()))
#     fys = _n.array((0.,seed.max() * slope))
#     ax.plot(fxs, fys, ':', color=colour)

# def test_unscaled_exp_means(ax, seed=None, size=(512,512), \
#                             npmax=200., nbins=40, low_thresh=0., \
#                             hi_thresh=0.7, pc=0.95, colour='g'):
#     if seed == None:
#         seed = _n.random.uniform(size=size)
#     imga = add_noise(seed, npmax, normalise=False, pc=pc)
#     imgb = add_noise(seed, npmax, normalise=False, pc=pc)
#     mean = (imga + imgb) / 2.
#     dif = imga - imgb

#     means = _n.zeros(nbins)
#     varis = _n.zeros(nbins)
#     bins = _n.linspace(0., mean.max(), num=nbins+1, endpoint=True)
#     for i in range(nbins):
#         bin = (mean > bins[i]) & (mean <= bins[i+1])
#         means[i] = mean[bin].mean()
#         varis[i] = var(dif[bin])/2.
#     scale = means.max()
#     #print "Scale factor", scale
#     ax.plot(means, varis, color=colour)

#     valids = ~_n.isnan(varis) & ~_n.isnan(means) & \
#              (means < hi_thresh * mean.max()) & \
#              (means > low_thresh * mean.max())
#     slope = get_slope(means[valids], varis[valids], nbins)
#     np = 1./slope
#     print "    np = %.1f, should be %d" % (np, npmax)
#     fxs = _n.array((0.,1.))
#     fys = _n.array((0.,slope))
#     ax.plot(fxs, fys, ':', color=colour)

# def test_unscaled_beads(ax, size=(512,512), npmax=200., \
#               nbins=30, low_thresh=0., hi_thresh=0.7):
#     img = build_data(imagesize=size, pixelsize=0.1, psfsize=0.3, \
#                       npmax=npmax, bg=5)
#     test_unscaled_exp_means(ax, seed=img, size=size, npmax=npmax, pc=0.9995, \
#                     nbins=nbins, low_thresh=low_thresh, hi_thresh=hi_thresh, \
#                             colour='r')

# def test_unscaled_shift(ax, seed=None, size=(512,512), npmax=200., \
#                         nbins=30, low_thresh=0., hi_thresh=0.7, \
#                         pc=0.9995, colour='m'):
#     if seed == None:
#         seed = build_data(imagesize=size, pixelsize=0.1, psfsize=0.4, \
#                           npmax=npmax, bg=5, beadsize=0.5)
#         seed = add_noise(seed, npmax, normalise=False, pc=pc)
#         _p.figure(2)
#         _p.imshow(seed)
#     imga = seed.copy()[:,1:]
#     imgb = seed.copy()[:,:-1]
#     mean = (imga + imgb) / 2.
#     dif = _n.abs(imga - imgb)

#     means = _n.zeros(nbins)
#     varis = _n.zeros(nbins)
#     bins = _n.linspace(0., mean.max(), num=nbins+1, endpoint=True)
#     for i in range(nbins):
#         bin = (mean > bins[i]) & (mean <= bins[i+1])
#         means[i] = mean[bin].mean()
#         varis[i] = var(dif[bin])/2.
#     scale = means.max()
#     #print "Scale factor", scale
#     ax.plot(means, varis, color=colour)

#     valids = ~_n.isnan(varis) & ~_n.isnan(means) & \
#              (means < hi_thresh * mean.max()) & \
#              (means > low_thresh * mean.max())
#     slope = get_slope(means[valids], varis[valids], nbins)
#     np = 1./slope
#     print "    np = %.1f, should be %d" % (np, npmax)
#     fxs = _n.array((0.,1.))
#     fys = _n.array((0.,slope))
#     ax.plot(fxs, fys, ':', color=colour)
    
# def run_unscaled_tests(npmax=200):
#     f = _p.figure(1)
#     ax = f.add_subplot(111)
#     ax.set_xlabel('mean (binned)')
#     ax.set_ylabel('variance (binned)')
#     print "Single deviate:"
#     test_unscaled_single(ax, npmax = npmax)
#     print "Double deviate:"
#     test_unscaled_double(ax, npmax = npmax)
#     print "Experimental means:"
#     test_unscaled_exp_means(ax, npmax = npmax)
#     print "Bead-like data:"
#     test_unscaled_beads(ax, npmax = npmax)
#     print "Shift comparison:"
#     test_unscaled_shift(ax, npmax = npmax)
#     _p.draw()

# -- scaled tests ------------------------------------------------------------

def test_single_dev(ax, img=None, size=(512,512), npmax=200., \
                    nbins=40, low_thresh=0., hi_thresh=0.7):
    '''Test (1) one poidev vs means'''
    if img == None:
        img = _n.random.uniform(size=size)
        noisy = add_noise(img, npmax)
        dif = img - noisy

    means = _n.zeros(nbins)
    varis = _n.zeros(nbins)
    bins = _n.linspace(0., 1., num=nbins+1, endpoint=True)
    for i in range(nbins):
        bin = (img > bins[i]) & (img <= bins[i+1])
        means[i] = img[bin].mean()
        varis[i] = var(dif[bin])
    scale = means.max()
    #print "Scale factor", scale
    ax.plot(means, varis, color='k')

    valids = ~_n.isnan(varis) & ~_n.isnan(means) & \
             (means < hi_thresh) & (means > low_thresh)
    slope = get_slope(means[valids], varis[valids], nbins)
    np = 1./slope
    print "    np = %.1f, should be %d" % (np, npmax)
    fxs = _n.array((0.,1.))
    fys = _n.array((0.,slope))
    ax.plot(fxs, fys, ':', color='k')

def test_single_alt(ax, img=None, size=(512,512), npmax=200., \
                    nbins=40, low_thresh=0., hi_thresh=0.7):
    '''Test (1) one poidev vs means'''
    if img == None:
        img = _n.random.uniform(size=size)
        noisy = add_noise(img, npmax)
        dif = _n.abs(img - noisy)

    means = _n.zeros(nbins)
    varis = _n.zeros(nbins)
    bins = _n.linspace(0., 1., num=nbins+1, endpoint=True)
    for i in range(nbins):
        bin = (img > bins[i]) & (img <= bins[i+1])
        means[i] = img[bin].mean()
        varis[i] = var(dif[bin])
    scale = means.max()
    #print "Scale factor", scale
    ax.plot(means, varis, color='y')

    valids = ~_n.isnan(varis) & ~_n.isnan(means) & \
             (means < hi_thresh) & (means > low_thresh)
    slope = get_slope(means[valids], varis[valids], nbins)
    np = 1./slope
    print "    np = %.1f, should be %d" % (np, npmax)
    fxs = _n.array((0.,1.))
    fys = _n.array((0.,slope))
    ax.plot(fxs, fys, ':', color='y')

def test_double_dev(ax, seed=None, size=(512,512), npmax=200., \
                    nbins=40, low_thresh=0., hi_thresh=0.7, pc=0.95, \
                    colour='b'):
    '''test (2) difference of two poidevs vs means'''
    if seed == None:
        seed = _n.random.uniform(size=size)
    imga = add_noise(seed, npmax, pc=pc)
    imgb = add_noise(seed, npmax, pc=pc)
    dif = imga - imgb    

    means = _n.zeros(nbins)
    varis = _n.zeros(nbins)
    bins = _n.linspace(0., 1., num=nbins+1, endpoint=True)
    for i in range(nbins):
        bin = (seed > bins[i]) & (seed <= bins[i+1])
        means[i] = seed[bin].mean()
        varis[i] = var(dif[bin])/2.
    scale = means.max()
    #print "Scale factor", scale
    ax.plot(means, varis, color=colour)

    valids = ~_n.isnan(varis) & ~_n.isnan(means) & \
             (means < hi_thresh) & (means > low_thresh)
    slope = get_slope(means[valids], varis[valids], nbins)
    np = 1./slope
    print "    np = %.1f, should be %d" % (np, npmax)
    fxs = _n.array((0.,1.))
    fys = _n.array((0.,slope))
    ax.plot(fxs, fys, ':', color=colour)

def test_exp_means(ax, seed=None, size=(512,512), npmax=200., \
                    nbins=40, low_thresh=0., hi_thresh=0.7, \
                   pc=0.95, colour='g', plot_fig=None, plot_title='', \
                   useabs=False):
    '''test (2) difference of two poidevs vs means'''
    if seed == None:
        seed = _n.random.uniform(size=size)
    imga = add_noise(seed, npmax, pc=pc)
    imgb = add_noise(seed, npmax, pc=pc)
    mean = (imga + imgb) / 2.
    if not useabs:
        dif = imga - imgb
    else:
        dif = _n.abs(imga - imgb)

    if plot_fig != None:
        plot_dists(imga[406-10:406+10,163-10:163+10], \
                   imgb[406-10:406+10,163-10:163+10], fignum=plot_fig, \
                   tit=plot_title)
    means = _n.zeros(nbins)
    varis = _n.zeros(nbins)
    bins = _n.linspace(0., 1., num=nbins+1, endpoint=True)
    for i in range(nbins):
        bin = (mean > bins[i]) & (mean <= bins[i+1])
        means[i] = mean[bin].mean()
        varis[i] = var(dif[bin])/2.
    scale = means.max()
    #print "Scale factor", scale
    ax.plot(means, varis, color=colour)

    valids = ~_n.isnan(varis) & ~_n.isnan(means) & \
             (means < hi_thresh) & (means > low_thresh)
    slope = get_slope(means[valids], varis[valids], nbins)
    np = 1./slope
    print "    np = %.1f, should be %d" % (np, npmax)
    fxs = _n.array((0.,1.))
    fys = _n.array((0.,slope))
    ax.plot(fxs, fys, ':', color=colour)
    return np

def test_straight_var(ax, seed=None, size=(512,512), npmax=200., \
                      nbins=40, low_thresh=0., hi_thresh=0.7, \
                      pc=0.95, colour='m'):
    '''test (2) difference of two poidevs vs means'''
    if seed == None:
        seed = _n.random.uniform(size=size)
    imga = add_noise(seed, npmax, pc=pc)
    imgb = add_noise(seed, npmax, pc=pc)
    mean = (imga + imgb) / 2.

    means = _n.zeros(nbins)
    varis = _n.zeros(nbins)
    bins = _n.linspace(0., 1., num=nbins+1, endpoint=True)
    for i in range(nbins):
        bin = (mean > bins[i]) & (mean <= bins[i+1])
        means[i] = mean[bin].mean()
        varis[i] = var(_n.concatenate((imga[bin], imgb[bin])))
    scale = means.max()
    #print "Scale factor", scale
    ax.plot(means, varis, color=colour)

    valids = ~_n.isnan(varis) & ~_n.isnan(means) & \
             (means < hi_thresh) & (means > low_thresh)
    slope = get_slope(means[valids], varis[valids], nbins)
    np = 1./slope
    print "    np = %.1f, should be %d" % (np, npmax)
    fxs = _n.array((0.,1.))
    fys = _n.array((0.,slope))
    ax.plot(fxs, fys, ':', color=colour)

def test_beads(ax, seed=None, size=(512,512), npmax=200., \
              nbins=30, low_thresh=0., hi_thresh=0.7):
    if seed == None:
        seed = build_data(imagesize=size, pixelsize=0.1, psfsize=gpsfsize, \
                         npmax=npmax, bg=5, beadsize=0.5)
    return test_exp_means(ax, seed=seed, size=size, npmax=npmax, pc=0.9995, \
                          nbins=nbins, low_thresh=low_thresh, \
                          hi_thresh=hi_thresh, colour='r')

def test_beads_abs(ax, size=(512,512), npmax=200., \
              nbins=30, low_thresh=0., hi_thresh=0.7):
    img = build_data(imagesize=size, pixelsize=0.1, psfsize=gpsfsize, \
                      npmax=npmax, bg=5, beadsize=0.5)
    return test_exp_means(ax, seed=img, size=size, npmax=npmax, pc=0.9995, \
                          nbins=nbins, low_thresh=low_thresh, \
                          hi_thresh=hi_thresh, colour='y', useabs=True)
    

def test_shift(ax, seed=None, size=(512,512), npmax=200., \
               nbins=30, low_thresh=0., hi_thresh=0.7, \
               pc=0.9995, colour='c', do_plot=True):
    if seed == None:
        seed = build_data(imagesize=size, pixelsize=0.1, psfsize=gpsfsize, \
                          npmax=npmax, bg=5, beadsize=0.5)
        seed = add_noise(seed, npmax, normalise=True, pc=pc)
        #_p.figure(2)
        #_p.imshow(seed)
    imga = seed.copy()[:,1:]
    imgb = seed.copy()[:,:-1]
    
    dif = _n.abs(imga - imgb)
    #print "diff max", dif.max()
    mean = (imga + imgb) / 2.

    if do_plot:
        plot_dists(imga[406-10:406+10,163-10:163+10], \
                   imgb[406-10:406+10,163-10:163+10], fignum=4, tit='Shifted')

    bins = _n.linspace(0., 1., num=nbins+1, endpoint=True)
    means = _n.zeros(nbins)
    varis = _n.zeros(nbins)
    for i in range(nbins):
        bin = (mean > bins[i]) & (mean <= bins[i+1])
        means[i] = mean[bin].mean()
        varis[i] = var(dif[bin])/2.

    ax.plot(means, varis, color=colour)
    ax.set_ylabel('variance (binned)')
    ax.set_xlabel('mean (binned)')

    valids = ~_n.isnan(varis) & ~_n.isnan(means) & \
             (means < hi_thresh) & (means > low_thresh)
    slope = get_slope(means[valids], varis[valids], nbins)
    np = 1./slope
    print "    np = %.1f, should be %d" % (np, npmax)
    fxs = _n.array((0.,1.))
    fys = _n.array((0.,slope))
    ax.plot(fxs, fys, ':', color=colour)
    return np

def run_tests(npmax=200.):
    f = _p.figure(1)
    f.clf()
    ax = f.add_subplot(111)
    ax.set_xlabel('mean (binned)')
    ax.set_ylabel('variance (binned)')
    print "Single deviate: (black)"
    test_single_dev(ax, npmax = npmax)
#    print "Single abs diff: (yellow)"
#    test_single_alt(ax, npmax = npmax)
    print "Double deviate: (blue)"
    test_double_dev(ax, npmax = npmax)
#    print "Straight variance: (magenta)" 
#    test_straight_var(ax, npmax = npmax)
    print "Experimental means: (green)"
    test_exp_means(ax, npmax = npmax)
    print "Bead-like data: (red)"
    test_beads(ax, npmax = npmax)
    print "Bead-like data, abs. diff.: (yellow)"
    test_beads_abs(ax, npmax = npmax)
    print "Shift comparison: (cyan)"
    test_shift(ax, npmax = npmax)
    _p.draw()

def run_comparison(npmax=200):
    f = _p.figure(1)
    f.clf()
    ax = f.add_subplot(111)
    ax.set_xlabel('mean (binned)')
    ax.set_ylabel('variance (binned)')
    print "Bead-like data, abs. diff.: (yellow)"
    test_beads_abs(ax, npmax = npmax)
    print "Shift comparison: (cyan)"
    test_shift(ax, npmax = npmax, do_plot = False)
    _p.draw()

def make_smoothed_real(filename, smoothsize=5):
    print "Loading %s..." % filename
    img = _n.array(Image.open(filename), dtype=float)
    return median_filter(img, size=smoothsize)

def test_accuracies():
    repeats = 5
    f0 = _p.figure(1)
    f0.clf()
    ax0 = f0.add_subplot(111)
    ax0.set_ylabel('accuracy')
    ax0.set_xlabel('np')
    f1 = _p.figure(2)
    ax1 = f1.add_subplot(111)
    ax1.set_ylabel('variance (binned)')
    ax1.set_xlabel('mean (binned)')

    nps = (_n.arange(16,dtype=float) + 1) * 25
    compaccuracies = _n.empty((len(nps), repeats))
    shiftaccuracies = _n.empty((len(nps), repeats))
    real = make_smoothed_real('urr/URR051.TIF')
    pc = 0.9995
    for j in range(repeats):
        for i, np in enumerate(nps):
            #beads = build_data(imagesize=(512,512), pixelsize=0.1, \
            #                   psfsize=gpsfsize, npmax=np, \
            #                   bg=5, beadsize=0.5)
            nptwo = test_beads(ax1, npmax=np, seed=real)
            noisy = add_noise(real, np, normalise=True, pc=pc)
            npshift = test_shift(ax1, npmax=np, seed=noisy, do_plot=False)
            compaccuracies[i,j] = nptwo / np
            shiftaccuracies[i,j] = npshift / np
            
    compmeans = compaccuracies.mean(axis=1)
    compsems = sem(compaccuracies, axis=1)
    shiftmeans = shiftaccuracies.mean(axis=1)
    shiftsems = sem(shiftaccuracies, axis=1)
    ax0.errorbar(nps, compmeans, yerr=compsems, fmt='b', label='comparison')
    ax0.errorbar(nps, shiftmeans, yerr=shiftsems, fmt='r', label='shift')

# -- analysis calculations --

def calculate_robustmax(img, pc=0.95):
    lim = img.size * pc
    return img.flatten()[img.argsort(None)][lim:].mean()

def analyse_noise(imga=None, imgb=None, noise=400., noisemax=400., \
                  low_thresh=0., hi_thresh=0.75, nbins=20, \
                  colour=None, show_np=False, smidgen=None, pc=0.999, \
                  ax=None, fignum=1, verbose=False):
    if imga == None:
        lenx = get_img()
        imga = add_noise(lenx, noise)
        imgb = add_noise(lenx, noise)
    else:
        if type(imga) == str:
            print "Loading", imga
            imga = _n.array(Image.open(imga),dtype=float)
        elif type(imga) == _n.ndarray:
            imga = _n.asarray(imga, dtype=float)
        if type(imgb) == str:
            print "Loading", imgb
            imgb = _n.array(Image.open(imgb),dtype=float)
        elif type(imgb) == _n.ndarray:
            imgb = _n.asarray(imgb, dtype=float)
        elif imgb == None:
            if verbose: print "Shifting image a..."
            imgb = imga.copy()[:,:-1]
            imga = imga.copy()[:,1:]
        
    imga = imga.copy() / calculate_robustmax(imga, pc=pc)
    #print "imga shape", imga.shape
    imgb = imgb.copy() / calculate_robustmax(imgb, pc=pc)
    diff = imga - imgb
    #print "diff max", diff.max()
    mean = (imga + imgb) / 2.

    if colour == None:
        colour = mpl.cm.Blues(noise/noisemax + 0.2)
    
    bins = _n.linspace(0., 1., nbins + 1, endpoint=True)
    means = _n.zeros(nbins)
    varis = _n.zeros(nbins)
    for i in range(nbins):
        bin = (mean > bins[i]) & (mean <= bins[i+1])
        means[i] = mean[bin].mean()
        varis[i] = var(diff[bin])/2. # /2 to account for -ing 2 poisson dists
    if ax == None:
        f = _p.figure(fignum)
        ax = f.add_subplot(111)
    ax.plot(means, varis, color=colour)
    ax.set_ylabel('variance (binned)')
    ax.set_xlabel('mean (binned)')

    valids = ~_n.isnan(varis) & ~_n.isnan(means) & \
             (means < hi_thresh) & (means > low_thresh)
    slope = get_slope(means[valids], varis[valids], nbins)
    np = (1./slope)
    #    print "slope = %.3g" %  slope
    if verbose: print "np= %.1f" % np
    fxs = _n.array((0., 1.))
    fys = _n.array((0., slope))
    ax.plot(fxs, fys, ':', color=colour, alpha=0.75)
    if show_np:
        smidgen = 0.002
        ax.text(0.85, slope * 0.85 + smidgen, \
                    "%.1f" % np, color=colour)
    return np

def plot_dists(a=None, b=None, fignum=3, lng=1000, tit=''):
    if a == None:
        a = _n.random.normal(loc=1.0, size=lng)
        b = _n.random.normal(loc=1.0, size=lng)
    else:
        lng = a.size
    if _n.rank(a) > 1:
        a = a.flatten()
    if _n.rank(b) > 1:
        b = b.flatten()
    dif = a - b
    _p.figure(fignum)
    _p.subplot(311)
    _p.title(tit)
    lw = 1.5
    _p.vlines(_n.arange(lng)[dif > 0], a[dif > 0], b[dif > 0], \
              color='b', linewidth=lw)
    _p.vlines(_n.arange(lng)[dif <= 0], a[dif <= 0], b[dif <= 0], \
              color='r', linewidth=lw)
    _p.subplot(312)
    _p.vlines(_n.arange(lng)[dif > 0], 0., dif[dif > 0], \
              color='b', linewidth=lw)
    _p.vlines(_n.arange(lng)[dif <= 0], 0., dif[dif <= 0], \
              color='r', linewidth=lw)
    _p.subplot(313)
    _p.vlines(_n.arange(lng)[dif > 0], 0., _n.abs(dif)[dif > 0], \
              color='b', linewidth=lw)
    _p.vlines(_n.arange(lng)[dif <= 0], 0., _n.abs(dif)[dif <= 0], \
              color='r', linewidth=lw)
    
def simple_compare(img, noise_levels=None, pv=None, bg=0.01, \
                   smooth='gauss', smsz=10):
    
    if type(img) == str:
        print "Loading", img
        img = _n.array(Image.open(img),dtype=float)
    if noise_levels == None:
        noise_levels=[100,200,400,800]
    if pv == None:
        pv = pyvis()
    pv.AddImg(img, 'original')
    
    sh = img.shape
    if smooth == 'gauss':
        k = gauss2d(smsz,smsz)
        smoothed = convolve(img, k).astype(float)
    elif smooth == 'median':
        smoothed = median_filter(img, size=smsz).astype(float)
    # ranges between mini and maxi
    smoothed -= bg
    # ranges between mini-bg and maxi-bg
    pv.AddImg(smoothed, 'smoothed')
    #noise_level = noise_levels[0]
    noisies = []
    for i, noise_level in enumerate(noise_levels):
        scaled = _n.clip(smoothed.flatten() / smoothed.max() * noise_level, \
                         1e-20, noise_level)
        noisy = poisson.rvs(scaled, size=scaled.size).reshape(sh).astype(float)
        pv.AddImg(noisy, str(noise_level))
        noisies.append( noisy )
    return noisies, pv
