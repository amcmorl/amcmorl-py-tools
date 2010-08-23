from curves import *
test = np.testing

def test_kerr2size():
    a, err = (1., 0.1)
    r = 2 * np.sqrt(np.log(1. / 0.1))
    test.assert_almost_equal(kerr2size(a,err), r)

def test_k2fwhm_fwhm2k():
    a = 1
    r = 2 * a * np.sqrt(np.log(2))
    test.assert_almost_equal(k2fwhm(a), r)
    test.assert_almost_equal(fwhm2k(r), a)
    test.assert_almost_equal(a, fwhm2k(k2fwhm(a)))

def test_gauss1d():
    x = np.linspace(0, 2, 101)
    fwhm = 0.5
    k = fwhm2k(fwhm)
    gs = gauss1d(x, 1., fwhm, 1.)
    test.assert_almost_equal(gs[50], 1.)
    test.assert_almost_equal(gs[0], 0., decimal=4)
    test.assert_almost_equal(gs[-1], 0., decimal=4)
    
def test_dbl_boltzman():
    x = np.linspace(0, 2, 101)
    t1, t2 = 0.6, 1.4
    tau, h = 0.05, 1
    y = dbl_boltzman(x, t1, tau, t2, tau, h)
    test.assert_almost_equal(y[0], 0, decimal=4)
    test.assert_almost_equal(y[-1], 0, decimal=4)
    ymax = np.max(y)
    test.assert_equal(y[50], ymax)
    hymax = ymax / 2.
    test.assert_equal(np.argmin(np.abs(y[:50] - hymax)), 30)
    test.assert_equal(np.argmin(np.abs(y[50:] - hymax)), 20)
    
    
    
