import scipy.io
import multiregress

def test_multiregress():
    data = scipy.io.read_array('/home/amcmorl/lib/biometry/'.\
                               'box16.1_air_pollution.txt',
                               columns=(1,-1))
    Y = data[...,0]
    X = data[...,1:3]
    coefs, ses = multiregress(X,Y)
    supplied_coefs = np.array((77.231, -1.0480, 0.0243))
    assert np.all((coefs - supplied_coefs) < 1e-2)
    supplied_ses = np.array((0.37327, 0.0047880))
    assert np.all((ses - supplied_ses) < 1e-4)
    print "multiregress ok"
