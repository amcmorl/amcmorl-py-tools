import numpy as np

def sphere(radius=1.0, npts=(100,100)):
    '''shamelessly borrowed from FPs test_mesh_sphere replacement'''
    pi = np.pi
    sin = np.sin
    cos = np.cos
    np_phi = npts[0]*1j
    np_theta = npts[1]*1j
    phi, theta = np.mgrid[0:pi:np_phi,0:2*pi:np_theta]
    x = radius * sin(phi) * cos(theta)
    y = radius * sin(phi) * sin(theta)
    z = radius * cos(phi)
    return x, y, z
