import numpy as np

def sphere(radius=1.0, npts=(100,100)):
    '''shamelessly borrowed from FPs test_mesh_sphere replacement'''
    np_phi = npts[0]*1j
    np_theta = npts[1]*1j
    phi, theta = np.mgrid[0:np.pi:np_phi,0:2*np.pi:np_theta]
    x = radius * np.sin(phi) * np.cos(theta)
    y = radius * np.sin(phi) * np.sin(theta)
    z = radius * np.cos(phi)
    return x, y, z
