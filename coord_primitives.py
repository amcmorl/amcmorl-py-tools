import numpy as np

def sphere(radius=1.0, npts=(100,100)):
    '''
    Notes
    -----
    Shamelessly borrowed from FPs test_mesh_sphere replacement'''
    np_phi = npts[0]*1j
    np_theta = npts[1]*1j    
    phi, theta = np.mgrid[0:2*np.pi:np_phi,0:np.pi:np_theta]
    x = radius * np.sin(theta) * np.cos(phi)
    y = radius * np.sin(theta) * np.sin(phi)
    z = radius * np.cos(theta)
    return x, y, z
