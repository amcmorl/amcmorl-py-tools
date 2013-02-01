import numpy as np
from amcmorl_py_tools.vecgeom.rotations import rotate_by_angles

def plot_pts(pts=None, mu=None):
    p3d = mlab.pipeline

    # plot unit sphere
    sphere_pts = coord_primitives.sphere()
    sphere1 = mlab.mesh(*sphere_pts)
    sphere1.actor.mapper.scalar_visibility = False
    sphere1.actor.property.color = (241/255., 233/255., 199/255.)
    sphere1.actor.property.opacity = 0.25

    # plot pts
    if pts != None:
        x = pts[...,0]
        y = pts[...,1]
        z = pts[...,2]
        src = p3d.scalar_scatter(x,y,z)
        glyphs = p3d.glyph(src)
        glyphs.actor.property.color = (0.0, 0.0, 1.0)
        glyphs.glyph.glyph_source.glyph_source = \
            glyphs.glyph.glyph_source.glyph_list[4]
        glyphs.glyph.glyph.scale_factor = 0.1

    # plot mean
    if mu != None:
        zeros = np.zeros_like(mu[0])[..., np.newaxis]
        mlab.quiver3d(zeros, zeros, zeros,
                      mu[0, np.newaxis], mu[1, np.newaxis], mu[2, np.newaxis])

def generate_cone_circle(theta, phi, angle, resolution=50.):
    x, y, z = 0, 1, 2
    phi_prime = np.linspace(0, 2 * np.pi, resolution)
    # start with defining cone around (0,0,1)
    origin = np.array((0.,0.,1.))
    P_j = np.array([rotate_by_angles(origin, angle, q)
                    for q in phi_prime]).T
    return rotate_by_angles(P_j, theta, phi).T
        
#~ def plot_circle(mu, angle, scalars=None, scalar_max=None,
        #~ color=None, radius=0.01, alpha=1.,
        #~ resolution=50.):
    #~ x, y, z = 0, 1, 2
    #~ theta, phi = convert_cartesian_to_polar(mu)
    #~ P_k = generate_cone_circle(theta, phi, angle, resolution).T
    #~ 
    #~ p3d = mlab.pipeline
    #~ if scalars == None:
        #~ tube = p3d.tube(p3d.line_source(P_k[x], P_k[y], P_k[z]))
    #~ else:
    #~ if type(scalars) == type(1) or type(scalars == type(1.)):
        #~ #'arraying %d' % (
        #~ scalars *= np.ones_like(P_k[0])
        #~ #'assigning scalars 
        #~ tube = p3d.tube(p3d.line_source(P_k[x],
                        #~ P_k[y],
                        #~ P_k[z], scalars))
    #~ tube.filter.radius = radius
    #~ surf = p3d.surface(tube)
    #~ if color != None:
    #~ surf.actor.actor.property.color = color
    #~ surf.actor.actor.property.opacity = alpha
    #~ if scalar_max != None:
    #~ mm = surf.module_manager.scalar_lut_manager
    #~ mm.data_range = np.array([0, scalar_max])
