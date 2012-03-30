import coord_primitives
from mayavi import mlab

def draw_sphere(color=(241/255., 233/255., 199/255.),
                alpha=0.25, radius=1.):
    # plot unit sphere
    sphere_pts = coord_primitives.sphere()
    sphere_pts = tuple([radius * x for x in sphere_pts])
    sphere1 = mlab.mesh(*sphere_pts)
    sphere1.actor.mapper.scalar_visibility = False
    sphere1.actor.property.color = color
    sphere1.actor.property.opacity = alpha

def animate_scene(n_frames = 36, step = 10, fbase='movie',
                  size=(700.,700.), save=True):
    sc = mlab.gcf()
    #sc.scene.set_size(size)
    #sc.scene.x_minus_view()
    #sc.scene.camera.position = np.array([-1., 0., 0.])
    for i in xrange(n_frames):
        sc.scene.camera.azimuth(step)
        if save:
            mlab.savefig(fbase + '%03d.jpg' % i)
        else:
            mlab.show()
