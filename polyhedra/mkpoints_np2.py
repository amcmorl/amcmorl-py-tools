#!/usr/bin/env python

# Distribute a set of points randomly across a sphere, allow them
# to mutually repel and find equilibrium.

import sys
import string
import random
from math import pi, asin, atan2, cos, sin, sqrt
import numpy as np
from enthought.mayavi import mlab
import enthought.mayavi.tools.pipeline as mvp
#mlab.options.backend = 'envisage'

def generate_points(n):
    theta = np.random.random(size=(n)) * 2 * np.pi
    phi = np.arcsin(np.random.random(size=(n)) * 2 - 1)
    pt = np.array((np.cos(theta) * np.cos(phi),
                       np.sin(theta) * np.cos(phi),
                       np.sin(phi))).transpose()
    return pt

def mkpoints(n=20, pt=None, graphic=False):
    # Invent a randomly distributed point.
    #
    # To distribute points uniformly over a spherical surface, the
    # easiest thing is to invent its location in polar coordinates.
    # Obviously theta (longitude) must be chosen uniformly from
    # [0,2*pi]; phi (latitude) must be chosen in such a way that
    # the probability of falling within any given band of latitudes
    # must be proportional to the total surface area within that
    # band. In other words, the probability _density_ function at
    # any value of phi must be proportional to the circumference of
    # the circle around the sphere at that latitude. This in turn
    # is proportional to the radius out from the sphere at that
    # latitude, i.e. cos(phi). Hence the cumulative probability
    # should be proportional to the integral of that, i.e. sin(phi)
    # - and since we know the cumulative probability needs to be
    # zero at -pi/2 and 1 at +pi/2, this tells us it has to be
    # (1+sin(phi))/2.
    #
    # Given an arbitrary cumulative probability function, we can
    # select a number from the represented probability distribution
    # by taking a uniform number in [0,1] and applying the inverse
    # of the function. In this case, this means we take a number X
    # in [0,1], scale and translate it to obtain 2X-1, and take the
    # inverse sine. Conveniently, asin() does the Right Thing in
    # that it maps [-1,+1] into [-pi/2,pi/2].

    # For the moment, my repulsion function will be simple
    # inverse-square, followed by a normalisation step in which we pull
    # each point back to the surface of the sphere.
    if pt == None:
        theta = np.random.random(size=(n)) * 2 * np.pi
        phi = np.arcsin(np.random.random(size=(n)) * 2 - 1)
        pt = np.array((np.cos(theta) * np.cos(phi),
                       np.sin(theta) * np.cos(phi),
                       np.sin(phi))).transpose()
        
    if graphic:
        mvpts = mvp.scalarscatter(pt[:,0], pt[:,1], pt[:,2])
        glyph = mvp.glyph(mvpts)
        surf = mvp.surface(mvp.delaunay3d(mvpts))
        surf.actor.property.representation = 'wireframe'
        surf.actor.property.color = (0,0,1.)
        glyph.actor.actor.property.color = (0,0,1.)
        glyph.glyph.glyph.scale_factor = 0.25
        
    z = 0
    while 1:
        z += 1
        # Determine the total force acting on each point.
        dv = (pt[np.newaxis,...] - pt[:,np.newaxis,:]) # n x n x 3
        dl = np.sqrt(np.sum(dv ** 2, axis=-1)) # n x n
        dl3 = dl ** 3
        fv = dv / dl3[...,np.newaxis]  # n x n x 3
        fv[np.isnan(fv)] = 0
        f = np.sum(fv, axis=0) # n x 3
        if z < 3:
            print "first f", f[0]
        ftotal = np.sum(np.sqrt(np.sum(f ** 2, axis=-1))) # scalar

        # Scale the forces to ensure the points do not move too far in
        # one go. Otherwise there will be chaotic jumping around and
        # never any convergence.
        if ftotal > 0.25:
            fscale = 0.25 / ftotal
        else:
            fscale = 1

        # Move each point, and normalize. While we do this, also track
        # the distance each point ends up moving.
        p2 = pt + f * fscale # n x 3
        pl = np.sqrt(np.sum(p2 ** 2, axis=-1)) # n
        p2 /= pl[...,np.newaxis] # n x 3
        dv = pt - p2 # n x 3
        dl = np.sqrt(np.sum(dv ** 2, axis=-1)) # n
        dist = np.sum(dl) # scalar
        pt = p2
        if graphic:
            glyph.module_manager.source.data.points = pt
        
        # Done. Check for convergence and finish.
        sys.stderr.write(str(dist) + "\n")
        if dist < 1e-2:
            break
    return pt

def realprint(a):
    for i in range(len(a)):
	outfile.write(str(a[i]))
	if i < len(a)-1:
	    outfile.write(" ")
	else:
	    outfile.write("\n")

def pprint(*a):
    realprint(a)

# for use from command line
if __name__ == "__main__":
    args = sys.argv[1:]

    if len(args) > 0:
        n = string.atoi(sys.argv[1])
        args = args[1:]
    else:
        n = 7

    if len(args) > 0:
        outfile = open(args[0], "w")
        args = args[1:]
    else:
        outfile = sys.stdout

    points = mkpoints(n)
    # Output the points.
    for x, y, z in points:
        pprint(x, y, z)
