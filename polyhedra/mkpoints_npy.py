#!/usr/bin/env python

# Distribute a set of points randomly across a sphere, allow them
# to mutually repel and find equilibrium.

import sys
import string
import random
from math import pi, asin, atan2, cos, sin, sqrt
import numpy as np

def mkpoints(n=20, pt=None):
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
        points = np.array((np.cos(theta) * np.cos(phi),
                           np.sin(theta) * np.cos(phi),
                           np.sin(phi))).transpose()
    else:
        points = pt.copy()

    z = 0
    while 1:
        z += 1
        # Determine the total force acting on each point.
        forces = np.zeros_like(points)
        for i in range(n):
            p = points[i]
            f = np.zeros((3,))
            ftotal = 0
            for j in range(n):
                if j == i: continue
                q = points[j]

                # Find the distance vector, and its length.
                dv = p - q
                dl = np.sqrt(np.sum(dv**2))

                # The force vector is dv divided by dl^3. (We divide by
                # dl once to make dv a unit vector, then by dl^2 to
                # make its length correspond to the force.)
                dl3 = dl ** 3
                fv = dv / dl3

                # Add to the total force on the point p.
                f += fv

            # Stick this in the forces array.
            forces[i] = f

            # Add to the running sum of the total forces/distances.
            ftotal += np.sqrt(np.sum(f**2))

        if z < 3:
            print "first f", forces[0]
        # Scale the forces to ensure the points do not move too far in
        # one go. Otherwise there will be chaotic jumping around and
        # never any convergence.
        if ftotal > 0.25:
            fscale = 0.25 / ftotal
        else:
            fscale = 1

        # Move each point, and normalise. While we do this, also track
        # the distance each point ends up moving.
        dist = 0
        for i in range(n):
            p = points[i]
            f = forces[i]
            p2 = p + f * fscale
            pl = np.sqrt(np.sum(p2**2))
            p2 /= pl
            dv = p - p2
            dl = np.sqrt(np.sum(dv**2))
            dist += dl
            points[i] = p2

        # Done. Check for convergence and finish.
        sys.stderr.write(str(dist) + "\n")
        if dist < 1e-6:
            break
    return points

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
