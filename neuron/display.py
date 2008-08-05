import matplotlib
matplotlib.use('PDF') # or try PDF
import numpy as n
import neuron.path as np
import pylab as p
ptsa = np.read_points('P0-2/040624/global_coords.txt')
lda = np.read_ldens('P0-2/040624/ldens_nspines.txt')
ptsb = np.read_points('P9-11/050822/global_coords.txt')
ldb = np.read_ldens('P9-11/050822/ldens_nspines.txt')
xa, ya = ptsa[1][:,0:2].transpose()
print "A min:", xa.min(), "A max:", xa.max(),
xb, yb = ptsb[1][:,0:2].transpose()
print "B min:", xb.min(), "B max:", xb.max()
xb += -xb.min() + xa.max()
yb -= 100
xt = n.hstack((xa,xb))
yt = n.hstack((ya,yb))
aspect = (yt.max() - yt.min())/(xt.max() - xt.min())
f = p.figure(figsize=(8, 8*aspect), facecolor='w')
ax = f.add_axes([0,0,1,1])
ca = matplotlib.cm.jet(lda)
cb = matplotlib.cm.jet(ldb)
for tx,ty,tc,ld in zip(xa,ya,ca,lda):
    ax.plot((tx,),(ty,),'o',color=tc, markersize=ld*8+0.5, mew=0.1)
for tx,ty,tc,ld in zip(xb,yb,cb,ldb):
    ax.plot((tx,),(ty,),'o',color=tc, markersize=ld*8+0.5, mew=0.1)
ax.set_xlim(xt.min() * 1.1, xt.max() * 1.1)
ax.set_ylim(yt.min() * 1.1 * aspect, yt.max() * 1.1 * aspect)
ax.set_xticks([])
ax.set_yticks([])
f.savefig('ldens_both')
print "Done"
