import numpy as n, ncdf, os.path as path, sys
from Scientific.IO.NetCDF import NetCDFFile
from glob import glob
import neuron.path as npth, neuron.tile as ntl
import pylab as p


def isin(start, stop, pt):
    if n.all(pt > start) and n.all(pt < stop):
        #print "isin..."
        #print "start", start
        #print "stop", stop
        #print "pt", pt
        return True
    else:
        return False

def get_dir_info(cdir):
    if cdir[-1] == '/':
        cdir = cdir[:-1]
    f = open(cdir + '/info.txt')
    lines = f.readlines()
    subdirs = (lines[1].strip() == 'subdirs')
    compat = (lines[3].strip() == 'compat')
    f.close()
    return subdirs, lines[2].strip(), compat

def container_stack(global_offsets, sizes, ptpx):
    beststack = ''
    beststackdist = 1e12
    beststart = None
    for k,v in global_offsets.iteritems():
        # calculate upper and lower limits in microns
        start = v
        stop  = (sizes[k] + v)
        #print "Stack: %4s start = %18s stop = %s" % \
        #(k, str(start), str(stop)),
        #print "stack name:", k
        if isin(start, stop, ptpx):
            centre = (start + stop) / 2.
            dist = n.linalg.norm(centre - ptpx)
            if dist < beststackdist:
                beststack = k
                beststart = start
                beststackdist = n.linalg.norm(centre - ptpx)
                #print "*"
    if not (beststart == None):
        return beststack, beststart, sizes[beststack]
    else:
        raise UserWarning('No stack found to contain this point.')

    '''extract nspines randomly-selected spines from cell in celldir'''
    # load stack information for later
    daydirs, ending, compat = get_dir_info(celldir)
    stacks = ntl.get_stacks(celldir, daydirs=daydirs, ending=ending)
    sizes = {}
    for stack in stacks:
        #   - get stack sizes
        stackname = path.splitext(path.split(stack)[-1])[0].split('_')[0]
        sz = ncdf.get_size(stack, compat=compat)
        sizes[stackname] = sz
    global_offsets = ntl.calculate_global_offsets( \
        celldir + 'stack_offsets.txt')
        
    coordsfile = celldir + 'global_coords.txt'
    res = npth.read_resolution(coordsfile)
    pts = npth.read_points(coordsfile)
    nsps = (pts[0] == 'p').sum() + (pts[0] == 'f').sum()
    if nspines > nsps:
        nspines = nsps

    sps = n.arange(nsps)
    n.random.shuffle(sps)
    spinestack = []
    i = 0
    #j = 0
    while (i < nspines):# and (j < attempts):
        # pts, 
        isspine = (pts[0] == 'p') | (pts[0] == 'f')
        spinepts = npth.sub_points(pts, isspine)
        # pick a spine
        spine_loc = spinepts[1][sps[i]] / res
        #print '-------------'
        #print "spine_loc", ','.join(['%5.1f' % (x) for x in spine_loc])

        # - identify which stack(s) it's in
        stack, start, size = container_stack(global_offsets, sizes, spine_loc)
        
        local_px = spine_loc - start
        if daydirs:
            stackfile = glob(celldir + '*/' + stack + ending)[0]
        else:
            stackfile = celldir + stack + ending

        # - extract volume
        ncfile = NetCDFFile(stackfile, 'r')
        var = ncfile.variables['var0']
        ori = n.zeros(3)
        spstart = n.vstack((ori, local_px - delta)).max(0).astype(int)
        spstop = n.vstack((size, local_px + delta)).min(0).astype(int)
        #print "spstart:", spstart
        #print "spstop: ", spstop
        #print "compat?:", compat
        if compat:
            data = n.asarray(var[spstart[2] : spstop[2], \
                                 spstart[1] : spstop[1], \
                                 spstart[0] : spstop[0]])
        else:
            data = n.asarray(var[spstart[0] : spstop[0], \
                                 spstart[1] : spstop[1], \
                                 spstart[2] : spstop[2]])
        if compat:
            data = data.transpose()
        ncfile.close()
        spinestack.append(data)
        i += 1
    return spinestack, res

class SpineSorter:
    def __init__(self, celldir, nspines=1, delta=20, attempts=500):
        '''extract nspines randomly-selected spines from cell in celldir'''
        # load stack information for later
        daydirs, ending, compat = get_dir_info(celldir)
        stacks = ntl.get_stacks(celldir, daydirs=daydirs, ending=ending)
        sizes = {}
        for stack in stacks:
            #   - get stack sizes
            stackname = path.splitext(path.split(stack)[-1])[0].split('_')[0]
            sz = ncdf.get_size(stack, compat=compat)
            sizes[stackname] = sz
        global_offsets = ntl.calculate_global_offsets( \
            celldir + 'stack_offsets.txt')

        coordsfile = celldir + 'global_coords.txt'
        self.res = npth.read_resolution(coordsfile)
        pts = npth.read_points(coordsfile)
        nsps = (pts[0] == 'p').sum() + (pts[0] == 'f').sum()
        if nspines > nsps:
            nspines = nsps

        sps = n.arange(nsps)
        n.random.shuffle(sps)
        self.spinestack = []
        i = 0
        while (i < nspines):
            isspine = (pts[0] == 'p') | (pts[0] == 'f')
            spinepts = npth.sub_points(pts, isspine)
            # pick a spine
            spine_loc = spinepts[1][sps[i]] / self.res
            #print '-------------'
            #print "spine_loc", ','.join(['%5.1f' % (x) for x in spine_loc])

            # - identify which stack(s) it's in
            stack, start, size = container_stack(global_offsets, sizes, spine_loc)

            local_px = spine_loc - start
            if daydirs:
                stackfile = glob(celldir + '*/' + stack + ending)[0]
            else:
                stackfile = celldir + stack + ending

            # - extract volume
            ncfile = NetCDFFile(stackfile, 'r')
            var = ncfile.variables['var0']
            ori = n.zeros(3)
            spstart = n.vstack((ori, local_px - delta)).max(0).astype(int)
            spstop = n.vstack((size, local_px + delta)).min(0).astype(int)
            #print "spstart:", spstart
            #print "spstop: ", spstop
            #print "compat?:", compat
            if compat:
                data = n.asarray(var[spstart[2] : spstop[2], \
                                     spstart[1] : spstop[1], \
                                     spstart[0] : spstop[0]])
            else:
                data = n.asarray(var[spstart[0] : spstop[0], \
                                     spstart[1] : spstop[1], \
                                     spstart[2] : spstop[2]])
            if compat:
                data = data.transpose()
            ncfile.close()
            self.spinestack.append(data)
            i += 1
    
    def display_spines(self):
        from enthought.mayavi.tools import mlab
        from enthought.mayavi.sources.array_source import ArraySource
        from enthought.mayavi.modules.api import IsoSurface, Outline
        self.f3d = mlab.figure()
        self.f3d.scene.z_plus_view()
        self.src = ArraySource()
        self.src.spacing = self.res
        self.src._update_image_data_fired()
        self.src.scalar_data = self.spinestack[0]
        self.f3d.add_child(self.src)
        oline = Outline()
        iso = IsoSurface()
        self.src.add_child(oline)
        self.src.add_child(iso)
        iso.actor.mapper.scalar_visibility = False
        iso.actor.property.color = (0., 0.5, 1.0)
        
        self.f2d = p.figure(figsize=(2,2))
        self.toolbar = self.f2d.canvas.toolbar
        self.cid = self.f2d.canvas.mpl_connect( \
            'button_press_event', self.click_handler)
        ax = self.f2d.add_subplot(111)
        ax.imshow(self.spinestack[0].max(2))
        ax.set_xticks([])
        ax.set_yticks([])
        p.show()
        self.i = 0

    def save_current(self, name):
        f2d_save = p.figure(2, figsize=(6,6), facecolor='w')
        ax_save = f2d_save.add_axes((0.,0.,1.,1.))
        ax_save.imshow(self.spinestack[self.i].max(2))
        ax_save.set_xticks([])
        ax_save.set_yticks([])
        name2d = path.splitext(name)[0] + '_2d' + path.splitext(name)[1]
        name3d = path.splitext(name)[0] + '_3d' + path.splitext(name)[1]
        f2d_save.savefig(name2d)
        p.close(f2d_save)
        self.f3d.scene.save_png(name3d)

    def click_handler(self, evt):
        if self.toolbar.mode == "":
            if not evt.inaxes:
                del self.spinestack[self.i]
            else:
                self.i += 1
        if self.i >= (len(self.spinestack) - 1):
            print "Finished"
            self.f2d.canvas.mpl_disconnect(self.cid)
        else:
            self.src.scalar_data = self.spinestack[self.i]
            self.f2d.axes[0].imshow(self.spinestack[self.i].max(2))
            p.draw()

def save_vol(vol, name, ss):
    f2d_save = p.figure(2, figsize=(6,6), facecolor='w')
    ax_save = f2d_save.add_axes((0.,0.,1.,1.))
    ax_save.imshow(vol.max(2))
    ax_save.set_xticks([])
    ax_save.set_yticks([])
    name2d = path.splitext(name)[0] + '_2d' + path.splitext(name)[1]
    name3d = path.splitext(name)[0] + '_3d' + path.splitext(name)[1]
    f2d_save.savefig(name2d)
    p.close(f2d_save)
    ss.src.scalar_data = vol
    ss.f3d.scene.save_png(name3d)
