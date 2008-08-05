import numpy as n
import ncdf
import os.path

from vectors import perpz
from numpy.linalg import norm
import coordhandling
from neuron.point import fwhm_diam, trace_one, get_cnrs
from neuron import pdbg
import neuron
from neuron.tile import calculate_global_offsets
import wx
import matplotlib
matplotlib.interactive(False)
matplotlib.use("WXAgg")
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg, \
     FigureManager, NavigationToolbar2WxAgg
from matplotlib.figure import Figure
from matplotlib.axes import Subplot
import os.path

neuron.debug = 'debug terse'

def trace_stack(stack, pixel_spacing, offsets_file):
    offsets = calculate_global_offsets(offsets_file)
    stack_code = stack.split('/')[-1].split('_')[0]
    try:
        offset = offsets[stack_code]
    except KeyError:
        print "Offset not found - defaulting to 0"
        offset = [0,0,0]
    tl = ThreadMaster(stack, write_thread, cbargs=[pixel_spacing, offset])
    

def get_last_line_num(f):
    fread = open(f, 'r')
    lines = fread.readlines()
    lines.reverse()
    line_num = 0
    for line in lines:
        if not line[0] == '#' and not line.strip() == '':
            line_num = int(line.split()[0])
            break
    fread.close()
    return line_num



def write_thread(stack_name, thread, pxspc, offset):
    '''write thread information to file'''
    stack_code = stack_name.split('/')[-1].split('_')[0]
    foutn = stack_code + '_coords.txt'
    # thread is in format x,y,z,di,code
    print "Writing file:", foutn
    #print "with info:", pxspc, offset
    if os.path.isfile(foutn):
        line_num_offset = get_last_line_num(foutn) + 1
        fout = open(foutn, 'a')
    else:
        line_num_offset = 0
        fout = open(foutn, 'w')
    # add stack offset
    thread[:,0:3] += offset
    # convert xyz pixels into microns
    thread[:,0:3] *= n.array(pxspc)
    # convert diam into microns - assumes x and y resolution are =
    thread[:,3] *= pxspc[0]
    # assign parents
    npts = thread.shape[0]
    for i in xrange(npts):
        code = str(thread[i,-1])
        if code == 's':
            parent = -1
        elif code in ['c', 'p', 'f', 't', 'b']:
            parent = 0
        else:
            parent = -2
        xyzd = ' '.join(['%7.3f' % x for x in list(thread[i,0:4])])
        outstr = ('%04d' % (i + line_num_offset)) + ' ' \
                 + code + ' ' + xyzd + '  %04d' % parent + '\n'
        #print outstr.strip()
        fout.write(outstr)
    fout.close
    
    
def encode_pointtype(fullname):
    '''convert full name of point type into one char code'''
    pointtypes = {'soma': 's', \
                  'continuation':'c', \
                  'spine':'p', \
                  'filopodium':'f', \
                  'branch':'b', \
                  'join':'j', \
                  'terminus':'t'}
    if pointtypes.has_key(fullname):
        return pointtypes[fullname]
    else:
        return None
         

def get_pointtype(types):
    '''Display a dialog requesting user select a point type for
    currect point. Returns the one-char code for the point.'''
    dialog = wx.SingleChoiceDialog(None, "Classify point", \
                                   "Classify point", types)
    if dialog.ShowModal() == wx.ID_OK:
        return encode_pointtype(dialog.GetStringSelection())
    else:
        raise UserWarning("Assuming continuation point")
        return 'c'


def get_branch():
    '''Display a text dialog for entry of this branch location'''
    dialog = wx.TextEntryDialog( \
                        None, "Enter branch location", \
                        "Branch location", "-1", \
                        style=wx.OK|wx.CANCEL)
    if dialog.ShowModal() == wx.ID_OK:
        result = dialog.GetValue()
        # need some validation here
    else: result = None
    dialog.Destroy()
    return result


#-----------------------------------------------------------------------------
class ThreadWindow(wx.Frame):
    '''Frame class for ThreadMaster. Handles plotting of image and points
    and initial handling of user click events.'''
    def __init__(self, parent, pxspc, offset):
        self.parent = parent
        wx.Frame.__init__(self, None, -1, "Tracer")

        myPanel = wx.Panel(self, -1)
        self.fig = Figure((3,3), 75)
        self.canvas = FigureCanvasWxAgg(myPanel, -1, self.fig)
        self.toolbar = NavigationToolbar2WxAgg(self.canvas)
        self.toolbar.Realize()
        self.CreateStatusBar()
        self.figmgr = FigureManager(self.canvas, 1, self)

        panelSizer = wx.BoxSizer(wx.VERTICAL)
        panelSizer.Add(self.canvas, 1, wx.LEFT|wx.TOP|wx.GROW)
        panelSizer.Add(self.toolbar, 0, wx.GROW)
        myPanel.SetSizer(panelSizer)
        myPanel.Fit()
        
        btnStart = wx.Button(self, -1, 'Start', size=(100,30))
        btnStart.Bind(wx.EVT_BUTTON, self.btnStart_Click)
        btnEdit = wx.Button(self, -1, 'Edit', size=(100,30))
        btnEdit.Bind(wx.EVT_BUTTON, self.btnEdit_Click)
        btnView3d = wx.Button(self, -1, 'View 3-D', size=(100,30))
        btnView3d.Bind(wx.EVT_BUTTON, self.btnView3d_Click)
        btnDiscard = wx.Button(self, -1, 'Discard', size=(100,30))
        btnDiscard.Bind(wx.EVT_BUTTON, self.btnDiscard_Click)
        btnDone = wx.Button(self, -1, 'Done', size=(100,30))
        btnDone.Bind(wx.EVT_BUTTON, self.btnDone_Click)
        ctrlSizer = wx.BoxSizer(wx.HORIZONTAL)
        ctrlSizer.Add(btnStart, 1, wx.LEFT|wx.GROW)
        ctrlSizer.Add(btnEdit, 1, wx.LEFT|wx.GROW)
        ctrlSizer.Add(btnView3d, 1, wx.RIGHT|wx.GROW)
        ctrlSizer.Add(btnDiscard, 1, wx.RIGHT|wx.GROW)
        ctrlSizer.Add(btnDone, 1, wx.RIGHT|wx.GROW)
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(myPanel,  1, wx.LEFT|wx.TOP|wx.GROW)
        sizer.Add(ctrlSizer, 0, wx.LEFT|wx.TOP|wx.GROW)
        
        self.SetSizer(sizer)
        self.Fit()
        self.canvas.mpl_connect('button_press_event', self.ClickHandler)
        self.canvas.mpl_connect('motion_notify_event', self.SetMouseStatus)

        self.pxspc = pxspc
        self.offset = offset


    def btnStart_Click(self, evt):
        self.parent.__call__()

    def btnDone_Click(self, evt):
        self.parent.Done()

    def btnEdit_Click(self, evt):
        self.parent.edit_catalog()

    def btnView3d_Click(self, evt):
        self.parent.view_3d()

    def btnDiscard_Click(self, evt):
        self.parent.Done(discard=True)

    def show_img(self, img):
        ax = self.fig.add_subplot(111)
        ax.imshow(n.flipud(img.transpose()))
#        self.toolbar.update()
        self.Refresh()
 

    def plot_pt(self, pt, plotstr):
        '''Plot points with given formatting string on the current axis
        keeps current x and y limits'''
        ax = self.fig.gca()
        (xmin, xmax), (ymin, ymax) = ax.get_xlim(), ax.get_ylim()
        line = ax.plot((pt[0] + 0.5,), (pt[1] + 0.5,), plotstr)
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
#        self.toolbar.update()
        self.Refresh()
        return line


    def GetToolBar(self):
        return self.toolbar


    def ClickHandler(self, evt):
        '''Handles the click events from the user, and passes them
        on to ThreadMaster'' click_to_fn routine, if not None'''
        #if not (self.mainframe.toolbar.panmode or \
        #        self.mainframe.toolbar.zoomode):
        if self.toolbar.mode == '':
            if evt.inaxes:
                x = evt.xdata
                y = evt.ydata
                statstr = "Click at pt (%f, %f)." % (x,y)
                self.SetStatusText(statstr)
                if not self.parent.click_to_fn == None:
                    self.parent.click_to_fn( (x,y) )
            else:
                self.SetStatusText("Click was outside axes.")
                self.parent.click_to_fn = None


    def delete_line(self, n):
        '''Deletes the n-th to last line from the current axes.'''
        ax = self.fig.gca()
        if ax.lines:
            del ax.lines[-1 - n]
            self.Refresh()


    def SetMouseStatus(self, evt):
        if evt.xdata <> None and evt.ydata <> None:
            xpos, ypos = (evt.xdata + self.offset[0]) * self.pxspc[0], \
                         (evt.ydata + self.offset[1]) * self.pxspc[1]
            try:
                status_str = "x: %4.2f  y: %4.2f" \
                             % (xpos, ypos)
                self.SetStatusText(status_str)
            except:
                pass


#-----------------------------------------------------------------------------
class ThreadMaster:
    
    def __init__(self, stackname, callback, cbargs=None):
        self.click_to_fn = None
        self.init_cube_l = 2
        self.init_cube_w = 20
        self.cb = callback
        self.cbargs = cbargs
        
        self.as = []
        self.dis = []
        self.codes = []
        
        if os.path.isfile(stackname):
            self.stackname = stackname
        else:
            raise ValueError(stackname + " is not a valid file.")
        self.vol = ncdf.r(stackname)
        volmax = self.vol.max(2)

        self.frame = ThreadWindow(self, self.cbargs[0], self.cbargs[1])
        self.frame.Show(True)
        self.frame.show_img(volmax)

        #self.mv = None
        self.sc = None
        self.mvsrc = None


    def __call__(self, keep_existing=True):
        '''called to start tracing - waits for user to click start point
        then continues with with_start_pt'''
        # get starting point
        self.frame.SetStatusText("Select a starting point.")
        if keep_existing == False:
            del self.as[:]
            del self.dis[:]
        self.click_to_fn = self.with_start_pt


    def with_start_pt(self, pt):
        # get z point
        z = coordhandling.zmax(self.vol, pt)
        self.a = n.array((pt[0], pt[1], z))

        # get second pt to calculate D
        self.frame.SetStatusText("Select a second point along the dendrite.")
        self.click_to_fn = self.with_second


    def with_second(self, pt):
        # get z point
        z = coordhandling.zmax(self.vol, pt)
        second = n.array((pt[0], pt[1], z))
        self.D = (second - self.a) / norm(second - self.a)
        self.click_to_fn = None

        # calculate diam
        self.frame.SetStatusText("Calculating diameter")

        lc, uc = get_cnrs(self.a, self.D, self.init_cube_l, \
                          self.init_cube_w, (0,0,0), self.vol.shape)
        bk = self.vol[lc[0]:uc[0], lc[1]:uc[1], lc[2]:uc[2]]
        self.di, self.hpk = fwhm_diam( bk, self.a - lc, \
                                              perpz(self.D) )
        # sanity check
        if self.di == None or self.di > 200:
            self.frame.SetStatusText("Cannot estimate diameter. " \
                                     "Try another point.")
            self.click_to_fn = self.with_start_pt
            return
        # save first point
        self.as.append(self.a)
        self.dis.append(self.di)
        pointtypes = ['soma', 'continuation', 'branch', \
                      'join', 'terminus']
        if self.codes:
            if self.codes[-1] == 'c':
                # implies we're continuing from a break
                # => assume continuation
                code = 'c'
            else:
                code = get_pointtype(pointtypes)
                #branch = get_branch()
        else:
            code = get_pointtype(pointtypes)
            #branch = get_branch()
        self.codes.append(code)
        self.frame.plot_pt(self.a, 'go')
        self.frame.SetStatusText("Idle")

        # begin loop:
        iternum = 0
        while True:
            pdbg('debug verbose', "**", iternum, '*')
            self.a, self.D, self.di, self.hpk = \
                    trace_one(self.vol, self.a, self.D, \
                              self.di, self.hpk, frame=self.frame)
            if self.a is None:
                pdbg('debug terse', "Finished")
                break
            if self.is_stuck():
                pdbg('debug terse', "Stuck")
                break
            if n.allclose(self.a, self.as[-1], atol=1e-1):
                pdbg('debug terse', "Repeating")
                break
            pdbg('debug terse', 'Next point: %0.2f %0.2f %0.2f' % \
                 (self.a[0], self.a[1], self.a[2]))

            # handle valid points here
            self.as.append(self.a)
            self.dis.append(self.di)
            self.codes.append('c')
            self.frame.plot_pt(self.a, 'yo')
            iternum += 1
            #ipbreak()

        # handle end of thread here
        pointtypes = ['soma', 'continuation', 'branch', \
                      'join', 'terminus']
        #end_loc = get_branch()
        code = get_pointtype(pointtypes)
        self.codes[-1] = code # change last point to selected value
                

    def is_stuck(self, back_pts=5, dist=3.0):
        '''check total path length from n points back is at least
        dist'''
        num_pts = len(self.as)
        if num_pts < back_pts:
            return False
        else:
            return ((self.as[-back_pts] - self.a)**2).sum() < dist**2

    def edit_catalog(self):
        '''manually edit the catalog of points'''
        self.frame.SetStatusText("Click near to a point to catalogue, " + \
                                 "or outside axes to finish.")
        self.click_to_fn = self._edit_catalog


    def _edit_catalog(self, pt):
        '''handle user clicks for editing catalog points'''
        self.click_to_fn = None

        dists = n.sqrt( ((n.array((self.as))[:,0:2] \
                              - n.array((pt)))**2).sum(1) )
        whichpt = n.where(dists == dists.min())[0][0]
        pdbg('debug verbose', "Idx of closest point", whichpt)
        #thept = n.array(self.as)[whichpt].squeeze()
        thept = self.as[whichpt]
        self.frame.plot_pt(thept, 'wp')

        pointtypes = ['soma', 'continuation', 'branch', \
                      'spine', 'filopodium', \
                      'join', 'terminus', 'delete']
        dialog = wx.SingleChoiceDialog(None, "Pick an action or label", \
                                       "Classify point", pointtypes)
        skip_pts = 0
        if dialog.ShowModal() == wx.ID_OK:
            opt = dialog.GetStringSelection()
            if opt == 'delete':
                del self.as[whichpt]
                del self.dis[whichpt]
                del self.codes[whichpt]
                self.frame.plot_pt(thept, 'kx')
                skip_pts += 1
            elif opt == 'continuation':
                self.codes[whichpt] = 'c'
            elif opt == 'spine':
                self.codes[whichpt] = 'p'
                self.frame.plot_pt(thept, 'co')
                skip_pts += 1
            elif opt == 'filopodium':
                self.frame.plot_pt(thept, 'cs')
                self.codes[whichpt] = 'f'
                skip_pts += 1
            elif opt == 'branch':
                self.frame.plot_pt(thept, 'mo')
                self.codes[whichpt] = 'b'
                skip_pts += 1
            elif opt == 'terminus':
                self.frame.plot_pt(thept, 'ro')
                self.codes[whichpt] = 't'
                skip_pts += 1
            elif opt == 'join':
                self.codes[whichpt] = 'j'
        self.frame.delete_line(skip_pts)
        dialog.Destroy()
        self.frame.SetStatusText("Click near to a point to catalogue," + \
                                 "or outside axes to finish.")
        self.click_to_fn = self._edit_catalog


    def view_3d(self):
        '''view a sub-volume around a selected point'''
        self.frame.SetStatusText("Click near to a point for the volume" \
                                 " centre, or outside axes to finish.")
        self.click_to_fn = self._view_3d

        
    def _select_pt(self, pt):
        dists = n.sqrt( ((n.array((self.as))[:,0:2] \
                          - n.array((pt)))**2).sum(1) )
        whichpt = n.where(dists == dists.min())[0][0]
        pdbg('debug verbose', "Idx of closest point", whichpt)
        thept = self.as[whichpt]
        return thept


    def _view_3d(self, pt):
        #self.click_to_fn = None
        thept = self._select_pt(pt)
        vol_size = (60,60,60)
        bk = coordhandling.extract_box(self.vol, thept, vol_size)
        from enthought.mayavi.sources.array_source import ArraySource
        from enthought.mayavi.modules.iso_surface import IsoSurface
        from enthought.mayavi.modules.outline import Outline
        from enthought.mayavi.tools import mlab
        if self.sc == None: # was self.mv
            self.sc = mlab.figure() # was self.mv
            #self.sc = self.mv.engine.current_scene
            self.sc.scene.background = (1.,1.,1.,)
        if self.mvsrc == None:
            self.mvsrc = ArraySource()
            self.sc.add_child(self.mvsrc) # was self.sc
            o = Outline()
            isos = IsoSurface()
            self.mvsrc.add_module(isos) # was self.mv
            self.mvsrc.add_module(o) # was self.mv
            o.actor.property.color = (0.,0.,0.)
        self.mvsrc.scalar_data = bk
        
                                             
        
    def Done(self, discard=False):
        thread = None
        if not len(self.as) == 0:
            thread = n.hstack((
                n.array(self.as, dtype=n.object), \
                n.array(self.dis,dtype=n.object), \
                n.array(self.codes, \
                        dtype=n.object)[...,n.newaxis]))
        if (not self.cb is None) and (not discard) and (not thread == None):
            self.cb(self.stackname, thread, *self.cbargs)
        self.frame.Destroy()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


def add_offset(in_file, out_file, offsets_file, resolution):
    offsets = calculate_global_offsets(offsets_file)
    stack_code = in_file.split('/')[-1].split('_')[0]
    try:
        offset = n.array(offsets[stack_code])
    except KeyError:
        print "Offset not found - defaulting to 0"
        offset = n.array([0,0,0])
    print "offset:", offset

    fin = open(in_file, 'r')
    fout = open(out_file, 'w')
    lines = fin.readlines()
    for line in lines:
        if line.strip() == '' or line[0] == '#':
            fout.write(line)
        else:
            bits = line.split()
            xyz = n.array(bits[2:5]).astype(float)
            xyz += (resolution * offset)
            xyzstr = str(xyz).strip('[]')
            strlist = bits[0:2] + [xyzstr] + bits[5:]
            lineout = ' '.join(strlist)
            fout.write(lineout + '\n')
    fin.close()
    fout.close()


def convert_diams(in_file, out_file, pxspc):
    fin = open(in_file, 'r')
    fout = open(out_file, 'w')
    lines = fin.readlines()
    for line in lines:
        if line.strip() == '' or line[0] == '#':
            fout.write(line)
        else:
            bits = line.split()
            diamstr = str(float(bits[5]) * pxspc)
            strlist = bits[0:5] + [diamstr] + [bits[6]]
            lineout = ' '.join(strlist)
            fout.write(lineout + '\n')
    fin.close()
    fout.close()
    

def unshrink_coords(in_file, percent=15):
    factor = 1/(1 - percent/100.)
    fin = open(in_file, 'r')
    out_file = in_file.split('.')[0] + '_unshrunk.txt'
    fout = open(out_file, 'w')
    lines = fin.readlines()
    for line in lines:
        if line.strip() == '' or line[0] == '#':
            fout.write(line)
        else:
            bits = line.split()
            pos = n.array([float(x) for x in bits[2:5]])
            newpos = pos * factor
            newline = ('%s %s %7.3f %7.3f %7.3f   %s  %s\n' % \
                      (bits[0], bits[1], \
                       newpos[0], newpos[1], newpos[2], \
                       bits[5], bits[6]))            
            fout.write(newline)
    fin.close()
    fout.close()

def write_dxt(in_file):
    fin = open(in_file, 'r')
    out_file = in_file.split('.')[0] + '.dxt'
    fout = open(out_file, 'w')
    fout.write('x y z data\n')
    lines = fin.readlines()
    for line in lines:
        if line.strip() == '' or line[0] == '#':
            pass
        else:
            bits = line.split()
            xyz = n.array(bits[2:5]).astype(float)
            xyzstr = str(xyz).strip('[]')
            strlist = [xyzstr] + [bits[5]]
            lineout = ' '.join(strlist)
            fout.write(lineout + '\n')
    fin.close()
    fout.close()
    

def renumber_file(in_file, out_file):
    fin = open(in_file, 'r')
    fout = open(out_file, 'w')
    lines = fin.readlines()
    i = 0
    for line in lines:
        if line.strip() == '' or line[0] == '#':
            fout.write(line)
        else:
            rest = line[4:]
            lid = '%04d' % i
            newline = lid + rest
            fout.write(newline)
            i += 1
    fin.close()
    fout.close()
