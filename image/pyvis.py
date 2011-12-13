import numpy
import wx
import matplotlib

matplotlib.interactive(False)
matplotlib.use("WXAgg")
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg, \
     NavigationToolbar2WxAgg, FigureManager
from matplotlib.figure import Figure
from matplotlib.axes import Subplot

ID_FILE_EXIT         =   101
ID_WINDOW_VIEWS      =   102
ID_WINDOW_DATA       =   103
ID_TOOLS_MAXPROJ     =   104
ID_TOOLS_DELITEM     =   105
ID_TOOLS_LOSEFOCUS   =   106
ID_CMAP_GRAY         =   107
ID_CMAP_JET          =   108


version = 0.100

class MyNavToolbar2WxAgg(NavigationToolbar2WxAgg):
    def __init__(self, canvas, parent):
        self.parent = parent
        return NavigationToolbar2WxAgg.__init__(self, canvas)
    
    def home(self, *args):
        self.parent.ResetLims()
        NavigationToolbar2WxAgg.home(self, *args)
        
#     '''rendered obsolete by use of toolbar.mode (which I didn''t
#     realise existed).'''
#     def __init__(self, canvas):
#         self.zoomode = False
#         self.panmode = False
#         return NavigationToolbar2WxAgg.__init__(self, canvas)
    
#     def zoom(self, *args):
#         self.zoomode = not self.zoomode
#         if self.panmode:
#             self.panmode = False
#         NavigationToolbar2WxAgg.zoom(self, *args)

#     def pan(self, *args):
#         self.panmode = not self.panmode
#         if self.zoomode:
#             self.zoomode = False
#         NavigationToolbar2WxAgg.pan(self, *args)


# ---------------------------------------------------------------
class ViewsFrame(wx.Frame):
    '''
    frame for Views controls'''


    def __init__(self,parent):
        wx.Frame.__init__(self,None,-1,"Views")
        self.parent = parent
        wx.EVT_CLOSE(self, self.OnCloseViews)

        self.sizerouter = wx.BoxSizer(wx.VERTICAL)
        self.sizerinner = wx.BoxSizer(wx.HORIZONTAL)
        self.sizerinner.Add((10,1))
        self.zslider = \
                         wx.Slider( self, -1, 0, 0, 100, \
                                    wx.DefaultPosition, wx.DefaultSize, \
                                    wx.SL_HORIZONTAL | wx.SL_AUTOTICKS | \
                                    wx.SL_LABELS )
        self.sizerinner.Add( wx.StaticText( self, -1, 'z'  ), 1, \
                             wx.ALIGN_CENTER )
        self.sizerinner.Add( self.zslider, 9, wx.EXPAND )
        self.sizerouter.Add( self.sizerinner, 1, wx.EXPAND )
        self.Bind(wx.EVT_SLIDER, self.SliderUpdate)
        self.SetSizer(self.sizerouter)
        self.SetAutoLayout(1)


    def OnCloseViews(self,evt):
        if not evt.CanVeto():
            self.Destroy()
        else:
            evt.Veto()
            self.parent.mainframe.windowmenu.Check(ID_WINDOW_VIEWS, False)
            self.Show(False)


    def SliderUpdate(self,evt):
        if self.zslider.GetValue() != self.parent.zpos:
            self.parent.zpos = self.zslider.GetValue()
            self.parent.UpdateImage(keepaxes=True)


# ---------------------------------------------------------------
class DataFrame(wx.Frame):
    '''
    frame for Data list'''


    def __init__(self,parent):
        wx.Frame.__init__(self,None,-1,"Data")
        self.parent = parent
        wx.EVT_CLOSE(self, self.OnCloseViews)

        id = wx.NewId()
        self.list = wx.ListCtrl(self,id,style=wx.LC_REPORT | wx.LC_SINGLE_SEL)
        self.list.InsertColumn(0, "name")
        self.list.InsertColumn(1, "shape")
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.list,1,wx.EXPAND)
        self.Bind(wx.EVT_LIST_ITEM_SELECTED, self.OnSelectItem)
        self.SetSizer(self.sizer)
        self.SetAutoLayout(1)


    def OnCloseViews(self,evt):
        if not evt.CanVeto():
            self.Destroy()
        else:
            evt.Veto()
            self.parent.mainframe.windowmenu.Check(ID_WINDOW_DATA, False)
            self.Show(False)


    def ShapeString(self, id):
        ndims = numpy.rank( self.parent.ims[id] )
        str = ''
        for i in range( ndims ):
            dimsize = self.parent.ims[id].shape[i]
            if dimsize > 1:
                if i > 0:
                    str = str + ' x '
                str = str + "%d" % dimsize
        return str


    def AddItem(self, name):
        id = self.list.GetItemCount()
        self.list.InsertStringItem(id, name)
        self.list.SetStringItem(id, 1, self.ShapeString(name) )
    

    def DelItem(self, name):
        id = self.list.FindItem(-1, name)
        self.list.DeleteItem(id)
        # select last remaining item in list
        self.list.Select(self.list.GetItemCount() - 1)
        

    def OnSelectItem(self, evt):
        self.parent._pdbg(3, '->dataframe.OnSelectItem')
        if self.parent.curname != evt.GetText():
            self.parent.SwitchToData(evt.GetText())


    def SelectItem(self, name):
        self.parent._pdbg(3, '->dataframe.SelectItem')
        id = self.list.FindItem(-1, name)
        self.list.Select(id)


    def GetCurrentItemText(self):
        sid = self.list.GetFirstSelected()
        if sid > 0:
            return self.list.GetItemText(sid)
        else:
            return None


# ---------------------------------------------------------------
class PlotFigure(wx.Frame):

    def __init__(self, parent):
        wx.Frame.__init__(self, None, -1, "PyVis")

        self.parent = parent
        self.fig = Figure((5,4), 75)

        # canvas stuff
        self.canvas = FigureCanvasWxAgg(self, -1, self.fig)
        self.canvas.mpl_connect('motion_notify_event', self.SetMouseStatus)

        
        self.toolbar = MyNavToolbar2WxAgg(self.canvas, self)
        #self.toolbar = NavigationToolbar2WxAgg(self.canvas)
        self.toolbar.Realize()
        self.CreateStatusBar()
        self.cmap = None

        # file menu
        filemenu = wx.Menu()
        filemenu.Append(ID_FILE_EXIT, "E&xit", "Terminate the program")
        self.Connect(ID_FILE_EXIT, -1, wx.wxEVT_COMMAND_MENU_SELECTED, \
                     self.OnMenuExit)

        # window menu
        self.windowmenu = wx.Menu()
        self.windowmenu.AppendCheckItem(ID_WINDOW_DATA, "&Data",
                                        "Open/close Data window")
        self.windowmenu.AppendCheckItem(ID_WINDOW_VIEWS, "&Views",
                                        "Toggle open state of Views window")
        self.Connect(ID_WINDOW_VIEWS, -1, wx.wxEVT_COMMAND_MENU_SELECTED, \
                     self.ShowHideViews)

# tools menu
        toolsmenu = wx.Menu()
        toolsmenu.Append(ID_TOOLS_MAXPROJ, "&Maximum projection", \
                        "Show maximum projection over selected axes" )
        toolsmenu.Append(ID_TOOLS_DELITEM, "&Delete", \
                         "Delete the currently selected item" )
        self.Connect(ID_TOOLS_MAXPROJ, -1,  wx.wxEVT_COMMAND_MENU_SELECTED, \
                     self.ShowMaxProj)
        self.Connect(ID_TOOLS_DELITEM, -1, wx.wxEVT_COMMAND_MENU_SELECTED, \
                     self.DelCurItem)

        # colourmap menu
        self.cmapmenu = wx.Menu()
        self.cmapmenu.AppendCheckItem(ID_CMAP_GRAY, "&Gray", \
                        "Use Gray colour map")
        self.cmapmenu.AppendCheckItem(ID_CMAP_JET, "&Jet", \
                        "Use Jet colour map")
        self.Connect(ID_CMAP_GRAY, -1, wx.wxEVT_COMMAND_MENU_SELECTED, \
                    self.CMapGray)
        self.Connect(ID_CMAP_JET, -1, wx.wxEVT_COMMAND_MENU_SELECTED, \
                    self.CMapJet)
        
        
        # menu bar
        menuBar = wx.MenuBar()
        menuBar.Append(filemenu, "&File")
        menuBar.Append(self.windowmenu, "&Window")
        menuBar.Append(toolsmenu, "&Tools")
        menuBar.Append(self.cmapmenu, "&Colourmaps")
        self.SetMenuBar(menuBar)

        self.Bind(wx.EVT_CLOSE, self.OnCloseWindow)
        wx.EVT_MENU(self, ID_WINDOW_DATA, self.ShowHideData)

        self.figmgr = FigureManager(self.canvas, 1, self)
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.canvas, 1, wx.LEFT|wx.TOP|wx.GROW)
        sizer.Add(self.toolbar, 0, wx.GROW)
        self.SetSizer(sizer)
        self.Fit()


    def plot_data(self, array, vmin=None, vmax=None, keepaxes=False):
        a = self.fig.add_subplot(111)
        if keepaxes:
            xlims = a.get_xlim()
            ylims = a.get_ylim()
        a.imshow(array, interpolation='nearest', \
                 vmin=vmin, vmax=vmax, cmap=self.cmap)
        if keepaxes:
            a.set_xlim(xlims)
            a.set_ylim(ylims)
        self.toolbar.update()
        self.Refresh()


    def GetLims(self):
        a = self.fig.add_subplot(111)
        xlims = a.get_xlim()
        ylims = a.get_ylim()
        return xlims, ylims


    def ResetLims(self):
        a = self.fig.add_subplot(111)
        maxes = self.parent.ims[self.parent.curname].shape
        a.set_xlim((0,maxes[0]))
        a.set_ylim((0,maxes[1]))
        self.Refresh()

    def GetToolbar(self):
        return self.toolbar


    def OnCloseWindow(self, evt):
        if not self.IsBeingDeleted():
            self.parent.viewsframe.Close(True)
            self.parent.dataframe.Close(True)
            self.Destroy()


    def OnMenuExit(self, evt):
        self.parent.viewsframe.Close(True)
        self.parent.dataframe.Close(True)
        self.Destroy()


    def ShowHideData(self, evt):
        self.parent.dataframe.Show(not self.parent.dataframe.IsShown())

       
    def ShowHideViews(self, evt):
        self.parent.viewsframe.Show(not self.parent.viewsframe.IsShown())


    def SetMouseStatus(self, evt):
        xpos, ypos = evt.xdata, evt.ydata
        if xpos != None and ypos != None:
            xpos = numpy.round(xpos)
            ypos = numpy.round(ypos)
            zpos = self.parent.zpos
            # should probably do proper check that exact upper value isn't met
            try:
                val = self.parent.im[xpos, ypos]
                status_str = "x: %d  y: %d  z: %d  I:%4.2f" \
                             % (xpos, ypos, zpos, val)
                self.SetStatusText(status_str)
            except:
                pass


    def ShowMaxProj(self, evt):
        self.parent.MaxProj(self.parent.curname)

        
    def DelCurItem(self, evt):
        self.parent.DelItem( self.parent.curname )


    def CMapGray(self, evt):
        self.cmap = matplotlib.cm.gray
        self.cmapmenu.Check(ID_CMAP_GRAY, True)
        self.cmapmenu.Check(ID_CMAP_JET, False)
        self.parent.UpdateImage(keepaxes=True)


    def CMapJet(self, evt):
        self.cmap = matplotlib.cm.jet
        self.cmapmenu.Check(ID_CMAP_GRAY, False)
        self.cmapmenu.Check(ID_CMAP_JET, True)
        self.parent.UpdateImage(keepaxes=True)
        
        
# ---------------------------------------------------------------
class pyvis(wx.App):

    verbosity = 0

    def OnInit(self):
        '''Usage pv = pyvis.pyvis()
        pv.AddImage(arr, 'arname')

        PyVis array visualization tool, based loosely after the
        excellent KVis by Richard Gooch, which connects to pdl'''

        self._pdbg(3, '->init')
        
        self.nextnum = 0
        self.zpos = 0
        self.ims = {}
        self.mainframe = PlotFigure(self)
        self.mainframe.Show(True)
        self.dataframe = DataFrame(self)
        self.viewsframe = ViewsFrame(self)
        self.SetTopWindow(self.mainframe)
        self.cid = []
        self.curname = None
        return True


    def _pdbg(self, level, str):
        if self.verbosity > level:
            print str


    def AddImg(self, stack, name=None):
        '''Usage:  AddImg(data, name=name)
        Add an image (2 or 3-D) to PyVis stack, with an optional
        name for reference.'''

        # check we're dealing with a numpy array
        if type(stack) == numpy.ndarray:

            self._pdbg(3, '->AddImg')
            # set appropriate data name
            if name == None:
                name = "image %d" % self.nextnum
                self.nextnum += 1
            if name in self.ims.keys():
                # ensure unique image name
                name = name + '%d' % self.nextnum
                self.nextnum += 1

            # load data into memory
            self.ims[name] = stack
            self.dataframe.AddItem(name)
            self.SwitchToData(name)

        else:
            #probably should do better than silently ignore
            pass
        

    def UpdateImage(self, keepaxes=False):
        '''Internal usage mainly. Used to redraw the current image
        (2-D). Keepaxes is used when switching z-position to maintain
        the axis position (zoom level).'''

        self._pdbg(3, '->UpdateImage')

        # for coding convenience
        data = self.ims[self.curname]

        # load data (plane if 3-D) into image
        if numpy.rank(data) == 2:
            self.im = data
        elif numpy.rank(data) == 3:
            self.im = data[:,:,self.zpos]
        else:
            raise TypeError("Data has too many dimensions")

        # display image
        self.mainframe.plot_data(self.im.transpose(), \
                                 vmin=self.imin, vmax=self.imax, \
                                 keepaxes=keepaxes)


    def SwitchToData(self, name):
        '''Usage:  SwitchToData(name)
        
        Switch to viewing the data named <name>.'''

        self._pdbg(3, '->SwitchToData')

        if name == None:
            self.mainframe.fig.clf()
            self.curname == None
        else:
            
            # working out if should keep same x-y position/zoom
            # as long as some old image is present...
            if not self.curname is None:
                old_lims = self.mainframe.GetLims()
                self._pdbg( 3, "in SwitchToData: old lims" + str( old_lims ) )
            else:
                self._pdbg(3, "in SwitchToData: no old lims")
                old_lims = None

            # for coding convenience
            data = self.ims[name]

            # for colour-table preservation throughout stack
            self.imin = data.min()
            self.imax = data.max()

            # if 2-D set maximum z to be 0
            if numpy.rank(data) == 2:
                zmax = 0
                # if 3-D go to same (current) plane or max available
            elif numpy.rank(data) == 3:
                if self.zpos > data.shape[2]:
                    self.zpos = data.shape[2] - 1
                zmax = data.shape[2] - 1
        
            # for later reference
            self.curname = name

            # make sure current item selected in Data window is correct
            if self.dataframe.GetCurrentItemText() != name:
                self.dataframe.SelectItem(name)

            # update maximum for z-slider in Views window
            self.viewsframe.zslider.SetMax( zmax )

            # continue working out if old x-y position
            # is appropriate for new data
            keepaxes = False
            if not old_lims is None:
                if numpy.all( numpy.array(old_lims).max(1) < data.shape[0:2] ):
                    keepaxes = True

            # update the screen image
            self.UpdateImage(keepaxes)


    def DelItem(self, name):
        '''Usage:  DelItem(name)
        
        Deletes item <name> from stack (and switches to last item if
        deleted item was the current one).'''
        if self.curname == name:
            lastname = self.ims.keys()[-1]
            if lastname != name:
                # if we're not deleting last (most recent) item
                # - switch to that one
                nextname = lastname
            else:
                # we are deleting last item
                if len( self.ims ) > 1:
                    # there is another one to switch to
                    nextname = self.ims.keys()[-2]
                else:
                    # we're deleting the only item
                    nextname = None
            self.SwitchToData( nextname  )
        self.dataframe.DelItem(name)
        del self.ims[name]
        

    def MaxProj(self, name):
        '''Usage:  MaxProj(name)

        Calculates a maximum z-projection of the data <name>
        and switches to make it the current data.'''
        maxprojx = self.ims[name].max(2)
        self.AddImg(maxprojx, name + ' (max projection)')


    def ClickHandler(self, evt):
        #if not (self.mainframe.toolbar.panmode or \
        #        self.mainframe.toolbar.zoomode):
        if self.mainframe.toolbar.mode == '':
            x = numpy.round(evt.xdata)
            y = numpy.round(evt.ydata)
            z = self.zpos
            if not self.clickcallback is None:
                #print "To callback."
                self.clickcallback(x,y,z, *self.clickcallbackargs)
        
        
    def ConnectClickCallback(self, func, *args):
        '''Usage:  ConnectClickCallback(func, *args)

        Connects function <func> to respond to non-zoom, non-pan click
        events in the current canvas. func must have the signature:

        def <callbackfunc>(x,y,z, *args):
            ...

        where x,y,z are the click co-ordinates when inside the
        axes and None (x and y anyway) when click is outside and
        args are any other parameters the callback needs.'''
        print "Connecting", func.__name__
        self.clickcallback = func
        self.clickcallbackargs = args
        self.cid.append( self.mainframe.canvas.mpl_connect(\
            'button_press_event', self.ClickHandler) )


    def DisconnectClickCallback(self):
        '''Usage: DisconnectClickCallback()

        Disconnects all connected callback routines from click event.'''
        while len(self.cid) > 0:
            i = self.cid.pop()
            self.mainframe.canvas.mpl_disconnect(i)
                
