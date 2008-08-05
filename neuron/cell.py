import numpy as n
import os.path
import wx
from neuron.thread import ThreadMaster
from neuron import pdbg
import sqlite
from copy import copy

debug = 'debug terse'

ID_OPEN_DB      = 101
ID_NEW_DB       = 102
ID_ADD_STACK    = 103
ID_EDIT_STACK   = 104
ID_DELETE_STACK = 105
ID_SAVE_DB      = 106
ID_SET_BASEDIR  = 107
ID_TREE_WINDOW  = 108

#-----------------------------------------------------------------------------
# notes:
#
# decided to independently assign parents (even though I could potentially
# calculate it from the offsets) because in the long run it would be good to
# automate calculation of the offsets, via correlation alignment, from the
# parent information
#
#-----------------------------------------------------------------------------

def _create_db_file(filename):
    '''Creates a new database file. Presumes that no file of that
    name exists currently.'''
    if os.path.isfile(filename):
        print "File already exists!"
        # maybe should delete it here - since has checked in selection dialog
        # definitely shouldn't handle like this...
        raise ValueError("Existing filename given")
    else:
        conn = sqlite.connect(filename, isolation_level=None)
        cur = conn.cursor()
        cur.execute('create table points (idx INTEGER,' \
                    'parent INTEGER, code TEXT, branch TEXT' \
                    'ax REAL, ay REAL, az REAL,' \
                    'di REAL, cdist REAL)')
        cur.execute('create table stacks (idx INTEGER, filename TEXT, ' \
                    'xoffs INTEGER, yoffs INTEGER, zoffs INTEGER, ' \
                    'xspc REAL, yspc REAL, zspc REAL, parent INTEGER)')
        cur.execute('create table basedir (filename TEXT)')
        conn.commit()
        cur.close()
        return conn


def _open_db_file(filename):
    '''Opens an existing database file.
    Presumes file''s existence has already been verified.'''
    if os.path.isfile(filename):
        conn = sqlite.connect(filename, isolation_level=None)
        return conn
    else:
        raise ValueError("An invalid filename supplied")
        

def _read_stack_entries(conn):
    '''Reads stack entries from database file and returns as a
    list of lists'''
    cur = conn.cursor()
    cur.execute('select * from stacks')
    data = []
    for row in cur:
        data.append(row)
    cur.close()
    return data


def _read_point_entries(conn):
    '''Reads point entries from database file and returns as an object array'''
    cur = conn.cursor()
    cur.execute('select * from points')
    data = []
    for row in cur:
        data.append(row)
    cur.close()
    return n.array(data,dtype=object)


def _write_stack_entries(conn, data):
    '''Writes an object array with rows of:

    (filename, offset-x,y,z, spacing-x,y,z)

    into database'''
    cur = conn.cursor()
    cur.execute("DELETE FROM stacks")
    for i,t in enumerate(list(data)):
        cur.execute("INSERT INTO stacks VALUES (%s,%s,%s,%s,"
                    "%s,%s,%s,%s,%s)", [i]+list(t))
    conn.commit()
    cur.close()


def _write_point_row(row,cur):
    '''Writes a single point row to the db'''
    t = list(row)
    cur.execute("INSERT INTO points VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s)", \
                t)

def _write_point_entries(conn, data):
    '''Writes a list of lists (index, parent, code, branch, x,y,z, di)
    into points table'''
    cur = conn.cursor()
    # start with a fresh slate
    cur.execute("DELETE FROM points")
    apply_along_axis(_write_point_row,1,data,cur)
    conn.commit()
    cur.close()
    

def _write_basedir(conn, basedir):
    '''Writes basedir into db'''
    cur = conn.cursor()
    cur.execute("DELETE FROM basedir")
    cur.execute("INSERT INTO basedir VALUES (%s)", basedir)
    conn.commit()
    cur.close()


def _read_basedir(conn):
    '''Reads basedir from db'''
    cur = conn.cursor()
    cur.execute("SELECT * FROM basedir")
    i = -1
    for i,f in enumerate(cur):
        basedir = f
    if i > 0:
        raise ValueError("There should only be one basedir entry in db!")
    elif i == -1:
        return ""
    else:
        return basedir[0]
        

def _get_closest_pt(loc, locs):
    '''Return index of closest point in locs to loc'''
    diststo = ((locs - loc)**2).sum(1)
    return n.where(diststo == diststo.min())


#-----------------------------------------------------------------------------
# class Branch:
#     '''Abstract class describing one dendritic branch from any possible
#     level of the dendritic tree'''
#     def __init__(self, parent, points_ids):
#         self.name = name
#         self.parent = parent
#         self.points = points_ids
#         self.order = parent.order + 1
#         self.finished = None
#-----------------------------------------------------------------------------

class StackWindow(wx.Frame):
    '''Frame class for StackHandler'''
    def __init__(self, parent):
        self.parent = parent
        wx.Frame.__init__(self, None, -1, "Stacks", size=(700,300))
        self.CreateStatusBar()

        # buttons
        btnSizer = wx.BoxSizer(wx.HORIZONTAL)
        btnAddStack = wx.Button(self, -1, 'Add Stack')
        btnAddStack.Bind(wx.EVT_BUTTON, self.btnAddStack_Click)
        btnSizer.Add(btnAddStack, 1, wx.ALIGN_CENTRE_VERTICAL)
        btnEditStack = wx.Button(self, -1, 'Edit Stack')
        btnEditStack.Bind(wx.EVT_BUTTON, self.btnEditStack_Click)
        btnSizer.Add(btnEditStack, 1, wx.ALIGN_CENTRE_VERTICAL)
        btnTraceStack = wx.Button(self, -1, 'Trace Stack')
        btnTraceStack.Bind(wx.EVT_BUTTON, self.btnTraceStack_Click)
        btnSizer.Add(btnTraceStack, 1, wx.ALIGN_CENTRE_VERTICAL)

        # base directory
#        self.txtBaseDir = wx.StaticText(self, -1, "Base Directory: ", \
#                                        size=(-1,20))

        # stacklist
        self.stacklist = wx.ListView(self, -1, \
                                     style=wx.LC_REPORT | wx.LC_SINGLE_SEL, \
                                     size=(400,100))
        self.stacklist.InsertColumn(0, 'index')
        self.stacklist.SetColumnWidth(0, 50)
        self.stacklist.InsertColumn(1, 'file name')
        self.stacklist.SetColumnWidth(1, 300)
        self.stacklist.InsertColumn(2, 'offset')
        self.stacklist.SetColumnWidth(2, 150)
        self.stacklist.InsertColumn(3, 'spacing')
        self.stacklist.SetColumnWidth(3, 150)
        self.stacklist.InsertColumn(4, 'parent')
        self.stacklist.SetColumnWidth(4, 50)

        # window sizer
        sizer = wx.BoxSizer(wx.VERTICAL)
#        sizer.Add(self.txtBaseDir, 0, wx.EXPAND)
        sizer.Add(self.stacklist, 2, wx.EXPAND)
        sizer.Add(btnSizer, 0, wx.EXPAND)
        self.SetSizer(sizer)

        # menus
        # database menu
        DBMenu = wx.Menu()
        DBMenu.Append(ID_NEW_DB, "&New", "Create a new database file.")
        self.Connect(ID_NEW_DB, -1, wx.wxEVT_COMMAND_MENU_SELECTED, \
                     self.OnNewDB)
        DBMenu.Append(ID_OPEN_DB, "&Open", "Open an existing database file.")
        self.Connect(ID_OPEN_DB, -1, wx.wxEVT_COMMAND_MENU_SELECTED, \
                     self.OnOpenDB)
        DBMenu.Append(ID_SAVE_DB, "&Save", "Save stacks to open database.")
        self.Connect(ID_SAVE_DB, -1, wx.wxEVT_COMMAND_MENU_SELECTED, \
                     self.OnSave)

        # stack menu
        StackMenu = wx.Menu()
        StackMenu.Append(ID_ADD_STACK, "&Add", "Add a stack.")
        self.Connect(ID_ADD_STACK, -1, wx.wxEVT_COMMAND_MENU_SELECTED, \
                     self.btnAddStack_Click)
        StackMenu.Append(ID_EDIT_STACK, "&Edit", "Edit a stack.")
        self.Connect(ID_EDIT_STACK, -1, wx.wxEVT_COMMAND_MENU_SELECTED, \
                     self.btnEditStack_Click)
        StackMenu.Append(ID_DELETE_STACK, "&Delete", "Delete a stack")
        self.Connect(ID_DELETE_STACK, -1, wx.wxEVT_COMMAND_MENU_SELECTED, \
                     self.DeleteStack)
        StackMenu.AppendSeparator()
        StackMenu.Append(ID_SET_BASEDIR, "&Base Directory", \
                         "Set base directory")
        self.Connect(ID_SET_BASEDIR, -1, wx.wxEVT_COMMAND_MENU_SELECTED, \
                     self.SetBaseDir)

        # window menu
        WindowMenu = wx.Menu()
        WindowMenu.Append(ID_TREE_WINDOW, "&Tree Viewer", \
                          "Toggle Visibility of Tree Viewer")
        self.Connect(ID_TREE_WINDOW, -1, wx.wxEVT_COMMAND_MENU_SELECTED, \
                     self.ToggleTreeWindow)

        # menu bar
        menuBar = wx.MenuBar()
        menuBar.Append(DBMenu, "&Database")
        menuBar.Append(StackMenu, "&Stack")
        menuBar.Append(WindowMenu, "&Window")
        self.SetMenuBar(menuBar)
        

    def ToggleTreeWindow(self, evt):
        if self.parent.treeFrame.IsShown():
            self.parent.treeFrame.Show(False)
        else:
            self.parent.treeFrame.Show(True)


    def btnAddStack_Click(self, evt):
        self.parent.AddStack()

    
    def btnEditStack_Click(self, evt):
        print "Edit stack clicked"


    def DeleteStack(self, evt):
        print "Delete stack clicked"


    def btnTraceStack_Click(self, evt):
        if self.stacklist.GetSelectedItemCount() != 1:
            self.SetStatusText("Exactly one item must be selected.")
        else:
            selectedItem = self.stacklist.GetFirstSelected()
            item = self.stacklist.GetItem(selectedItem, 1)
            fname = item.m_text
            self.parent.tm = ThreadMaster(self.parent.basedir + fname, \
                              self.parent.ProcessThread)
            print "Trace stack clicked, with %s selected." \
                  % (self.parent.basedir + fname)


    def OnOpenDB(self, evt):
        self.parent.OpenDB()
        

    def OnNewDB(self, evt):
        self.parent.NewDB()


    def OnSave(self, evt):
        self.parent.SaveToDB()


    def ExtractData(self):
        data = []
        n_items = self.stacklist.GetItemCount()
        for i in xrange(n_items):
            item = self.stacklist.GetItem(i,1)
            fname = item.m_text
            item = self.stacklist.GetItem(i,2)
            offsets = item.m_text
            item = self.stacklist.GetItem(i,3)
            spacings = item.m_text
            item = self.stacklist.GetItem(i,4)
            parent = item.m_text
            data.append([fname] + offsets.split(" ") \
                        + spacings.split(" ") + [parent])
        return data


    def SetBaseDir(self, evt):
        '''Sets the base directory for stack paths'''
        dirDialog = wx.DirDialog(self, "Choose a base directory:")
        if dirDialog.ShowModal() == wx.ID_OK:
            self.parent.basedir = dirDialog.GetPath() + '/'
            print "New base directory: %s" % self.parent.basedir
        dirDialog.Destroy()
        n_items = self.stacklist.GetItemCount()
        for i in xrange(n_items):
            item = self.stacklist.GetItem(i,1)
            fname = item.m_text
            if fname.startswith(self.parent.basedir):
                self.stacklist.SetStringItem(
                    i,1,fname[len(self.parent.basedir):])
        self.SetStatusText("Basedir is %s" % self.parent.basedir)
        # parse stacklist and remove basedir where present
        


    def AddRow(self, row):
        '''Adds a row to stack list'''
        nid = self.parent.nextStack
        self.stacklist.InsertStringItem(nid, str(nid))
        self.stacklist.SetStringItem(nid, 1, row[1])
        offset_str = " ".join( \
            str( row[2:5] ).strip("()").split(", ") )
        self.stacklist.SetStringItem(nid, 2, offset_str)
        spacing_str = " ".join( \
            str( row[5:8] ).strip("()").split(", ") )
        self.stacklist.SetStringItem(nid, 3, spacing_str)
        self.stacklist.SetStringItem(nid, 4, str(row[8]))


    def GetRowByName(self, name):
        '''Returns the stack row containing the given name'''
        n_items = self.stacklist.GetItemCount()
        i = 0
        done = False
        while i > n_items and not done:
            item = self.stacklist.GetItem(i,1)
            if item.m_text == name:
                done = True
            else:
                i += 1

        item = self.stacklist.GetItem(i,0)
        idx = int(item.m_text)
        
        item = self.stacklist.GetItem(i,1)
        name = item.m_text

        item = self.stacklist.GetItem(i,2)
        offset = [int(x) for x in item.m_text.split(" ")]

        item = self.stacklist.GetItem(i,3)
        spacing = [float(x) for x in item.m_text.split(" ")]

        item = self.stacklist.GetItem(i,4)
        parent = int(item.m_text)
        return [idx, name] + offset + spacing + [parent]
        
        

#-----------------------------------------------------------------------------
class TreeViewer(wx.Frame):
    '''Viewer window showing tree outline view of processed points'''
    def __init__(self, parent):
        self.parent = parent
        wx.Frame.__init__(self, None, -1, "Tree Viewer", size=(400,400))

        self.tree = wx.TreeCtrl(self, -1, \
                           style=wx.TR_HAS_BUTTONS| \
                                wx.SUNKEN_BORDER| \
                                wx.TR_FULL_ROW_HIGHLIGHT)
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.tree, 1, wx.ALL|wx.GROW, border=10)
        self.SetSizer(sizer)
        self.root = self.tree.AddRoot('soma')
        

#-----------------------------------------------------------------------------
class StackHandler(wx.App):
    def OnInit(self):
        # class variables
        self.stackFrame = StackWindow(self)
        self.stackFrame.Show(True)

        self.treeFrame = TreeViewer(self)
        #self.treeFrame.Show(True)

        self.nextPoint = -1
        self.nextStack = 0

        opts = ['New', 'Existing']
        dialog = wx.SingleChoiceDialog(None, "Select file type to open:", \
                                    "Open...", opts)
        if dialog.ShowModal() == wx.ID_OK:
            dialog.Destroy()

            # create a new db
            if dialog.GetStringSelection() == 'New':
                self.NewDB()
            # open an existing db
            elif dialog.GetStringSelection() == 'Existing':
                self.OpenDB()
        else:
            self.stackFrame.SetStatusText("No database open.")
            # user cancelled
            dialog.Destroy()

        return True
    

    def OpenDB(self):
        file_dialog = wx.FileDialog( \
            None, "Pick an existing file:", \
            defaultDir=os.getcwd(), defaultFile="cell.db", \
            wildcard="*.db",  style=wx.OPEN, \
            pos=wx.DefaultPosition)
        if file_dialog.ShowModal() == wx.ID_OK:
            fname = file_dialog.GetPath()
            self.stackFrame.SetStatusText("Connected to " + fname)
            #print "Selected file name: " + fname
            self.dbcon = _open_db_file(fname)

            # read stack
            all_stacks = _read_stack_entries(self.dbcon)
            for stack in all_stacks:
                self.stackFrame.AddRow(stack)
            self.basedir = _read_basedir(self.dbcon)

            # read points
            # points are stored as idx, parent, code, x,y,z, di, cdist
            all_points = n.array(_read_point_entries(self.dbcon), \
                                        dtype=object)
            if len(all_points.shape) > 1:
                self.nextPoint = all_points[-1][0] + 1
            else:
                self.nextPoint = 0
            pdbg('debug terse', 'Next point set to %d' % self.nextPoint)
            return True
            
        file_dialog.Destroy()
        return False


    def NewDB(self):
        file_dialog = wx.FileDialog( \
            None, "Pick a file name:", \
            defaultDir=os.getcwd(), \
            wildcard="*.db", \
            style = \
            wx.SAVE | wx.OVERWRITE_PROMPT, \
            pos=wx.DefaultPosition)
        if file_dialog.ShowModal() == wx.ID_OK:
            fname = file_dialog.GetPath()
            self.stackFrame.SetStatusText("Connected to " + fname)
            #print "Selected file name: " + fname
            self.dbcon = _create_db_file(fname)
            self.all_stacks = None
            self.basedir = ""
            self.nextPoint = 0
            return True

        file_dialog.Destroy()
        return False


    def SaveToDB(self):
        '''Save stacks, basedir and points to database,
        overriding old values'''
        # format data correctly
        data = self.stackFrame.ExtractData()
        _write_stack_entries(self.dbcon, data)
        _write_basedir(self.dbcon, self.basedir)
        #_write_point_entries(self.dbcon, self.all_points)

    def ProcessThread(self, stackName, thread):
        '''Process co-ordinates returned from ThreadMaster'''

        thread = copy(thread)

        n_points = thread.shape[0]
        points = thread[:,0:3].astype(n.float)
        dis = thread[:,3].astype(n.float)
        codes = thread[:,4]

        # Remove basedir from stackName if needed
        if stackName.startswith(self.basedir):
            stackName = stackName[len(self.basedir):]
            
        # convert co-ordinates to global co-ordinates
        stackRow = self.stackFrame.GetRowByName(stackName)
        offset = tuple(stackRow[2:5])
        global_points = points + n.array(offset)

        # convert to microns
        resolution = tuple(stackRow[5:8])
        real_points = global_points * n.array(resolution)
        # (ideally this would calculate local microns and then
        # add a global micron offset - but I'm not sure how
        # accurate this would be, so I'll do it the other (simpler)
        # way first to compare

        # Assign appropriate parent code to first co-ordinate
        if codes[0] == 's':
            print "first point is soma"  
            parent = -1
            
        elif codes[0] == 'j':
            print "first point is join"
            # need to identify parent join

            # there will be only one matching join
            # - location may vary

            # matching code will be same branch number, letter - 1
            this_br_num = codes[1,0]
            this_br_part = codes[1,0]
            # what happens if this_br_part is 'a'?
            if this_br_part == 'a':
                raise ValueError("Sorry, I haven't implemented this option" \
                                 "Try tracing from proximal to distal.")
            prev_br_part = chr(ord(this_br_part) - 1)
            print "Will be looking for %s%d" % (this_br_num, prev_br_part)
            
        elif codes[0] == 'b':
            print "first point is branch"
            # need to identify parent branch
            
            # there may be multiple matching branches
            # - location should be similar
            loc = global_points[0]
            locs = self.all_locs[self.all_codes[0] == 'b']
            closest = _get_closest_pt(loc, self.all_locs)
            print "Closest point is at %s" % (closest)

        # Assign correct id and parent values
        # points will append to end of list
        idxs = []
        parents = []
        idx = n.arange(n_points) + self.nextPoint
        self.nextPoint += n_points
        parents = idx - 1 # needs changing
        
        # construct all_points array
        1/0

        # Save in database'''
        del self.tm


    def AddStack(self):
        # add item with
        file_dialog = wx.FileDialog( \
            None, "Select stack file:", defaultDir=os.getcwd(), \
            wildcard="NetCDF (*.nc)|*.nc|" \
            "TIF stack (*.TIF or *.tif)|*.[Tt][Ii][Ff]|" \
            "Bitmap stack (*.bmp)|*.bmp", \
            style=wx.OPEN, pos=wx.DefaultPosition)
        if file_dialog.ShowModal() == wx.ID_OK:
            fname = file_dialog.GetPath()
            # remove basedir bit of path, if present
            if fname.startswith(self.basedir):
                fname = fname[len(self.basedir):]
            
            offset_dialog = wx.TextEntryDialog( \
                None, "Enter stack's global offset, in the format " \
                "x y z:", "Offset", "0 0 0", style=wx.OK|wx.CANCEL)
            if offset_dialog.ShowModal() == wx.ID_OK:
                offset = offset_dialog.GetValue()
                # need some validation here
                
                space_dialog = wx.TextEntryDialog( \
                    None, "Enter stack's pixel spacing, in the format " \
                    "x y z:", "Pixel Spacing", "0.104 0.104 0.4", \
                    style=wx.OK|wx.CANCEL)
                if space_dialog.ShowModal() == wx.ID_OK:
                    spacing = space_dialog.GetValue()
                    # need some validation here

                    parent_dialog = wx.TextEntryDialog( \
                        None, "Enter stack's parent index", \
                        "Stack parent", "-1", \
                        style=wx.OK|wx.CANCEL)
                    if parent_dialog.ShowModal() == wx.ID_OK:
                        parent = parent_dialog.GetValue()
                        # need some validation here
                        stack_row = n.array(
                            [self.nextStack] + [fname] + \
                            [float(x) for x in offset.split(" ")] + \
                            [float(x) for x in spacing.split(" ")] + \
                            [parent], dtype=n.object)
                        self.stackFrame.AddRow(stack_row)
                        self.nextStack += 1
                    parent_dialog.Destroy()
                space_dialog.Destroy()
            offset_dialog.Destroy()
        file_dialog.Destroy()        

        
if __name__ == "__main__":
    app = StackHandler()
    app.MainLoop()
