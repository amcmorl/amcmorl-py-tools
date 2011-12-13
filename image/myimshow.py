from pylab import get_current_fig_manager, gca, imshow

class MyFormatter:

    def __init__(self, img):
        self.img = img
        self.ypos = None
        self.xpos = None
        self.val = None

    def on_move(self, event):
        "save the last co-ords and get the index"
        if event.inaxes:
            # get the index into x,y closest to
            self.xpos = int(event.xdata)
            self.ypos = int(event.ydata)
            self.val = self.img[self.xpos,self.ypos]
        else:
            self.xpos = None
            self.ypos = None
            self.val = None
            
    def fmty(self, y):
        '''ignore y and format val of img at ind'''
        if (self.xpos is None) or (self.ypos is None): return ''
        return "value: %0.3f" % self.val

    def fmtx(self, x):
        '''ignore x - could use - and print x and y vals'''
        if (self.xpos is None) or (self.ypos is None): return ''
        return "x: %d y: %d" % (self.xpos, self.ypos)

def convert_ndarr( img ):
    return img[:,::-1].transpose()

def imshow_val( img, *args, **kwargs ):
    o = MyFormatter( img )
    imshow( convert_ndarr(img), *args, **kwargs )
    canvas = get_current_fig_manager().canvas
    ax = gca()
    canvas.mpl_connect( "motion_notify_event", o.on_move )
    ax.fmt_xdata = o.fmtx
    ax.fmt_ydata = o.fmty

    
