import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid import AxesGrid

class Margins:
    def __init__(self, left=0, bottom=0, right=0, top=0, hgap=0, vgap=0):
        self.left = left
        self.bottom = bottom
        self.right = right
        self.top = top
        self.hgap = hgap
        self.vgap = vgap

    def __str__(self):
        str = "Left: %s\n" % (self.left)
        str += "Bottom: %s\n" % (self.bottom)
        str += "Right: %s\n" % (self.right)
        str += "Top: %s\n" % (self.top)
        str += "Horizontal gap: %s\n" % (self.hgap)
        str += "Vertical gap: %s" % (self.vgap)
        return str 

    def get_rect(self):
        return [self.left, self.bottom,
                1 - self.right - self.left, 1 - self.top - self.bottom]

    def set_from_rect(self, rect, total_width=1, total_height=1):
        self.left = rect[0]
        self.bottom = rect[1]
        width = rect[2]
        height = rect[3]
        self.right = total_width - width - self.left
        self.top = total_height - height - self.bottom

def make_margin_from_rect(rect, hgap=0, vgap=0,
                          total_width=1, total_height=1):
    left = rect[0]
    bottom = rect[1]
    width = rect[2]
    height = rect[3]
    right = total_width - width - left
    top = total_height - height - bottom
    return Margins(left, bottom, right, top, hgap, vgap)

def make_enough_rows(total, cols):
    rows = int(total / cols)
    extra = int(total % cols)
    if extra > 0:
        rows += 1
    return rows

def get_ax_rect(i_ax, ncols, nrows, margin=Margins(), direction='row'):
    '''Calculate rect values for axis.
    '''
    i_axs = np.array(i_ax)
    axs = get_ax_rect(i_axs, ncols, nrows, margin=margin, direction=direction)
    return axs[0]

def get_col_row(i, ncols, nrows, direction):
    if direction == 'row':
        ncol = i % ncols
        nrow = i / ncols
    if direction == 'col':
        ncol = i / nrows
        nrow = i % nrows
    return ncol, nrow

def get_ax_rects(i_axs, ncols, nrows, margin=Margins(), direction='row'):
    '''Calculate rect values for several axes.
    '''
    i_axs = np.asarray(i_axs)    
    assert direction in ['row', 'col']
    w = (1 - (margin.left + margin.right + \
                  (ncols - 1) * margin.hgap)) / float(ncols)
    h = (1 - (margin.bottom + margin.top + \
                  (nrows - 1) * margin.vgap)) / float(nrows)
    ncol, nrow = get_col_row(i_axs, ncols, nrows, direction)
    l = margin.left + ncol * (w + margin.hgap)
    b = margin.bottom + (nrows - nrow - 1) * (h + margin.vgap)
    ax_rects = np.vstack((l, b, np.ones_like(l) * w, np.ones_like(b) * h))
    return ax_rects.T

def create_plot_grid(n_axes, ncols=1, margin=Margins(), fig=None,
                     direction='row', sharex='none', sharey='none',
                     yspines='left', xspines='bottom'):
    '''Create a grid of axes suitable for plots.

    Parameters
    ----------
    n_axes, ncols : int
      number of axes and columns
    margin : Margins instance
    fig : mpl figure
      figure to use; if None, creates a new figure
    direction : {'row', 'col'}
      numbering direction for plots
    sharex : {'col', 'row', 'all', 'none'}
      how, if at all, to share x axes
    yspines : {'left', 'all'}
      where to draw y spines, relative to grid, 'all' means on all columns
    xspines : {'bottom', 'all'}
      where to draw x spines, relative to grid, 'all' mean on all rows
      
    Returns
    -------
    axes : list of mpl Axes objects
    '''
    assert(direction in ['row', 'col'])
    assert(sharex in ['col', 'row', 'all', 'none'])
    assert(yspines in ['left', 'all', 'none'])
    assert(xspines in ['bottom', 'all'])
    
    nrows = n_axes / ncols if n_axes % ncols == 0 \
        else n_axes / ncols + 1
    axes = []
    i_axs = np.arange(n_axes)
    axrects = get_ax_rects(i_axs, ncols, nrows, margin=margin,
                           direction=direction)
    if fig == None:
        fig = plt.figure()
    col_leader = None
    row_leader = None
    for i, axrect in enumerate(axrects):
        ncol, nrow = get_col_row(i, ncols, nrows, direction)
        if sharex == 'col':
            if (nrow == 0) and (ncol > 0):
                # reset at the top of new columns
                col_leader = None
        if sharey == 'row':
            if (ncol == 0) and (nrow > 0):
                # reset at the beginning of new rows
                row_leader = None
        ax = fig.add_axes(axrect, sharex=col_leader, sharey=row_leader)

        # axis sharing
        if sharex == 'col':
            if nrow == 0:
                col_leader = ax
        elif sharex == 'all':
            if (nrow == 0) and (ncol == 0):
                col_leader = ax
        
        if sharey == 'row':
            if ncol == 0:
                row_leader = ax
        elif sharey == 'all':
            if (ncol == 0) and (nrow == 0):
                row_leader = ax
                
        # spines
        which = []
        if yspines == 'left':
            if ncol == 0:
                which.append('left')
        elif yspines == 'all':
            which.append('left')
        if xspines == 'all':
            which.append('bottom')
        elif xspines == 'bottom':
            if nrow == (nrows - 1):
                which.append('bottom')
        format_spines(ax, which)
        axes.append(ax)
    return axes

def format_spines(ax, which=[], hidden_color='none'):
    '''
    Convenience function for formatting spines of a plot.

    Parameters
    ----------
    ax : matplotlib axis
    which : list of strings
      names of spines (defined by matplotlib axes object) to format
    hidden_color : string
      any matplotlib color specification
    '''
    for loc, spine in ax.spines.iteritems():
        if loc in which:
            spine.set_visible(True) # in case it was hidden previously
            spine.set_position(('outward', 5))
        else:
            if hidden_color != 'none':
                spine.set_color(hidden_color)
            else:
                spine.set_visible(False)

    if 'bottom' in which:
        ax.xaxis.set_ticks_position('bottom')
    else:
        ax.xaxis.set_ticks_position('none')
        for label in ax.get_xticklabels():
            label.set_visible(False)

    if 'left' in which:
        ax.yaxis.set_ticks_position('left')
    elif 'right' in which:
        ax.yaxis.set_ticks_position('right')
    else:
        ax.yaxis.set_ticks_position('none')
    if ('left' in which) or ('right' in which):
        for label in ax.get_yticklabels():
            label.set_visible(True)
    else:
        for label in ax.get_yticklabels():
            label.set_visible(False)

dashes = [[20, 2],
          [2, 2],
          [20, 2, 2, 2],
          [20, 2, 20, 2],
          [20, 2, 2, 2, 2, 2],
          [20, 2, 20, 2, 2, 2],
          [8, 2, 2, 2, 2, 2],
          [8, 2, 2, 2]]

colours = [[0., 0., 0.],
           [0., 0., 0.],
           [0.5, 0.5, 0.5],
           [0.5, 0.5, 0.5]]

class LineCycler():
    def __init__(self):
        self.c = 0
        self.d = 0

    def __call__(self, what='dashes'):
        if what == 'dashes' or what == 'd':
            n_styles = len(dashes)
            style = dashes[self.d % n_styles]
            self.d += 1
        else:
            n_styles = len(colours)
            style = colours[self.c % n_styles]
            self.c += 1
        return style

class FigNumer():
    def __init__(self):
        self.next_num = 0

    def __call__(self):
        next_num = self.next_num
        self.next_num += 1
        return next_num

def plot_panels(array, fig=None, nrows=1, panel_labels=None, extra_col=0.2, share_axes=True):
    '''
    Plot an array as a series of image panels.
    
    Parameters
    ----------
    array : array_like
      shape = n_panels, n_rows, n_cols
    nrows : int
      number of rows to plot panels in
    panel_labels : list of strings
      labels to give each panel
    extra_col : float
      "pad" around colour range
    '''
    if fig == None:
        fig = plt.figure()
    n_panels = array.shape[0]
    ncols = n_panels / nrows if ((n_panels % nrows) == 0) \
        else (n_panels / nrows) + 1
    grid = AxesGrid(fig, 111, nrows_ncols = (nrows, ncols),
                    axes_pad = 0.05, share_all=share_axes,
                    cbar_mode='single', cbar_location='right', cbar_size='15%')
    for i in xrange(n_panels):
        vmin = np.nanmin(array[0]) * (1 - extra_col)
        vmax = np.nanmax(array[0]) * (1 + extra_col)
        im = grid[i].imshow(array[i], cmap=mpl.cm.jet, vmin=vmin, vmax=vmax)
        if panel_labels != None:
            grid[i].set_title(panel_labels[i], size='small')
    plt.colorbar(im, cax=grid.cbar_axes[0])
    cax = grid.cbar_axes[0]
    cax.axis["right"].toggle(ticks=True, ticklabels=True, label=True)
    #cax.set_ylabel("")
    return grid
