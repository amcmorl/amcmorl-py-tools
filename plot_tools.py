import numpy as np
import matplotlib.pyplot as plt

class Margins:
    def __init__(self, left=0, bottom=0, right=0, top=0, hgap=0, vgap=0):
        self.left = left
        self.bottom = bottom
        self.right = right
        self.top = top
        self.hgap = hgap
        self.vgap = vgap

def get_ax_rect(i_ax, ncols, nrows, margin=Margins(), order='across-first'):
    '''Calculate rect values for axis.
    '''
    i_axs = np.array(i_ax)
    axs = get_ax_rect(i_axs, ncols, nrows, margin=margin, order=order)
    return axs[0]

def get_col_row(i, ncols, nrows, order):
    if order == 'across-first':
        ncol = i % ncols
        nrow = i / ncols
    if order == 'down-first':
        ncol = i / nrows
        nrow = i % nrows
    return ncol, nrow

def get_ax_rects(i_axs, ncols, nrows, margin=Margins(), order='across-first'):
    '''Calculate rect values for several axes.
    '''
    i_axs = np.asarray(i_axs)    
    assert order in ['across-first', 'down-first']
    w = (1 - (margin.left + margin.right + \
                  (ncols - 1) * margin.hgap)) / float(ncols)
    h = (1 - (margin.bottom + margin.top + \
                  (nrows - 1) * margin.vgap)) / float(nrows)
    ncol, nrow = get_col_row(i_axs, ncols, nrows, order)
    l = margin.left + ncol * (w + margin.hgap)
    b = margin.bottom + (nrows - nrow - 1) * (h + margin.vgap)
    ax_rects = np.vstack((l, b, np.ones_like(l) * w, np.ones_like(b) * h))
    return ax_rects.T

def create_plot_grid(n_axes, ncols=1, margin=Margins(), fig=None,
                     order='across-first', sharex='none',
                     yspines='left', xspines='bottom'):
    '''Create a grid of axes for suitable for plots.

    Parameters
    ----------
    n_axes, ncols : int
      number of axes and columns
    margin : Margins instance
    fig : mpl figure
      figure to use; if None, creates a new figure
    order : {'across-first', 'down-first'}
      numbering order for plots
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
    assert(order in ['across-first', 'down-first'])
    assert(sharex in ['col', 'row', 'all', 'none'])
    assert(yspines in ['left', 'all'])
    assert(xspines in ['bottom', 'all'])
    
    nrows = n_axes / ncols if n_axes % ncols == 0 \
        else n_axes / ncols + 1
    axes = []
    i_axs = np.arange(n_axes)
    axrects = get_ax_rects(i_axs, ncols, nrows, margin=margin, order=order)
    if fig == None:
        fig = plt.figure()
    col_leader = None
    for i, axrect in enumerate(axrects):
        ncol, nrow = get_col_row(i, ncols, nrows, order)
        if sharex == 'col':
            if (nrow == 0) and (ncol > 0):
                # reset at the top of new columns
                col_leader = None
        ax = fig.add_axes(axrect, sharex=col_leader)

        # axis sharing
        if sharex == 'col':
            if nrow == 0:
                col_leader = ax
        elif sharex == 'all':
            if (nrow == 0) and (ncol == 0):
                col_leader = ax

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
            spine.set_position(('outward', 5))
        else:
            spine.set_color(hidden_color)
    if 'bottom' in which:
        ax.xaxis.set_ticks_position('bottom')
    else:
        ax.xaxis.set_ticks_position('none')
        ax.xaxis.set_ticklabels([])
    if 'left' in which:
        ax.yaxis.set_ticks_position('left')
    elif 'right' in which:
        ax.yaxis.set_ticks_position('right')
    else:
        ax.yaxis.set_ticks_position('none')
        ax.yaxis.set_ticklabels([])
