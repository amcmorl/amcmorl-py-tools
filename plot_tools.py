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

def get_ax_rect(i_ax, n_cols, n_rows, margin=Margins(), order='across-first'):
    '''Calculate rect values for axis.
    '''
    i_axs = np.array(i_ax)
    axs = get_ax_rect(i_axs, n_cols, n_rows, margin=margin, order=order)
    return axs[0]
        
def get_ax_rects(i_axs, n_cols, n_rows, margin=Margins(), order='across-first'):
    '''Calculate rect values for several axes.
    '''
    i_axs = np.asarray(i_axs)    
    assert order in ['across-first', 'down-first']
    w = (1 - (margin.left + margin.right + \
                  (n_cols - 1) * margin.hgap)) / float(n_cols)
    h = (1 - (margin.bottom + margin.top + \
                  (n_rows - 1) * margin.vgap)) / float(n_rows)
    if order == 'across-first':
        n_col = i_axs % n_cols
        n_row = i_axs / n_cols
    if order == 'down-first':
        n_col = i_axs / n_rows
        n_row = i_axs % n_rows
    l = margin.left + n_col * (w + margin.hgap)
    b = margin.bottom + (n_rows - n_row - 1) * (h + margin.vgap)
    ax_rects = np.vstack((l, b, np.ones_like(l) * w, np.ones_like(b) * h))
    return ax_rects.T

def create_plot_grid(n_axes, n_cols=1, margin=Margins(), fig=None,
                     order='across-first' ):
    '''Create a grid of axes for suitable for plots.

    Parameters
    ----------
    n_axes, n_cols : int
      number of axes and columns
    margin : Margins instance
    fig : mpl figure
      figure to use; if None, creates a new figure
    order : {'across-first', 'down-first'}
      numbering order for plots

    Returns
    -------
    axes : list of mpl Axes objects
    '''
    n_rows = n_axes / n_cols if n_axes % n_cols == 0 \
        else n_axes / n_cols + 1
    axes = []
    i_axs = np.arange(n_axes)
    axrects = get_ax_rects(i_axs, n_cols, n_rows, margin=margin, order=order)
    if fig == None:
        fig = plt.figure()
    for axrect in axrects:
        ax = fig.add_axes(axrect)
        which = []
        if n_col == 0:
            which.append('left')
        if n_row == (n_rows - 1):
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
