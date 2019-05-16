""" Utility code for Talk Figures
"""

def set_fontsize(ax,fsz):
    '''
    Generate a Table of columns and so on
    Restrict to those systems where flg_clm > 0

    Parameters
    ----------
    ax : Matplotlib ax class
    fsz : float
      Font size
    '''
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                     ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(fsz)

def set_spines(ax, lw):
    """ Set axis line widths
    Parameters
    ----------
    ax
    lw

    Returns
    -------

    """
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(lw)
