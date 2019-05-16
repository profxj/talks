# Module for the Spectra figures for talks
# Imports
from __future__ import print_function

import numpy as np
import glob, os, sys

import matplotlib as mpl
mpl.rcParams['font.family'] = 'stixgeneral'
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec

from astropy import units as u
from astropy import constants as const

from linetools.spectra import io as lsio

#from xastropy.spec.lines_utils import AbsLine
#from xastropy.igm.abs_sys.lls_utils import LLS_System, LLS_Survey
#from xastropy.igm.abs_sys import abssys_utils as abssys
#from xastropy import spec as xsp
#from xastropy.plotting import simple as xplots
#from xastropy.plotting import utils as xputils
#from xastropy.xutils import xdebug as xdb

# Local
#sys.path.append(os.path.abspath("../Analysis/py"))
#import lls_sample as lls_s


####
#  Simple plot showing EW
def fig_q1009(outfil=None):
    """ 
    :param outfil: 
    :return: 
    """
    datfil = '~/Keck/HIRES/RedData/Q1009+2956/Q1009+2956a_f.fits'
    FWHM = 4.
    zhabs = 2.50369

    spec = lsio.readspec(datfil)

    plt.figure(figsize=(7, 5))
    plt.clf()
    gs = gridspec.GridSpec(1, 1)
    ax = plt.subplot(gs[0,0])

    # Lya first
    #x_pixminmax, wave, 1215.6701
    #d, zhabs, -8000, 8000, PIXMIN = pmn, $
    #PIXMAX = pmx, VELO = vel


    xr = [4256.5, 4263]
    yrange = [-0.05, 1.1]
    ax.plot(spec.wavelength, spec.flux, 'k', drawstyle='steps-mid')
    ax.set_xlabel('Wavelength (Ang)')
    ax.set_ylabel('Normalized Flux')
    ax.set_xlim(xr)
    ax.set_ylim(yrange)

    # Horizontal lines
    for yy in [0., 1.]:
        ax.axhline(yy, color='gray', linestyle='--')

    # Label
    lsz = 18.
    ax.text(4259.5, 0.25, r'HI Ly$\alpha$', color='blue', ha='center', size=lsz)
    ax.text(4257.95, 0.80,r'DI Ly$\alpha$', color='darkgreen', size=lsz, ha='left')
    set_fontsize(ax, 16.)
    ax.arrow(4258.15, 0.76, 0., -0.07, linewidth=2,
             head_width=0.1, head_length=0.03, fc='darkgreen', ec='darkgreen')

    #plotsym, 1, 2.5, thick = 6
    #oplot, [4258.15], [0.8], psym = 8, color = clr.darkgreen
    # Write
    plt.tight_layout(pad=0.2,h_pad=0.,w_pad=0.1)
    outfil = 'fig_q1009.png'
    plt.savefig(outfil, dpi=700)
    plt.close()
    print("Wrote: {:s}".format(outfil))

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


#### ########################## #########################
#### ########################## #########################
#### ########################## #########################

# Main
def main(flg_fig):
    # Init
    if flg_fig == 'all':
        flg_fig = np.sum( np.array( [2**ii for ii in range(1)] )) # Skipping fig_chk_nhi
    else:
        flg_fig = int(flg_fig)

    # EW
    if (flg_fig % 2**1) >= 2**0:
        fig_q1009(outfil='fig_EW.pdf')


# Command line execution
if __name__ == '__main__':

    if len(sys.argv) == 1: # Figs
        flg_fig = 0 
        flg_fig += 2**0 # D in Q1009
        #flg_fig += 2**1 # COG
    else:
        flg_fig = sys.argv[1]

    main(flg_fig)
