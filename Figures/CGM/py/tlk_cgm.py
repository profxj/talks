# Module for the Spectra figures for talks
# Imports
from __future__ import print_function

import numpy as np
import glob, os, sys

import matplotlib as mpl
mpl.rcParams['font.family'] = 'stixgeneral'
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec

from astropy.io import ascii
from astropy.relativity import velocities as arv
from astropy import units as u
from astropy import constants as const

from xastropy.spec.lines_utils import AbsLine
from xastropy.igm.abs_sys.lls_utils import LLS_System, LLS_Survey
from xastropy.igm.abs_sys import abssys_utils as abssys
from xastropy import spec as xsp
from xastropy.plotting import simple as xplots
from xastropy.plotting import utils as xputils
from xastropy.xutils import xdebug as xdb

# Local
#sys.path.append(os.path.abspath("../Analysis/py"))
#import lls_sample as lls_s


####
#  Simple plot showing EW
def fig_ew(outfil=None):


    # Init
    wvmnx = (1214.5, 1216.8)*u.AA
    wave = np.linspace(wvmnx[0].value,wvmnx[1].value,100)
    wave = wave * u.AA

    # Single line (Lya)
    zabs = 0.0
    line = AbsLine(1215.6701*u.AA)
    line.z = zabs
    line.attrib['N'] = 16.0
    line.attrib['b'] = 30.0

    vmodel = xsp.voigt.voigt_model(wave, line, Npix=None)  # No smoothing
    line.spec = vmodel

    # EW
    line.analy['WVMNX'] = wvmnx
    EW, sigEW = line.ew()

    # Plot boundaries
    ymnx = (-0.05, 1.1)

    if outfil == None:
        outfil='fig_EW.pdf'

    # Start the plot
    if outfil != None:
        pp = PdfPages(outfil)

    # Make the plots
    for ss in range(2):

        plt.figure(figsize=(8, 5))
        plt.clf()
        gs = gridspec.GridSpec(1, 1)

        # Axes
        ax = plt.subplot(gs[0,0])
        ax.xaxis.set_minor_locator(plt.MultipleLocator(0.5))
        ax.xaxis.set_major_locator(plt.MultipleLocator(1.))
        ax.get_xaxis().get_major_formatter().set_useOffset(False)
        ax.set_xlim(wvmnx.value)
        ax.set_ylim(ymnx)
        ax.set_ylabel('Normalized Flux')
        ax.set_xlabel(r'Wavelength ($\AA$)')

        # Zero line
        ax.plot(wvmnx.value, (0.,0.), 'g--')

        # Data
        lines = ax.step(vmodel.dispersion, vmodel.flux, 'k', where='mid', 
            linewidth=1.5)
        #lines = ax.step(vmodel.dispersion, vmodel.flux, 'k', where='mid')#, linewidth=1.5)

        # Fill?
        if ss > 0:
            #ax.fill_between( vmodel.dispersion, vmodel.flux, [1.]*len(wave),
            #                        color='blue', alpha=0.3)
            ax.fill_between(lines[0].get_xdata(orig=False), 1, 
                lines[0].get_ydata(orig=False), alpha=0.4)
            csz = 20.
            ax.text(1215.670, 0.5, 'EW', color='black', ha='center',
                size=csz)
            ax.text(1214.8, 0.2, 
                r'$W_\lambda = \, $'+'{:.2f}'.format(EW.value)+r' $\AA$',
                color='blue', ha='left', size=csz)

        # Fonts
        xputils.set_fontsize(ax,15.)

        # Write
        plt.tight_layout(pad=0.2,h_pad=0.,w_pad=0.1)
        pp.savefig()
        plt.close()

    # Finish
    print('tlk_spectra.fig_ew: Wrote {:s}'.format(outfil))
    pp.close()

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
        fig_ew(outfil='fig_EW.pdf')

    # pLLS
    if (flg_fig % 2**2) >= 2**1:
        fig_cog_abs_web()


# Command line execution
if __name__ == '__main__':

    if len(sys.argv) == 1: # Figs
        flg_fig = 0 
        #flg_fig += 2**0 # EW
        flg_fig += 2**1 # COG
        #flg_fig += 2**1 # pLLS
        #flg_fig += 2**2 # Ambig
        #flg_fig += 2**3 # Summary
        #flg_fig += 2**4 # Check NHI
    else:
        flg_fig = sys.argv[1]

    main(flg_fig)
