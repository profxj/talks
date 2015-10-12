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
from astropy import units as u
from astropy import constants as const

from xastropy.cgm.cos_halos import COSHalos
from xastropy.spec.lines_utils import AbsLine
from xastropy.obs import finder
from xastropy.plotting import utils as xputils
from xastropy.xutils import xdebug as xdb

# Local
#sys.path.append(os.path.abspath("../Analysis/py"))
#import lls_sample as lls_s


####
#  Series of plots illustrating COS-Halos experiment
def fig_experiment(outfil=None):

    # Init COS-Halos sightline
    cos_halos = COSHalos()
    cos_halos.load_single( ('J0950+4831','177_27'))
    cgm_abs = cos_halos.cgm_abs[0]

    # ########################################
    # Finder (out of order to avoid PDF issues)
    finder.main([cgm_abs.name, cgm_abs.galaxy.coord], imsize=2.*u.arcmin,
        show_circ=False)

    # Start the plot
    if outfil is None:
        outfil='fig_experiment.pdf'
    pp = PdfPages(outfil)

    # Lya spec
    lclr = 'blue'
    for ss in range(3):

        plt.figure(figsize=(8, 5))
        plt.clf()
        gs = gridspec.GridSpec(1, 1)

        # Axes
        if ss == 0:
            spec = cos_halos.load_spec(0, 1215.6700*u.AA)
            wvmnx = np.array([1455., 1490.])*u.AA
            scl_wv = 1.
        elif ss == 1:
            scl_wv = (1+cgm_abs.galaxy.z)
        elif ss == 2:
            spec = cos_halos.load_spec(0, 1025.7222*u.AA)
            scl_wv = (1+cgm_abs.galaxy.z)
            wvmnx = np.array([1015., 1037.])*u.AA * scl_wv # convoluted, yes

        ax = plt.subplot(gs[0,0])
        #ax.xaxis.set_minor_locator(plt.MultipleLocator(0.5))
        #ax.xaxis.set_major_locator(plt.MultipleLocator(1.))
        #ax.get_xaxis().get_major_formatter().set_useOffset(False)
        ax.set_xlim(wvmnx.value/scl_wv)
        ax.set_ylim(-0.1, 1.3)
        ax.set_ylabel('Normalized Flux')

        # Zero line
        ax.plot(wvmnx.value/scl_wv, (0.,0.), 'g--')

        # Data
        ax.plot(spec.dispersion/scl_wv, spec.flux, 'k')
        ax.plot(spec.dispersion/scl_wv, spec.sig, 'r:')

        # Label
        if ss == 0:
            ax.set_xlabel(r'Wavelength ($\AA$)')
            ax.text(1480., 0.4, cgm_abs.field, ha='left', fontsize=21., color=lclr)
        elif ss == 1:
            ax.set_xlabel(r'Rest Wavelength ($\AA$)')
            ax.text(1220., 0.4, r'HI Ly$\alpha$'+' \n '+r'$z=${:0.4f}'.format(cgm_abs.galaxy.z),
                ha='left', fontsize=21., multialignment='center', color=lclr)
        elif ss == 2:
            ax.set_xlabel(r'Rest Wavelength ($\AA$)')
            ax.text(1017., 0.4, r'HI Ly$\beta$'+' \n '+r'$z=${:0.4f}'.format(cgm_abs.galaxy.z),
                ha='left', fontsize=21., multialignment='center', color=lclr)

        # Fonts
        xputils.set_fontsize(ax,17.)

        # Write
        plt.tight_layout(pad=0.2,h_pad=0.,w_pad=0.1)
        pp.savefig()
        plt.close()


    # ########################################
    # Simple stack plot
    Zion = (1,1) # HI
    trans = np.array([1215.6701,  1025.7223])*u.AA
    lbls = [r'HI Ly$\alpha$', r'HI Ly$\beta$']

    plt.figure(figsize=(8, 5))
    plt.clf()
    gs = gridspec.GridSpec(2, 1)

    for qq,itran,lbl in zip(range(len(trans)),trans,lbls):
        # Load
        spec = cos_halos.load_spec(0, itran)
        # Axes
        ax = plt.subplot(gs[qq,0])
        ax.set_xlim(-1000., 1000.)
        ax.set_ylim(-0.1, 1.3)
        ax.set_ylabel('Normalized Flux')
        if qq == 0:
            ax.xaxis.set_ticklabels([])
        else:
            ax.set_xlabel('Relative Velocity (km/s)')
        # Velo
        velo = spec.relative_vel(itran*(cgm_abs.galaxy.z+1))
        ax.plot(velo, spec.flux, 'k', drawstyle='steps', lw=1.3)
        # Label
        ax.text(-800., 0.2, lbl, ha='left', fontsize=21., color=lclr)
        # Fonts
        xputils.set_fontsize(ax,17.)
    # Write
    plt.tight_layout(pad=0.2,h_pad=0.,w_pad=0.1)
    pp.savefig()
    plt.close()

    # Finish
    print('tlk_coshalos: Wrote {:s}'.format(outfil))
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
        fig_experiment()

    # pLLS
    if (flg_fig % 2**2) >= 2**1:
        fig_cog_abs_web()


# Command line execution
if __name__ == '__main__':

    if len(sys.argv) == 1: # Figs
        flg_fig = 0 
        flg_fig += 2**0 # Experiment
        #flg_fig += 2**1 # COG
        #flg_fig += 2**1 # pLLS
        #flg_fig += 2**2 # Ambig
        #flg_fig += 2**3 # Summary
        #flg_fig += 2**4 # Check NHI
    else:
        flg_fig = sys.argv[1]

    main(flg_fig)
