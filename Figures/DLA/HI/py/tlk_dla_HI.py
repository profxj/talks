# Module for DLA HI figures
# Imports
from __future__ import print_function


import numpy as np
import glob, os, sys
import pdb

import matplotlib as mpl
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['mathtext.fontset'] = 'cm'
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec

from astropy import units as u
from astropy import constants as const

from linetools.lists.linelist import LineList
from linetools.analysis.voigt import voigt_from_abslines
from linetools.spectralline import AbsLine

from pyigm.surveys.dlasurvey import DLASurvey

# Local
sys.path.append(os.path.abspath("../../py"))
from tlk_fig_utils import set_fontsize, set_spines

def fig_lya_line(lw=1.5, csz=15.):
    """  Generate a DLA in optical depth and flux space
    Parameters
    ----------
    """
    llist = LineList('ISM')
    # Lya
    lya = AbsLine('HI 1215', z=0., llist=llist)

    # Wavelength
    wave = np.linspace(1100., 1310., 100000)*u.AA

    outfile = 'fig_lya_line.png'

    # Figure
    plt.figure(figsize=(8, 4))
    plt.clf()
    gs = gridspec.GridSpec(1, 2)


    # Tau plot
    ax1 = plt.subplot(gs[0])

    # Optical depth
    lya.attrib['N'] = 1e21 / u.cm**2
    lya.attrib['b'] = 20 * u.km/u.s
    tau = voigt_from_abslines(wave, lya, ret='tau')

    # Plot
    ax1.plot(wave, tau, 'g')

    # Axes
    ax1.set_xlim(1200., 1230.)
    ax1.set_ylim(1e-2, 5e7)
    ax1.set_yscale("log", nonposy='clip')
    ax1.set_ylabel(r'Optical Depth for $N_{\rm HI} = 10^{21} \, \rm cm^{-2}$')
    ax1.set_xlabel(r'Wavelength ($\AA$)')
    ax1.xaxis.set_major_locator(plt.MultipleLocator(10.))
    #
    set_spines(ax1, 2.)
    set_fontsize(ax1,csz)

    # Flux space
    ax2 = plt.subplot(gs[1])

    for logN in [19., 20., 21., 22.]:
        lya.attrib['N'] = 10**logN / u.cm**2
        spec = voigt_from_abslines(wave, lya)
        # Plot
        ax2.plot(spec.wavelength, spec.flux, label=r'$\log N_{\rm HI} = $'+'{:d}'.format(int(logN)))
    # Axes
    wvoff = 100.
    ax2.set_xlim(1215.-wvoff, 1215.+wvoff)
    ax2.set_ylim(-0.05, 1.1)
    #ax2.set_yscale("log", nonposy='clip')
    ax2.set_ylabel('Normalized Flux')
    ax2.set_xlabel(r'Wavelength ($\AA$)')
    ax2.xaxis.set_major_locator(plt.MultipleLocator(50.))
    #
    set_spines(ax2, 2.)
    set_fontsize(ax2,csz)
    legend = ax2.legend(loc='lower left', scatterpoints=1, borderpad=0.3,
                      handletextpad=0.3, fontsize='small', numpoints=1)

    # Write
    plt.tight_layout(pad=0.2,h_pad=0.,w_pad=0.1)
    plt.savefig(outfile, dpi=750)
    plt.close()
    print("Wrote {:s}".format(outfile))


def fig_rhoHI(lw=1.5, csz=15., lsz=14.):
    """  Generate a DLA in optical depth and flux space
    Parameters
    ----------
    """
    sdss = DLASurvey.load_SDSS_DR5()
    zbins = [2.2, 2.4, 2.75, 3., 3.5, 4.5]
    rho_HI, rho_HI_low, rho_HI_hi = sdss.binned_rhoHI(zbins)

    outfile = 'fig_rhoHI.png'

    # Figure
    plt.figure(figsize=(5, 5))
    plt.clf()
    gs = gridspec.GridSpec(1, 1)


    # Tau plot
    ax = plt.subplot(gs[0])

    # Plot
    for kk in range(len(zbins)-1):
        zcen = np.sum(zbins[kk:kk+2])/2.
        yerr= np.array([rho_HI_low[kk].value/1e8, rho_HI_hi[kk].value/1e8])
        ax.errorbar([zcen], [rho_HI[kk].value/1e8], xerr=zcen-zbins[kk], fmt='o', color='blue', capthick=2)
        ax.errorbar([zcen], [rho_HI[kk].value/1e8], yerr=[yerr], color='blue', capthick=2)
    # z=0
    xmnx = [2., 4.5]
    ax.fill_between(xmnx, 0.45, 0.6, color='green', alpha=0.5)

    # Axes
    ax.set_xlim(xmnx)
    #ax.set_ylim(1e-2, 5e7)
    ax.set_ylabel(r'$\rho_{\rm HI} \; (10^8 \, \rm M_\odot \, Mpc^{-3} \, h_{72})$')
    ax.set_xlabel(r'$z$')
    ax.text(0.1, 0.9, 'SDSS-DR5 (PW09)', color='blue', size=lsz, transform=ax.transAxes, ha='left')
    ax.text(0.9, 0.1, 'z~0 [21cm] \n (Zwaan+05)', color='green', size=lsz, transform=ax.transAxes, ha='right')
    #ax.xaxis.set_major_locator(plt.MultipleLocator(10.))
    #
    set_spines(ax, 2.)
    set_fontsize(ax,csz)

    # Write
    plt.tight_layout(pad=0.2,h_pad=0.,w_pad=0.1)
    plt.savefig(outfile, dpi=750)
    plt.close()
    print("Wrote {:s}".format(outfile))
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

    # Experimental design
    if flg_fig & 2**0:
        fig_lya_line()
    # rho_HI
    if flg_fig & 2**1:
        fig_rhoHI()

# Command line execution
if __name__ == '__main__':

    if len(sys.argv) == 1: # Figs
        flg_fig = 0
        #flg_fig += 2**0 # Lya line
        flg_fig += 2**1 # rho_HI
        #flg_fig += 2**1 # Image gallery
        #flg_fig += 2**2 # Abs gallery
        #flg_fig += 2**3  # Lya zoom-in
    else:
        flg_fig = sys.argv[1]

    main(flg_fig)
