# Module for the Spectra figures for talks
# Imports
from __future__ import print_function

import numpy as np
import glob, os, sys
from IPython import embed

import matplotlib as mpl
mpl.rcParams['font.family'] = 'stixgeneral'
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec

from astropy import units
from astropy import constants
from astropy.cosmology import Planck15

from linetools.spectra import io as lsio


# Local
sys.path.append(os.path.abspath("../py"))
import tlk_fig_utils


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


def fig_age_T_BBN(outfile='fig_age_T_BBN.png'):
    """
    """
    tlk_fig_utils.set_mplrc()

    z = 10 ** np.linspace(7., 14, 10000)
    # Temperature
    Tcmb = Planck15.Tcmb(z)

    # Ages
    age = 152 * (0.1 * units.MeV / constants.k_B / Tcmb) * units.s
    age = age.to('s')

    xrange = [age[-1].value, 2e3]
    yrange = [1e8, 1e15]

    # Figure
    plt.figure(figsize=(5, 5))
    plt.clf()
    gs = gridspec.GridSpec(1, 1)
    ax = plt.subplot(gs[0])

    # Temperature
    ax.plot(age.to('s'), Tcmb, 'k-')

    # Axes
    ax.set_xlabel(r'\textbf{Age (seconds)}')
    ax.set_ylabel(r'\textbf{Temperature (Kelvin)}')
    ax.set_xlim(xrange)
    ax.set_ylim(yrange)

    ax.set_yscale("log", nonposy='clip')
    ax.set_xscale("log", nonposx='clip')

    lsz = 15.

    # Quark soup
    iquark = np.argmin(np.abs(Tcmb.value-1e13))
    ax.fill_between(age[iquark:].value,
                    yrange[0], Tcmb[iquark:].value, color='green', alpha=0.4)
    ax.text(1.5e-3, 2e14, r'\textbf{Quark}'+'\n'+r'\textbf{Soup}', color='black',
            fontsize=lsz, ha='left', va='bottom')

    # Electrons decouple
    iBBN0 = np.argmin(np.abs(Tcmb.value-9e8))
    ax.fill_between(age[iBBN0:iquark].value,
                    yrange[0], Tcmb[iBBN0:iquark].value, color='red', alpha=0.4)
    ax.text(2e-2, 1e13, r'\textbf{p,n freeze out}', color='black',
            fontsize=lsz, ha='left', va='bottom')
    ax.text(8, 10**10.5, r'$\nu \, \rm decouple$', color='black',
            fontsize=lsz, ha='left', va='bottom')
    ax.text(20, 10**10, r'$\rm e^-$', color='black',
            fontsize=lsz, ha='left', va='bottom')

    # BBN
    iBBN1 = np.argmin(np.abs(Tcmb.value-2e8))
    ax.fill_between(age[iBBN1:iBBN0].value,
                    yrange[0], Tcmb[iBBN1:iBBN0], color='blue', alpha=0.4)
    ax.text(220, 10**9, r'\textbf{BBN}', color='black',
            fontsize=lsz, ha='left', va='bottom')

    # Label
    #lsz = 18.
    #ax.text(4259.5, 0.25, r'HI Ly$\alpha$', color='blue', ha='center', size=lsz)
    #ax.text(4257.95, 0.80,r'DI Ly$\alpha$', color='darkgreen', size=lsz, ha='left')
    set_fontsize(ax, 16.)

    # Write
    plt.tight_layout(pad=0.2,h_pad=0.,w_pad=0.1)
    plt.savefig(outfile, dpi=500)
    plt.close()
    print("Wrote: {:s}".format(outfile))



def fig_age_CMB(outfile='fig_age_CMB.png'):
    """
    Cartoon the T evolution into the CMB release
    """
    tlk_fig_utils.set_mplrc()

    z = 10 ** np.linspace(np.log10(6.), 4, 10000)

    # Temperature
    Tcmb = Planck15.Tcmb(z)

    # Ages
    age = Planck15.age(z)
    age = age.to('year')

    xrange = [1e4, np.max(age.value)]
    yrange = [10, 2e4]

    # Figure
    plt.figure(figsize=(5, 5))
    plt.clf()
    gs = gridspec.GridSpec(1, 1)
    ax = plt.subplot(gs[0])

    # Temperature
    ax.plot(age.to('year'), Tcmb, 'k-')

    # Axes
    ax.set_xlabel(r'\textbf{Age (year)}')
    ax.set_ylabel(r'\textbf{Temperature (Kelvin)}')
    ax.set_xlim(xrange)
    ax.set_ylim(yrange)

    ax.set_yscale("log", nonposy='clip')
    ax.set_xscale("log", nonposx='clip')

    lsz = 15.

    # Radiation dominated
    irad = np.argmin(np.abs(Tcmb.value-1e4))
    ax.fill_between(age[irad:].value,
                    yrange[0], Tcmb[irad:], color='green', alpha=0.4)
    ax.text(1.2e4, 1e4, r'\textbf{Radiation}'+'\n'+r'\textbf{Dominated}', color='black',
            rotation=90., fontsize=lsz, ha='left', va='top')

    # Photons Decouple
    iCMB = np.argmin(np.abs(Tcmb.value-10**3.5))
    ax.fill_between(age[iCMB:irad].value,
                    yrange[0], Tcmb[iCMB:irad], color='red', alpha=0.4)
    ax.text(3e5, 10**3.5, r'\textbf{CMB}'+'\n'+r'\textbf{Released}', color='black',
            rotation=90., fontsize=lsz, ha='right', va='top')

    # Dark ages
    ax.fill_between(age[0:iCMB].value, yrange[0], Tcmb[0:iCMB], color='blue', alpha=0.4)
    ax.text(1e7, 60, r'\textbf{Dark}'+'\n'+r'\textbf{Ages}', color='black',
            fontsize=lsz, ha='left', va='top')

    # Label
    #lsz = 18.
    #ax.text(4259.5, 0.25, r'HI Ly$\alpha$', color='blue', ha='center', size=lsz)
    #ax.text(4257.95, 0.80,r'DI Ly$\alpha$', color='darkgreen', size=lsz, ha='left')
    set_fontsize(ax, 16.)

    # Write
    plt.tight_layout(pad=0.2,h_pad=0.,w_pad=0.1)
    plt.savefig(outfile, dpi=500)
    plt.close()
    print("Wrote: {:s}".format(outfile))


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
    if (flg_fig & 2**0):
        fig_q1009(outfil='fig_EW.pdf')

    # T during BBN
    if (flg_fig & 2**1):
        fig_age_T_BBN()

    # T during CMB
    if (flg_fig & 2**2):
        fig_age_CMB()



# Command line execution
if __name__ == '__main__':

    if len(sys.argv) == 1: # Figs
        flg_fig = 0 
        #flg_fig += 2**0 # D in Q1009
        flg_fig += 2**1 # age vs T
        flg_fig += 2**2 # age vs T near CMB release
    else:
        flg_fig = sys.argv[1]

    main(flg_fig)
