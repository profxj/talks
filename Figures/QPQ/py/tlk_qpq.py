# Module for the CGM Kinematics
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
from astropy.units import Quantity
from astropy import constants as const
from astropy.coordinates import SkyCoord

from linetools.lists.linelist import LineList
from linetools.spectra.xspectrum1d import XSpectrum1D

from xastropy.plotting import simple as xpsimp
from xastropy.plotting import utils as xputils
from xastropy.xutils import xdebug as xdb

from enigma.qpq import utils as equ
from enigma.qpq import spec as eqs

# Local
sys.path.append(os.path.abspath("../Analysis/py"))
#import lls_sample as lls_s


####
#  CGM Kinematic Measures
def fig_qpq_sample(qpq7=None, outfil=None):

    # Load QPQ7
    if qpq7 is None:
        qpq7 = equ.load_qpq(7)
    bg_coord = SkyCoord(ra=qpq7['RAD_BG'], dec=qpq7['DECD_BG'], unit=u.deg)
    fg_coord = SkyCoord(ra=qpq7['RAD'], dec=qpq7['DECD'], unit=u.deg)

    # Small sep
    gdsep = np.where(qpq7['R_PHYS'] < 300.)[0]
    ngd = len(gdsep)
    print('Nclose = {:d}'.format(ngd))

    # Get PA, x, y
    PA = []
    for ii,gds in enumerate(gdsep):
        PA.append(fg_coord[gds].position_angle(bg_coord[gds]))
    PA = Quantity(PA)
    x = qpq7['R_PHYS'][gdsep] * np.cos(PA)
    y = qpq7['R_PHYS'][gdsep] * np.sin(PA)

    if outfil is None:
        outfil='fig_qpq_sample.pdf'
    pp = PdfPages(outfil)

    fig = plt.figure(figsize=(6.5, 5))
    plt.clf()
    gs = gridspec.GridSpec(1, 1)
    xpsimp.dark_bkgd(plt)

    jet = cm = plt.get_cmap('jet')

    ax = plt.subplot(gs[0,0])
    ax.set_frame_on(False)
        #ax.xaxis.set_minor_locator(plt.MultipleLocator(0.5))
        #ax.xaxis.set_major_locator(plt.MultipleLocator(1.))
        #ax.get_xaxis().get_major_formatter().set_useOffset(False)
    ax.set_xlim(-300, 300)
    ax.set_ylim(-300, 300)
    #ax.set_xlabel('Relative Velocity (km/s)')
    #ax.set_ylabel('Normalized Flux')
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])

    mplt = ax.scatter(x, y, c=qpq7['L_BOL'][gdsep], cmap=jet,
                      edgecolor='none')

    mplt.set_clim(vmin=45, vmax=47.)
    cb = fig.colorbar(mplt)
    cb.set_label(r'$\log \; L_{\rm Bol}$')

    # Fonts

    # Write
    plt.tight_layout(pad=0.2,h_pad=0.,w_pad=0.1)
    pp.savefig()
    pp.close()
    plt.close()
    print('Genereated {:s}'.format(outfil))



def fig_cii_star(qpq7=None, outfil=None):
    """ Simple CII* plot for two systemsj
    """
    # Spec files

    # Load QPQ7
    #if qpq7 is None:
    #    qpq7 = equ.load_qpq(7)

    if outfil is None:
        outfil='fig_qpq_ciistar.pdf'
    pp = PdfPages(outfil)

    fig = plt.figure(figsize=(6.5, 5))
    plt.clf()
    gs = gridspec.GridSpec(2, 2)
    xpsimp.dark_bkgd(plt)

    qpq8_sys = ['J142758.73-012136.1', 'J120416.68+022110.9']
    #qpq8_sys = ['J120416.68+022110.9']
    qpq8_zfg = [2.2736, 2.4358]
    qpq8_vmnx = [(700,880), (600,770)]
    trans = [1334.5323, 1335.7077]
    lbls = ['CII 1334', 'CII* 1335']

    for qq,qsys in enumerate(qpq8_sys):
        # Load spectrum
        sdict = eqs.spec_wvobs(qsys, trans[0]*(1+qpq8_zfg[qq])*u.AA, high_res=2)
        spec = sdict['spec']
        spec.normalize(sdict['conti'])
        # Loop on transitions
        for jj,itrans in enumerate(trans):

            ax = plt.subplot(gs[jj,qq])
            ax.set_xlim(qpq8_vmnx[qq])
            ax.xaxis.set_major_locator(plt.MultipleLocator(50.))
            ax.set_ylim(-0.1, 1.2)
            ax.set_ylabel('Normalized Flux')
            if jj == 1:
                ax.set_xlabel('Relative Velocity (km/s)')
            else:
                ax.get_xaxis().set_ticks([])

            # Plot
            velo = spec.relative_vel(itrans*u.AA*(1+qpq8_zfg[qq]))
            ax.plot(velo, spec.flux, color='white', drawstyle='steps')
            ax.plot(qpq8_vmnx[qq],[1.]*2, '--', color='cyan')
            # Label
            ax.text(0.80, 0.10, lbls[jj], transform=ax.transAxes, size='large',
                    ha='center', color='yellow')

    # Write
    plt.tight_layout(pad=0.2,h_pad=0.,w_pad=0.1)
    pp.savefig()
    pp.close()
    plt.close()
    print('Genereated {:s}'.format(outfil))


def fig_molecules():
    """ H2, Lya, CII plots
    """
    # Spec files

    # Load QPQ7
    #if qpq7 is None:
    #    qpq7 = equ.load_qpq(7)

    # LineList
    H2 = LineList('H2')
    CO_waves = [1447.3521, 1477.5649, 1509.7480]
    CO_lbls = ['CO 2-0', 'CO 1-0', 'CO 0-0']

    outfil='fig_qpq_molecules.pdf'
    pp = PdfPages(outfil)

    fig = plt.figure(figsize=(8.5, 4))
    plt.clf()
    gs = gridspec.GridSpec(2, 5)
    xpsimp.dark_bkgd(plt)

    qpq8_sys = ['J114436.65+095904.9', 'J142758.73-012136.1']
    qpq8_yrng = [(0., 1.2), (0.,1.2)]
    wvobs = [(4160,4240.), (4735,4950)]  # Lyman-Werner, CO
    #qpq8_sys = ['J120416.68+022110.9']
    qpq8_zfg = [2.973, 2.27616]
    trans = [1215.6701, 1334.5323]
    tvmnx = [(-400,400), (-200,200)]
    lbls = ['HI Lya', 'CII 1334']

    for qq,qsys in enumerate(qpq8_sys):
        # Load spectrum
        sdict = eqs.spec_wvobs(qsys, wvobs[qq][0]*u.AA, high_res=2)
        spec = sdict['spec']
        spec.normalize(sdict['conti'])
        if qq == 1:  # Smooth for presentation
            tmpspec = XSpectrum1D.from_tuple((spec.wavelength, spec.flux))
            spec = tmpspec.box_smooth(3)

        # Loop on standard transitions
        for jj,itrans in enumerate(trans):

            ax = plt.subplot(gs[qq,jj])
            ax.set_xlim(tvmnx[jj])
            ax.xaxis.set_major_locator(plt.MultipleLocator(200.))
            ax.set_ylim(-0.05, 1.2)
            if jj == 0:
                ax.set_ylabel('Flux')
            ax.set_xlabel('Velocity (km/s)')

            # Plot
            velo = spec.relative_vel(itrans*u.AA*(1+qpq8_zfg[qq]))
            ax.plot(velo, spec.flux, color='white', drawstyle='steps')
            ax.plot(tvmnx[jj],[1.]*2, '--', color='lightgreen')
            # Label
            ax.text(0.50, 0.90, lbls[jj], transform=ax.transAxes, size='large',
                    ha='center', color='yellow')
        # Molecules
        ax = plt.subplot(gs[qq,len(trans):])
        ax.set_xlim(wvobs[qq])
        ax.set_ylim(qpq8_yrng[qq])
        ax.set_xlabel('Wavelength (Angstroms)')
        ax.plot(spec.wavelength, spec.flux, color='white', drawstyle='steps')
        if qq == 0: # Lyman-Werner
            wrlim = [iwv*u.AA/(1+qpq8_zfg[qq]) for iwv in wvobs[qq]]
            gdlin = np.where((H2._fulltable['Jk']<7) & (H2._fulltable['wrest']>wrlim[0]) & (H2._fulltable['wrest']<wrlim[1]))[0]
            for igd in gdlin:
                if H2._fulltable['Jk'][igd] <= 2:
                    clr = 'cyan'
                    do_lbl = True
                else:
                    clr = 'gray'
                    do_lbl = False
                ax.plot([H2._fulltable['wrest'][igd].value*(1+qpq8_zfg[qq])]*2,
                        qpq8_yrng[qq], '--', color=clr)
                if do_lbl:
                    ax.text(H2._fulltable['wrest'][igd].value*(1+qpq8_zfg[qq]),
                            0.4, H2._fulltable['name'][igd], rotation=90., color=clr,
                            size=12)
        else: # CO
            for jj,CO_wave in enumerate(CO_waves):
                ax.plot([CO_wave*(1+qpq8_zfg[qq])]*2, qpq8_yrng[qq], '--',
                        color='cyan')
                ax.text(CO_wave*(1+qpq8_zfg[qq]), 0.4, CO_lbls[jj],
                        rotation=90, color='cyan', size=12)


    # Write
    plt.tight_layout(pad=0.2,h_pad=0.,w_pad=0.1)
    pp.savefig()
    pp.close()
    plt.close()
    print('Genereated {:s}'.format(outfil))





# Main
def main(flg_fig):
    # Init
    if flg_fig == 'all':
        flg_fig = np.sum( np.array( [2**ii for ii in range(1)] )) # Skipping fig_chk_nhi
    else:
        flg_fig = int(flg_fig)

    # Experimental design
    if (flg_fig % 2**1) >= 2**0:
        fig_qpq_sample()

    # CII* figure
    if (flg_fig % 2**2) >= 2**1:
        fig_cii_star()

    # Molecules
    if (flg_fig % 2**3) >= 2**2:
        fig_molecules()


# Command line execution
if __name__ == '__main__':

    if len(sys.argv) == 1: # Figs
        flg_fig = 0 
        #flg_fig += 2**0 # QPQ sample
        #flg_fig += 2**1 # CII*
        flg_fig += 2**2 # Molecules
        #flg_fig += 2**2 # Abs gallery
    else:
        flg_fig = sys.argv[1]

    main(flg_fig)
