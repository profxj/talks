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

from astropy.io import ascii, fits
from astropy.relativity import velocities as arv
from astropy import units as u
from astropy import constants as const

from xastropy.spec.lines_utils import AbsLine
from xastropy.igm.abs_sys import abssys_utils as abssys
from xastropy import spec as xsp
from xastropy.plotting import simple as xplots
from xastropy.plotting import utils as xputils
from xastropy.sdss import quasars as xsq

from xastropy.xutils import xdebug as xdb

# Local
#sys.path.append(os.path.abspath("../Analysis/py"))
#import lls_sample as lls_s


####
#  Mean SDSS spectrum at z~3.6
def mean_sdss_spec(outfil=None):

    # Init
    mean_spec_fil = '/u/xavier/paper/LLS/taueff/Analysis/stack_DR7_z3.59_z3.63.fits'
    mean_hdu = fits.open(mean_spec_fil)
    head = mean_hdu[0].header
    zqso = head['AVG_Z']
    mean_wrest = mean_hdu[2].data*u.AA
    mean_flux = mean_hdu[0].data

    # Single QSO
    dr7 = xsq.SdssQuasars()
    qso = dr7[(2161,206)]
    qso.load_spec()

    # Plot boundaries
    print('zqso = {:g}'.format(zqso))
    #wvmnx = np.array((850., 1205.))*u.AA*(1+zqso)
    wvmnx = (3400.,5700)*u.AA
    mean_ymnx = (0.0, 2.0)

    if outfil == None:
        outfil='fig_mean_sdss_spec.pdf'

    elbl_dict = {r'Ly$\alpha$ (QSO)':[5450.,1.5], 
        r'Ly$\beta$ (QSO)':[4750.,1.6]}
    albl_dict = {r'Ly$\alpha$':[5200.,0.6], r'Ly$\beta$':[4600.,0.6],
        r'Ly$\gamma$':[4420.,0.50], r'Ly$\delta$':[4350.,0.4], 
        'LL':[912.*(1+zqso),0.25]}
    # Start the plot
    if outfil != None:
        pp = PdfPages(outfil)

    # Make the plots
    for ss in range(2):

        if ss ==0:
            ymnx = [0.,25.]
            wave = qso.spec.dispersion
            flux = qso.spec.flux.value
            title='Individual QSO at z=3.6'
        elif ss == 1:
            ymnx = mean_ymnx
            wave = mean_wrest*(1+zqso)
            flux = mean_flux
            title='150 QSOs at z=3.6'

        plt.figure(figsize=(8, 4))
        plt.clf()
        gs = gridspec.GridSpec(1, 1)

        # Axes
        ax = plt.subplot(gs[0,0])
        #ax.xaxis.set_minor_locator(plt.MultipleLocator(0.5))
        #ax.xaxis.set_major_locator(plt.MultipleLocator(1.))
        #ax.get_xaxis().get_major_formatter().set_useOffset(False)
        ax.set_xlim(wvmnx.value)
        ax.set_ylim(ymnx)
        ax.set_ylabel('Relative Flux')
        ax.set_xlabel(r'Observed Wavelength ($\AA$)')

        # Zero line
        ax.plot(wave,flux,'b')

        # Label
        if ss==1:
            for key in albl_dict.keys():
                ax.text(albl_dict[key][0],albl_dict[key][1],
                    key, color='green', ha='center')
            for key in elbl_dict.keys():
                ax.text(elbl_dict[key][0],elbl_dict[key][1],
                    key, color='red', ha='center')
        ax.text(3550.,ymnx[1]*0.85, title, color='black',ha='left',size=17.)


        # Data
        '''
        if ss > 0:
            #ax.fill_between( vmodel.dispersion, vmodel.flux, [1.]*len(wave),
            #                        color='blue', alpha=0.3)
            csz = 20.
            ax.text(1215.670, 0.5, 'EW', color='black', ha='center',
                size=csz)
            ax.text(1214.8, 0.2, 
                r'$W_\lambda = \, $'+'{:.2f}'.format(EW.value)+r' $\AA$',
                color='blue', ha='left', size=csz)
        '''

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
        mean_sdss_spec(outfil='mean_sdss_spec.pdf')


# Command line execution
if __name__ == '__main__':

    if len(sys.argv) == 1: # Figs
        flg_fig = 0 
        flg_fig += 2**0 # mean SDSS spec
        #flg_fig += 2**1 # COG
        #flg_fig += 2**1 # pLLS
        #flg_fig += 2**2 # Ambig
        #flg_fig += 2**3 # Summary
        #flg_fig += 2**4 # Check NHI
    else:
        flg_fig = sys.argv[1]

    main(flg_fig)
