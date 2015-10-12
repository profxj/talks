# Module for pLLS plots.  Mainly XQ-100 

from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import glob, copy, os, sys

from scipy.interpolate import interp1d
from scipy.interpolate import splev, splrep # bspline


import matplotlib as mpl
mpl.rcParams['font.family'] = 'stixgeneral'
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import json
from xastropy.spec import continuum as xspc
from xastropy.plotting import utils as xputils


from astropy import units as u
from astropy.units import Unit, Quantity
from astropy.io import ascii, fits
from astropy.table import QTable, Table, Column
from astropy.coordinates import SkyCoord
from astropy import constants as const

from linetools.spectra.xspectrum1d import XSpectrum1D

from xastropy.igm.abs_sys.lls_utils import LLSSystem
from xastropy.stats import basic as xsb
from xastropy.xutils import fits as xxf
from xastropy.xutils import xdebug as xdb
from xastropy.atomic import ionization as xai
from xastropy.xguis import utils as xxgu
from xastropy.spec import voigt as xsv

# Local 
#sys.path.append(os.path.abspath("./py"))
#import qpq_spec as qpqs

##
def load_files():
    '''Load up list of spectra and summary Table'''
    # Spectra
    #all_spec = glob.glob(os.getenv('DROPBOX_DIR')+'/XQ-100/data/*_uvb_*flx.fits')
    all_spec = glob.glob(os.getenv('DROPBOX_DIR')+'/XQ-100/data/*_uvb.fits')

    # Summary Table
    summary = ascii.read(os.getenv('DROPBOX_DIR')+'/XQ-100/data/XQ-100_summary.ascii', 
        format='fixed_width_no_header', comment='#', #data_start=1, 
        names=('Name', 'NED', 'RAS', 'DECS', 'zQSO', 'R_APM'), 
        col_starts=(0,13,42,56,70,77), col_ends=(9,38,53,67,77,84))

    # Return
    return all_spec, summary

def load_spectrum(spec_fil):
    '''Load X-Shooter spectra'''
    # Load Spectrum
    uvb_spec = XSpectrum1D.from_file(spec_fil)
    vis_specfil = spec_fil.replace('uvb','vis')
    vis_spec = XSpectrum1D.from_file(vis_specfil)
    comb_spec = uvb_spec.splice(vis_spec)
    comb_spec.filename = spec_fil
    # Return
    return comb_spec

#### ########################## #########################
def xq100_example(outfil='fig_xq100_example.pdf'):

    # Load spectrum
    spec_fil = os.getenv('DROPBOX_DIR')+'/XQ-100/data/J0030-5159_uvb.fits'
    spec = load_spectrum(spec_fil)

    # Start the plot
    if outfil is not None: 
        pp = PdfPages(outfil)

    plt.figure(figsize=(8, 5))
    plt.clf()
    gs = gridspec.GridSpec(1,1)
    fsz = 17.

    # Read and plot full spectrum
    ax_full = plt.subplot(gs[0])

    # Limits
    wvmnx = [3500., 7000.]
    gdpix = np.where(spec.dispersion<wvmnx[1]*u.AA)[0]
    perc = xsb.perc(spec.flux[gdpix])
    ax_full.set_xlim(wvmnx)
    #ax_full.set_ylim(-0.05*perc[1], 1.1*perc[1])
    ax_full.set_ylim(-1e-18, 7e-17)

    # Plot
    ax_full.plot(spec.dispersion, spec.flux, color='black', lw=1.0)
    ax_full.plot(wvmnx, [0.,0.], '--', color='green')

    # Label
    ax_full.set_xlabel('Wavelength')
    ax_full.set_ylabel('Relative Flux')
    ax_full.text(0.05, 0.9, 'XQ100 J0030-5159', transform=ax_full.transAxes, 
        color='black',
        size='x-large', ha='left', va='center', bbox={'facecolor':'white'})

    # Font size
    xputils.set_fontsize(ax_full,fsz)

    # Finish page
    plt.tight_layout(pad=0.2,h_pad=0.3,w_pad=0.0)
    pp.savefig()
    plt.close()
    print('Wrote: {:s}'.format(outfil))
    # Finish 
    pp.close()

#### ########################## #########################
def xq100_plls_ex(outfil='fig_xq100_plls_ex.pdf'):

    # Load spectrum
    spec_fil = os.getenv('DROPBOX_DIR')+'/XQ-100/data/J0030-5159_uvb.fits'
    spec = load_spectrum(spec_fil)

    # Generate model
    model_file = os.getenv('DROPBOX_DIR')+'/XQ-100/LLS/convg_J0030-5159_llsfit.json'
    with open(model_file) as data_file:    
        lls_dict = json.load(data_file)
    # Continuum
    basec = Quantity(lls_dict['conti'])
    telfer_spec = XSpectrum1D.from_tuple((spec.dispersion.value, 
        basec.value))
 

    # Start the plot
    if outfil is not None: 
        pp = PdfPages(outfil)

    plt.figure(figsize=(8, 4))
    plt.clf()
    gs = gridspec.GridSpec(1,1)
    fsz = 17.

    # Read and plot full spectrum
    ax_full = plt.subplot(gs[0])

    # Limits
    wvmnx = [3600., 5000.]
    gdpix = np.where(spec.dispersion<wvmnx[1]*u.AA)[0]
    perc = xsb.perc(spec.flux[gdpix])
    ax_full.set_xlim(wvmnx)
    #ax_full.set_ylim(-0.05*perc[1], 1.1*perc[1])
    ax_full.set_ylim(-1e-18, 4e-17)

    # Plot
    ax_full.plot(spec.dispersion, spec.flux, color='black', lw=1.0)
    ax_full.plot(wvmnx, [0.,0.], '--', color='green')

    # Label
    ax_full.set_xlabel('Wavelength')
    ax_full.set_ylabel('Relative Flux')
    ax_full.text(0.05, 0.9, 'XQ100 J0030-5159', transform=ax_full.transAxes, 
        color='black',
        size='x-large', ha='left', va='center', bbox={'facecolor':'white'})

    all_ills = []
    all_lls = []
    all_zlls = []
    for key in lls_dict['LLS'].keys():
        new_sys = LLSSystem(NHI=lls_dict['LLS'][key]['NHI'])
        new_sys.zabs = lls_dict['LLS'][key]['z']
        new_sys.fill_lls_lines(bval=lls_dict['LLS'][key]['bval']*u.km/u.s)
        #
        all_ills.append(new_sys)
        all_lls.append(new_sys)
        all_zlls.append(new_sys.zabs)
    # Model
    norm_flux = lls_model(spec.dispersion, all_ills, smooth=lls_dict['smooth'])
    continuum = telfer_spec.flux * lls_dict['conti_model']['Norm'] * (
        spec.dispersion.value/
        lls_dict['conti_model']['piv_wv'])**lls_dict['conti_model']['tilt']

    #  Model
    ax_full.plot(spec.dispersion, continuum*norm_flux)

    # Font size
    xputils.set_fontsize(ax_full,fsz)

    # Finish page
    plt.tight_layout(pad=0.2,h_pad=0.3,w_pad=0.0)
    pp.savefig()
    plt.close()
    print('Wrote: {:s}'.format(outfil))
    # Finish 
    pp.close()



def lls_model(wave, all_lls, smooth=0.):
    '''Generate an absorption model '''
    from linetools.spectra import convolve as lsc
    from xastropy.igm.abs_sys import lls_utils as xialu

    # Tau from LLS
    all_tau_model = xialu.tau_multi_lls(wave,all_lls)
    # Flux and smooth
    flux = np.exp(-1. * all_tau_model)
    if smooth > 0:
        lls_model = lsc.convolve_psf(flux, smooth)
    else:
        lls_model = flux

    # Finish
    norm_flux = lls_model 
    # Return
    return norm_flux

    
# ##################################################
# ##################################################
# ##################################################
# Command line execution for testing
# ##################################################
if __name__ == '__main__':


    if len(sys.argv) == 1: # TESTING

        flg_fig = 0 
        #flg_fig += 2**0  # XQ-100 example
        flg_fig += 2**1  # XQ-100 pLLS
        #flg_fig += 2**1  # XAbsID
        #flg_fig += 2**2  # XVelPlt Gui
        #flg_fig += 2**3  # XVelPlt Gui without ID list; Also tests select wave
        #flg_fig += 2**4  # XAODM Gui
        #flg_fig += 2**5  # Fit LLS GUI

        if (flg_fig % 2**1) >= 2**0:
            xq100_example()

        if (flg_fig % 2**2) >= 2**1:
            xq100_plls_ex()
