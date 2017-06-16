# Module for the Spectra figures for talks
# Imports
from __future__ import print_function

import numpy as np
import glob, os, sys
import pdb

import matplotlib as mpl
mpl.rcParams['font.family'] = 'stixgeneral'
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec

from astropy.io import ascii
from astropy import units as u
from astropy import constants as const

'''
from xastropy.abund import solar as xsolar
from xastropy.atomic.elements import ELEMENTS
from xastropy.atomic import ionization as xion
from xastropy.cgm.cos_halos import COSHalos

from xastropy.igm.cuba import CUBA
from xastropy import spec as xsp
from xastropy.plotting import simple as xplots
from xastropy.plotting import utils as xputils
from xastropy.xutils import xdebug as xdb
'''


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

####
#  Series of plots illustrating COS-Halos experiment
def fig_ionstate(outfil=None):

    # Init COS-Halos sightline
    cos_halos = COSHalos()
    cos_halos.load_single( ('J1016+4706', '274_6') )
    cgm_abs = cos_halos.cgm_abs[0]
    FeH = -1.

    # CUBA
    cuba = CUBA()
    npt = 30
    xeV = np.linspace(13.6, 130., npt)*u.eV
    phi = np.zeros(npt)
    for kk,ixeV in enumerate(xeV):
        phi[kk] = cuba.phi(0.2, min_energy=ixeV).value

    # Ions
    Zions = [ (1,1), # HI
        (6,2),  # CII
        (6,3),  # CIII
        (7,2),  # NII
        (7,3),  # NIII
        (7,5),  # NV
        (8,6),  # OVI
        (12,1), # MgI
        (12,2), # MgII
        (14,2), # SiII
        (14,3), # SiIII
        (14,4), # SiIV
    ]

    # Start the plot
    if outfil is None:
        outfil='fig_ionstate.pdf'
    pp = PdfPages(outfil)

    # Lya spec
    lclr = 'red'
    pclr = 'blue'
    xmnx = (0., 130.) # eV
    ymnx = (-1.5, 2.)


    for ss in range(5):
        plt.figure(figsize=(7, 5))
        plt.clf()
        gs = gridspec.GridSpec(1, 1)

        ax = plt.subplot(gs[0,0])
        #ax.xaxis.set_minor_locator(plt.MultipleLocator(0.5))
        ax.yaxis.set_major_locator(plt.MultipleLocator(1.))
        #ax.get_xaxis().get_major_formatter().set_useOffset(False)
        ax.set_xlim(xmnx)
        ax.set_ylim(ymnx)
        ax.set_xlabel('Ionization Potential (eV)')
        ax.set_ylabel(r'$\log N_{\rm ion} - \log N_{\rm HI} + 12 - \epsilon_\odot$ - [Fe/H]')
        ax.minorticks_on()

        # Zero line
        ax.plot(xmnx, (0.,0.), '--', color='gray')

        # Data
        HIclm = cgm_abs.abs_sys.ions[(1,1)]['CLM']
        for Zion in Zions:
            # Column and flag
            try:
                clm = cgm_abs.abs_sys.ions[Zion]['CLM']
            except KeyError:
                print('Zion = {:d},{:d} not analyzed'.format(Zion[0],Zion[1]))
                continue
            else:
                flg_clm = cgm_abs.abs_sys.ions[Zion]['FLG_CLM']
                sig_clm = cgm_abs.abs_sys.ions[Zion]['SIG_CLM']

            # IP
            elm = ELEMENTS[Zion[0]]
            IP = elm.ionenergy[Zion[1]-1]

            # Solar abund
            abund = xsolar.abund(Zion[0])

            yval = clm - HIclm + 12 - abund 
            if Zion != (1,1):
                yval = yval + FeH

            # Plot
            if flg_clm==1: 
                mark = 'o'
                pclr = 'blue'
                yoff = 0.1
            elif flg_clm==2: 
                mark = (3,0,0)
                pclr = 'red'
                yoff = -0.3
            elif flg_clm==3: 
                mark = (3,0,180)
                pclr = 'gray'
                yoff = 0.1
            ax.scatter(IP, yval, marker=mark, color=pclr, s=35.)
            # Label
            ion_name = xion.ion_name(Zion)
            ax.text(IP, yval+yoff, ion_name, color=pclr, fontsize=17., ha='center')

        #ax.text(1480., 0.4, cgm_abs.field, ha='left', fontsize=21., color=lclr)

        if ss>=1: # Add EUVB
            ax_euvb = ax.twinx()
            y0 = np.log10(phi[0])
            ax_euvb.plot(xeV,np.log10(phi)-y0, 'k-') # Normalize 1 Ryd to 1
            ax_euvb.set_ylabel(r'$\log \; \Phi_{\rm EUVB} ({\rm E > IP})$')
            xputils.set_fontsize(ax_euvb,17.)
            ax_euvb.set_xlim(xmnx)
            ax_euvb.set_ylim(-1.2,0.)
            ax_euvb.text(100., -1., 'EUVB', color='black', size=17.)
        if ss>=2: # Add shading
            xrng = np.array([10.,55])
            ax.fill_between( xrng, [ymnx[0]]*2, [ymnx[1]]*2, color='blue', alpha=0.3)
            ax.text(np.mean(xrng), -1.2, 'Photoionization\n'+r'$U \equiv \Phi/cn_{\rm H}$',
                color='black', size=19., ha='center')
        if ss>=3: # Add shading
            xrng = np.array([70., 120.])
            ax.fill_between( xrng, [ymnx[0]]*2, [ymnx[1]]*2, color='red', alpha=0.3)
            ax.text(np.mean(xrng), -1.2, 'Another phase', color='black', size=19., ha='center')

        if ss>=4: # EUV lines
            Zions2 = [ (10,8), (16,4), (16,5), 
                (8,2), (8,3), (8,5)]
            for Zion2 in Zions2:
                elm = ELEMENTS[Zion2[0]]
                IP = elm.ionenergy[Zion2[1]-1]
                yval = 0.2
                yoff = 0.1
                ax.scatter(IP, yval, marker='o', color='black', s=35.)
                # Label
                try:
                    ion_name = xion.ion_name(Zion2)
                except TypeError:
                    ion_name = 'NeVIII'
                ax.text(IP, yval+yoff, ion_name, color='black', fontsize=17., ha='center')
            #ax.set_xlim(0.,200)

        # Fonts
        xputils.set_fontsize(ax,17.)

        # Write
        plt.tight_layout(pad=0.2,h_pad=0.,w_pad=0.1)
        pp.savefig()
        plt.close()


    # Finish
    print('tlk_coshalos: Wrote {:s}'.format(outfil))
    pp.close()

def fig_pie_chart():
    """ Plot MCMC outputs for ne/nH from Neeleman+15
    """
    # Read COS-Halos

    lsz = 13.
    tsz = 18.

    # Plot
    fig = plt.figure(figsize=(8, 5))#, dpi=700)
    gs = gridspec.GridSpec(2,2)

    # Clusters
    ax_c_mass = plt.subplot(gs[0])

    labels = 'Galaxies', 'ICL', 'ICM'
    c_clrs = ['green', 'orange', 'red']
    stars = 0.0003 / (3e-4 + 0.00155)
    icl = 0.15
    icm = 1. - icl - stars
    wedges = [stars, icl, icm]
    explode = (0, 0, 0.1)

    _, texts, _ = ax_c_mass.pie(wedges, explode=explode, labels=labels, autopct='%1.0f%%',
        startangle=90, colors=c_clrs)
    for text in texts:
        text.set_fontsize(lsz)
    ax_c_mass.axis('equal')
    ax_c_mass.set_xlabel('Cluster Baryonic Mass', size=tsz)

    # Cluster metals
    ax_c_metals = plt.subplot(gs[2])
    Z_comp = (1., 0.5, 0.3) # galaxies, ICL, ICM
    Z_tot = np.sum(np.array(wedges)*Z_comp)
    c_Z_frac = np.array(wedges)*Z_comp / Z_tot
    _, texts, _ = ax_c_metals.pie(c_Z_frac, explode=explode, labels=labels, autopct='%1.0f%%',
                  startangle=60, colors=c_clrs)
    for text in texts:
        text.set_fontsize(lsz)
    ax_c_metals.axis('equal')
    ax_c_metals.set_xlabel('Cluster Metal Mass', size=tsz)

    # Galaxy mass
    ax_g_mass = plt.subplot(gs[1])

    labels = 'Galaxy', r'Cool CGM$^1$', r'Hot CGM$^2$'
    g_clrs = ['green', 'skyblue', 'red']
    mhalo_b = 2e11 # P17
    galaxy = 0.17   # check
    cool = 9e10 / mhalo_b
    warm = 1. - galaxy - cool
    wedges = [galaxy, cool, warm]
    explode = (0, 0, 0)

    _, texts, _ = ax_g_mass.pie(wedges, explode=explode, labels=labels, autopct='%1.0f%%',
                  startangle=40, colors=g_clrs)
    for text in texts:
        text.set_fontsize(lsz)
    ax_g_mass.axis('equal')
    ax_g_mass.set_xlabel('Galaxy Baryonic Mass', size=tsz)


    # Galaxy metals
    ax_g_metals = plt.subplot(gs[3])
    Z_comp = (1., 0.5, 0.5) # galaxy, cold, warm
    Z_tot = np.sum(np.array(wedges)*Z_comp)
    g_Z_frac = np.array(wedges)*Z_comp / Z_tot
    _, texts, _ = ax_g_metals.pie(g_Z_frac, explode=explode, labels=labels, autopct='%1.0f%%',
                    startangle=30, colors=g_clrs)
    for text in texts:
        text.set_fontsize(lsz)
    ax_g_metals.axis('equal')
    ax_g_metals.set_xlabel('Galaxy Metal Mass', size=tsz)

    # Write
    plt.tight_layout(pad=0.2,h_pad=0.,w_pad=0.1)

    outfil='fig_pie_chart.png'
    plt.savefig(outfil, dpi=800)
    plt.close()
    print('Generated {:s}'.format(outfil))

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
        fig_cog_abs_web()  # Did I accidentally remove this??

    # Ionization
    if (flg_fig % 2**3) >= 2**2:
        fig_ionstate()

    # Pie chart (with ICM)
    if flg_fig & (2**3):
        fig_pie_chart()


    # Command line execution
if __name__ == '__main__':

    if len(sys.argv) == 1: # Figs
        flg_fig = 0 
        #flg_fig += 2**0 # EW
        #flg_fig += 2**1 # COG
        #flg_fig += 2**2 # Ionization state
        flg_fig += 2**3 # Check NHI

    else:
        flg_fig = sys.argv[1]

    main(flg_fig)
