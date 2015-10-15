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
from astropy import constants as const

from linetools.lists import linelist as lll

from xastropy.plotting import simple as xpsimp
from xastropy.cgm.cos_halos import COSHalos
from xastropy.spec.lines_utils import AbsLine
from xastropy.obs import finder
from xastropy.plotting import utils as xputils
from xastropy.xutils import xdebug as xdb

# Local
#sys.path.append(os.path.abspath("../Analysis/py"))
#import lls_sample as lls_s


####
#  CGM Kinematic Measures
def fig_kin_measures(outfil=None):

    # Init COS-Halos sightline
    cos_halos = COSHalos()
    cos_halos.load_single( ('J1220+3853', '225_38') )
    cgm_abs = cos_halos.cgm_abs[0]
    #print('zem = {:g}'.format(cgm_abs.abs_sys.zem))

    HI = 972.5368*u.AA

    # ########################################
    # Finder (out of order to avoid PDF issues)
    #finder.main([cgm_abs.name, cgm_abs.galaxy.coord], imsize=2.*u.arcmin, show_circ=False)

    # Start the plot
    if outfil is None:
        outfil='fig_kin_measures.pdf'
    pp = PdfPages(outfil)

    # dv
    lclr = 'cyan'
    vmnx = np.array([-350.,200.])*u.km/u.s
    ymnx = (-0.05, 1.4)
    for ss in range(2):
        plt.figure(figsize=(5, 5))
        plt.clf()
        gs = gridspec.GridSpec(1, 1)
        xpsimp.dark_bkgd(plt)

        # Axes
        spec = cos_halos.load_bg_cos_spec(0, HI)

        ax = plt.subplot(gs[0,0])
        #ax.xaxis.set_minor_locator(plt.MultipleLocator(0.5))
        #ax.xaxis.set_major_locator(plt.MultipleLocator(1.))
        #ax.get_xaxis().get_major_formatter().set_useOffset(False)
        ax.set_xlim(vmnx.value)
        ax.set_ylim(ymnx)
        ax.set_xlabel('Relative Velocity (km/s)')
        ax.set_ylabel('Normalized Flux')

        # Zero line
        ax.plot(vmnx.value, (0.,0.), ':', color='lightgreen')
        ax.plot((0.,0.), ymnx, '--', color='lightgray')

        # Data
        velo = spec.relative_vel(HI*(cgm_abs.galaxy.z+1))
        ax.plot(velo, spec.flux, color='white', drawstyle='steps')

        # Label
        if ss == 0:
            ax.plot((-95,-95), ymnx, '--', color='cyan')
            ax.text(-50., 1.2, r'$\delta v$', ha='center', fontsize=21., color=lclr)
        elif ss == 1:
            ax.plot((-170,-170), ymnx, '-', color='cyan')
            ax.plot((60,60), ymnx, '-', color='cyan')
            ax.text(-50., 1.2, r'$\Delta v$', ha='center', fontsize=21., color=lclr)
        # Fonts
        xputils.set_fontsize(ax,17.)

        # Write
        plt.tight_layout(pad=0.2,h_pad=0.,w_pad=0.1)
        pp.savefig()
        plt.close()


    # ########################################
    # Simple HI stack plot
    trans = np.array([HI.value, 1031.9261])*u.AA
    lbls = [r'HI Ly$\gamma$', 'OVI 1031']

    plt.figure(figsize=(8, 5))
    plt.clf()
    gs = gridspec.GridSpec(2, 1)

    for qq,itran,lbl in zip(range(len(trans)),trans,lbls):
        # Load
        spec = cos_halos.load_bg_cos_spec(0, itran)
        # Axes
        ax = plt.subplot(gs[qq,0])
        ax.set_xlim(-500., 500.)
        ax.set_ylim(-0.1, 1.3)
        ax.set_ylabel('Normalized Flux')
        if qq == 0:
            ax.xaxis.set_ticklabels([])
        else:
            ax.set_xlabel('Relative Velocity (km/s)')
        # Velo
        velo = spec.relative_vel(itran*(cgm_abs.galaxy.z+1))
        ax.plot(velo, spec.flux, color='white', drawstyle='steps', lw=1.3)
        # Label
        ax.text(-400., 0.2, lbl, ha='left', fontsize=21., color=lclr)
        # Fonts
        xputils.set_fontsize(ax,17.)
    # Write
    plt.tight_layout(pad=0.2,h_pad=0.,w_pad=0.1)
    pp.savefig()
    plt.close()


    # Finish
    print('tlk_coshalos: Wrote {:s}'.format(outfil))
    pp.close()



####
#  CGM dv/Dv Histograms for COS-Halos
def fig_coshalo_hist(outfil=None, cos_halos=None):

    # Init COS-Halos sightline
    if cos_halos is None:
        cos_halos = COSHalos()
        cos_halos.load_mega(skip_ions=True)
        cos_halos.load_abskin(flg=1)

    # ########################################
    # Finder (out of order to avoid PDF issues)
    #finder.main([cgm_abs.name, cgm_abs.galaxy.coord], imsize=2.*u.arcmin, show_circ=False)

    # Start the plot
    if outfil is None:
        outfil='fig_coshalo_kin_hist.pdf'
    pp = PdfPages(outfil)

    # dv
    lclr = 'cyan'
    plt.figure(figsize=(8, 4))
    plt.clf()
    xpsimp.dark_bkgd(plt)
    gs = gridspec.GridSpec(1, 2)

    for ss in range(2):
        ax = plt.subplot(gs[ss])
        #ax.xaxis.set_minor_locator(plt.MultipleLocator(0.5))
        #ax.get_xaxis().get_major_formatter().set_useOffset(False)
        if ss == 0:
            xmnx = [-300., 300.] # km/s
            ymnx = (0, 8)
            ax.set_xlabel(r'$\delta v$ (km/s)')
            binsz = 25.
            gdk = cos_halos.abs_kin('Metal')['flg'] > 0
            val = cos_halos.abs_kin('Metal')['delta_v'][gdk]
            ax.text(-200., 7, 'Metals', ha='left', fontsize=21., color=lclr)
        elif ss == 1:
            xmnx = [0., 400.] # km/s
            ymnx = (0, 10)
            ax.set_xlabel(r'$\Delta v$ (km/s)')
            val = cos_halos.abs_kin('Metal')['Dv'][gdk]
        ax.set_ylabel('N')
        ax.set_xlim(xmnx)
        ax.set_ylim(ymnx)
        ax.minorticks_on()
        ax.xaxis.set_major_locator(plt.MultipleLocator(100.))

        # Zero line
        #ax.plot(vmnx.value, (0.,0.), ':', color='lightgreen')
        #ax.plot((0.,0.), ymnx, '--', color='lightgray')

        # Histogram
        bins = np.arange(xmnx[0], xmnx[1] + binsz, binsz)
        ax.hist(val, bins=bins, color='cyan', histtype='stepfilled', 
            edgecolor='white')

        # Labels
        #ax.text(-50., 1.2, r'$\Delta v$', ha='center', fontsize=21., color=lclr)
        # Fonts
        xputils.set_fontsize(ax,17.)

    # Write
    plt.tight_layout(pad=0.2,h_pad=0.,w_pad=0.1)
    pp.savefig()
    plt.close()

    # Finish
    print('tlk_coshalos: Wrote {:s}'.format(outfil))
    pp.close()


####
#  COS-Halos Kin correlations
def fig_coshalo_corr(outfil=None, cos_halos=None):

    # Init COS-Halos sightline
    if cos_halos is None:
        cos_halos = COSHalos()
        cos_halos.load_mega(skip_ions=True)
        cos_halos.load_abskin(flg=1)

    # ########################################
    # Finder (out of order to avoid PDF issues)
    #finder.main([cgm_abs.name, cgm_abs.galaxy.coord], imsize=2.*u.arcmin, show_circ=False)

    # Start the plot
    if outfil is None:
        outfil='fig_coshalo_kin_corr.pdf'
    pp = PdfPages(outfil)

    # dv
    lclr = 'cyan'
    plt.figure(figsize=(8, 4))
    plt.clf()
    xpsimp.dark_bkgd(plt)
    gs = gridspec.GridSpec(2, 3)

    dv_ymnx = [-300., 290.] # km/s
    Dv_ymnx = [0., 410.] # km/s
    gdk = cos_halos.abs_kin('Metal')['flg'] > 0

    for ss in range(6):
        ax = plt.subplot(gs[ss%2,ss/2])
        #ax.xaxis.set_minor_locator(plt.MultipleLocator(0.5))
        #ax.get_xaxis().get_major_formatter().set_useOffset(False)
        if (ss/2) == 0: # dv vs. rho
            xmnx = (0, 150.)
            xlbl = r'$\rho$ (kpc)'
            xval = cos_halos.rho.value[gdk]
        elif (ss/2) == 1:
            xmnx = (9, 12.)
            xlbl = r'$\log M*$'
            xval = cos_halos.stellar_mass[gdk]
        elif (ss/2) == 2:
            xmnx = (-13., -8)
            xlbl = 'sSFR'
            sfr = [float(sfr[1]) for sfr in cos_halos.sfr]
            xval = (np.log10(sfr)-cos_halos.stellar_mass)[gdk]
        #
        ax.set_xlim(xmnx)
        if ss == 0:
            ax.set_ylabel(r'$\delta v$ (km/s)')
        elif ss == 1:
            ax.set_ylabel(r'$\Delta v$ (km/s)')

        if (ss%2) == 0:
            ax.xaxis.set_ticklabels([])
            ymnx = dv_ymnx
            yval = cos_halos.abs_kin('Metal')['delta_v'][gdk]
        else:
            ymnx = Dv_ymnx
            yval = cos_halos.abs_kin('Metal')['Dv'][gdk]
            ax.set_xlabel(xlbl)
        ax.set_ylim(ymnx)
        ax.minorticks_on()
        #ax.xaxis.set_major_locator(plt.MultipleLocator(100.))

        # Zero line
        #ax.plot(vmnx.value, (0.,0.), ':', color='lightgreen')
        #ax.plot((0.,0.), ymnx, '--', color='lightgray')

        # Scatter plot
        ax.scatter(xval, yval, marker='o', color='cyan')

        # Labels
        #ax.text(-50., 1.2, r'$\Delta v$', ha='center', fontsize=21., color=lclr)
        # Fonts
        xputils.set_fontsize(ax,13.)

    # Write
    plt.tight_layout(pad=0.2,h_pad=0.,w_pad=0.1)
    pp.savefig()
    plt.close()

    # Finish
    print('tlk_coshalos: Wrote {:s}'.format(outfil))
    pp.close()

####
#  CGM dv/Dv Histograms for COS-Dwarfs
def fig_cosdwarf_hist(outfil=None, cos_dwarfs=None):

    # Init COS-Halos sightline
    if cos_dwarfs is None:
        cos_dwarfs = COSDwarfs()
        cos_dwarfs.load_mega(skip_ions=True)
        cos_dwarfs.load_abskin(flg=1)

    # ########################################
    # Finder (out of order to avoid PDF issues)
    #finder.main([cgm_abs.name, cgm_abs.galaxy.coord], imsize=2.*u.arcmin, show_circ=False)

    # Start the plot
    if outfil is None:
        outfil='fig_cosdwarf_kin_hist.pdf'
    pp = PdfPages(outfil)

    # dv
    lclr = 'yellow'
    plt.figure(figsize=(8, 4))
    plt.clf()
    xpsimp.dark_bkgd(plt)
    gs = gridspec.GridSpec(1, 2)

    for ss in range(2):
        ax = plt.subplot(gs[ss])
        #ax.xaxis.set_minor_locator(plt.MultipleLocator(0.5))
        #ax.get_xaxis().get_major_formatter().set_useOffset(False)
        if ss == 0:
            xmnx = [-400., 400.] # km/s
            ymnx = (0, 10)
            ax.set_xlabel(r'$\delta v$ (km/s)')
            binsz = 25.
            gdk = cos_dwarfs.abs_kin('HI')['flg'] > 0
            val = cos_dwarfs.abs_kin('HI')['delta_v'][gdk]
            ax.text(-200., 8, 'HI', ha='left', fontsize=21., color=lclr)
            ax.xaxis.set_major_locator(plt.MultipleLocator(200.))
        elif ss == 1:
            xmnx = [0., 400.] # km/s
            ymnx = (0, 10)
            ax.set_xlabel(r'$\Delta v$ (km/s)')
            val = cos_dwarfs.abs_kin('HI')['Dv'][gdk]
            ax.xaxis.set_major_locator(plt.MultipleLocator(100.))
        ax.set_ylabel('N')
        ax.set_xlim(xmnx)
        ax.set_ylim(ymnx)
        ax.minorticks_on()

        # Zero line
        #ax.plot(vmnx.value, (0.,0.), ':', color='lightgreen')
        #ax.plot((0.,0.), ymnx, '--', color='lightgray')

        # Histogram
        bins = np.arange(xmnx[0], xmnx[1] + binsz, binsz)
        ax.hist(val, bins=bins, color='yellow', histtype='stepfilled', 
            edgecolor='white')

        # Labels
        #ax.text(-50., 1.2, r'$\Delta v$', ha='center', fontsize=21., color=lclr)
        # Fonts
        xputils.set_fontsize(ax,17.)

    # Write
    plt.tight_layout(pad=0.2,h_pad=0.,w_pad=0.1)
    pp.savefig()
    plt.close()

    # Finish
    print('tlk_cgm_kin: Wrote {:s}'.format(outfil))
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

    # Experimental design
    if (flg_fig % 2**1) >= 2**0:
        fig_kin_measures()


# Command line execution
if __name__ == '__main__':

    if len(sys.argv) == 1: # Figs
        flg_fig = 0 
        flg_fig += 2**0 # Kinematic Measures
        #flg_fig += 2**1 # Image gallery
        #flg_fig += 2**2 # Abs gallery
    else:
        flg_fig = sys.argv[1]

    main(flg_fig)
