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
#  Series of plots illustrating COS-Halos experiment
def fig_experiment(outfil=None):

    # Init COS-Halos sightline
    cos_halos = COSHalos()
    cos_halos.load_single( ('J0950+4831','177_27'))
    cgm_abs = cos_halos.cgm_abs[0]
    print('zem = {:g}'.format(cgm_abs.abs_sys.zem))

    # ########################################
    # Finder (out of order to avoid PDF issues)
    #finder.main([cgm_abs.name, cgm_abs.galaxy.coord], imsize=2.*u.arcmin, show_circ=False)

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
            spec = cos_halos.load_bg_cos_spec(0, 1215.6700*u.AA)
            wvmnx = np.array([1455., 1490.])*u.AA
            scl_wv = 1.
        elif ss == 1:
            scl_wv = (1+cgm_abs.galaxy.z)
        elif ss == 2:
            spec = cos_halos.load_bg_cos_spec(0, 1025.7222*u.AA)
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
    # Galaxy spectrum
    plt.figure(figsize=(8, 5))
    plt.clf()
    gs = gridspec.GridSpec(1, 1)

    # Load
    spec = cos_halos.load_gal_spec(0)
    # Axes
    ax = plt.subplot(gs[0,0])
    ax.set_xlim(4000., 8000.)
    ax.set_ylim(-0.1, 25.)
    ax.set_xlabel(r'Wavelength ($\AA$)')
    ax.set_ylabel('Relative Flux')
    # Velo
    ax.plot(spec.dispersion, spec.flux, 'k')#, lw=1.3)
    ax.plot(spec.dispersion, spec.sig, 'r:')#, lw=1.3)
    # Label
    #ax.text(-800., 0.2, lbl, ha='left', fontsize=21., color=lclr)
    # Fonts
    xputils.set_fontsize(ax,17.)
    # Write
    plt.tight_layout(pad=0.2,h_pad=0.,w_pad=0.1)
    pp.savefig()
    plt.close()

    # ########################################
    # Simple HI stack plot
    trans = np.array([1215.6701,  1025.7223])*u.AA
    lbls = [r'HI Ly$\alpha$', r'HI Ly$\beta$']

    plt.figure(figsize=(8, 5))
    plt.clf()
    gs = gridspec.GridSpec(2, 1)

    for qq,itran,lbl in zip(range(len(trans)),trans,lbls):
        # Load
        spec = cos_halos.load_bg_cos_spec(0, itran)
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

    # ########################################
    # Greater stack plot
    trans = np.array([1215.6701,  1025.7223, 
        2852.9642, 2803.5315,
        1083.9937, 1206.5, 
        1393.755, 1037.6167])*u.AA
    lbls = [r'HI Ly$\alpha$', r'HI Ly$\beta$',
        'MgI 2852', 'MgII 2796',
        'NII 1083', 'SiIII 1206',
        'SiIV 1393', 'OVI 1037']


    plt.figure(figsize=(8, 5))
    plt.clf()
    gs = gridspec.GridSpec(2, 4)
    xmnx=(-500., 500.)
    ymnx = (-0.1, 1.3)

    wbox= {'facecolor':'white', 'edgecolor':'white'}
    for qq,itran,lbl in zip(range(len(trans)),trans,lbls):
        # Load
        spec = cos_halos.load_bg_cos_spec(0, itran)
        # Axes
        ax = plt.subplot(gs[qq % 2,qq // 2])
        ax.set_xlim(xmnx) 
        ax.xaxis.set_major_locator(plt.MultipleLocator(300.))
        ax.set_ylim(ymnx)
        if qq % 2 == 0:
            ax.xaxis.set_ticklabels([])
        if qq < 2:
            ax.set_ylabel('Normalized Flux')
        else:
            ax.yaxis.set_ticklabels([])
        #if qq == 3:
        #    ax.set_xlabel('Relative Velocity (km/s)')
        # Velo
        velo = spec.relative_vel(itran*(cgm_abs.galaxy.z+1))
        ax.plot(velo, spec.flux, 'k', drawstyle='steps', lw=1.3)
        ax.plot([0.,0.], ymnx, 'g:')
        # Label
        ax.text(10., 1.15, lbl, ha='left', fontsize=15., color=lclr, bbox=wbox)
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
#  Series of plots illustrating COS-Halos experiment
def fig_images(outfil=None):

    from xastropy.obs import x_getsdssimg as xgs

    # Init COS-Halos sightline
    cos_halos = COSHalos()
    cos_halos.load_mega(skip_ions=True)
    #xdb.set_trace()

    # Start the plot
    if outfil is None:
        outfil='fig_images.pdf'
    pp = PdfPages(outfil)

    dx = 0.2
    dy = 0.55
    imsize = 1.2

    #for ss in range(5):
    for ss in range(1):
        fig = plt.figure(figsize=(8, 5))
        fig.clf()
        # Setup for dark
        xpsimp.dark_bkgd(fig)
        gs = gridspec.GridSpec(1, 1)

        # Axes
        ax = fig.subplot(gs[0,0])
        ax.set_xlim(9.4, 11.7) 
        #ax.xaxis.set_major_locator(plt.MultipleLocator(300.))
        ax.set_ylim(-13,-9)
        ax.set_ylabel('sSFR')
        ax.set_xlabel('log M*')

        # Plot
        if ss == 0:
            mstar = [cgm_abs.galaxy.stellar_mass for cgm_abs in cos_halos.cgm_abs]
            sSFR = [np.log10(cgm_abs.galaxy.sfr[1])-cgm_abs.galaxy.stellar_mass for cgm_abs in cos_halos.cgm_abs]
            ax.scatter(mstar,sSFR)
        elif ss == 1: # Show our first one
            cgm_abs = cos_halos[('J0950+4831','177_27')]
            img, oBW = xgs.getimg(cgm_abs.galaxy.coord.ra.deg, 
                cgm_abs.galaxy.coord.dec.deg, imsize=imsize)
            # Figure out placement
            ximg = cgm_abs.galaxy.stellar_mass
            yimg = np.log10(cgm_abs.galaxy.sfr[1])-cgm_abs.galaxy.stellar_mass
            # Show
            ax.imshow(img,extent=(ximg-dx, ximg+dx, yimg-dy, yimg+dy),aspect=dx/dy)
        elif ss == 2: # Show two
            cgm_list = [('J0950+4831','177_27'), ('J1245+3356', '236_36')]
            for icgm in cgm_list:
                cgm_abs = cos_halos[icgm]
                img, oBW = xgs.getimg(cgm_abs.galaxy.coord.ra.deg, 
                    cgm_abs.galaxy.coord.dec.deg, imsize=imsize)
                # Figure out placement
                ximg = cgm_abs.galaxy.stellar_mass
                yimg = np.log10(cgm_abs.galaxy.sfr[1])-cgm_abs.galaxy.stellar_mass
                # Show
                ax.imshow(img,extent=(ximg-dx, ximg+dx, yimg-dy, yimg+dy),aspect=dx/dy)
        elif ss >= 3: # Show half
            if ss == 3:
                iend = 20
            else:
                iend = -1
            for cgm_abs in cos_halos.cgm_abs[0:iend]:
                img, oBW = xgs.getimg(cgm_abs.galaxy.coord.ra.deg, 
                    cgm_abs.galaxy.coord.dec.deg, imsize=1.2)
                # Figure out placement
                ximg = cgm_abs.galaxy.stellar_mass
                yimg = np.log10(cgm_abs.galaxy.sfr[1])-cgm_abs.galaxy.stellar_mass
                # Show
                scl=2.
                ax.imshow(img,extent=(ximg-dx/scl, ximg+dx/scl, yimg-dy/scl, yimg+dy/scl),aspect=dx/dy)

        xputils.set_fontsize(ax,17.)
        # Write
        fig.tight_layout(pad=0.2,h_pad=0.,w_pad=0.1)
        pp.savefig()
        fig.close()

    # Finish
    print('tlk_coshalos: Wrote {:s}'.format(outfil))
    pp.close()

def fig_abs_gallery(wrest=None, lbl=None, outfil=None, 
    cos_halos=None, passback=False):
    # import tlk_coshalos as tch
    # reload(tch)
    # cos_halos = tch.fig_abs_gallery(passback=True)
    '''Gallery of absorption lines
    '''
    # Linelist
    llist = lll.LineList('Strong')
    if wrest is None:
        wrest = 1215.670*u.AA
    Zion = (llist[wrest]['Z'], llist[wrest]['ion'])
    if lbl is None:
        lbl = r'HI Ly$\alpha$'
    if outfil is None:
        outfil='fig_HI_gallery.pdf'
    # Init COS-Halos (slow)
    if cos_halos is None:
        cos_halos = COSHalos()
        cos_halos.load_mega()#skip_ions=True)
        if passback:
            return cos_halos

    xmnx = (9.4, 11.7) 
    ymnx = (-13., -9.)
    vmnx = (-600., 600.)
    wbox= {'facecolor':'white', 'edgecolor':'white'}
    lclr = 'blue'
    # To help with looping
    def sngl_abs(ax,field_gal,wbox=wbox,lclr=lclr,vmnx=vmnx,
        xmnx=xmnx,ymnx=ymnx,cos_halos=cos_halos, show_xylbl=True,
        srect=0.2, show_ewlbl=True):
        #
        cgm_abs = cos_halos[field_gal]
        xplt = cgm_abs.galaxy.stellar_mass
        yplt = np.log10(cgm_abs.galaxy.sfr[1])-cgm_abs.galaxy.stellar_mass
        # Spec
        spec = cos_halos.load_bg_cos_spec(field_gal, wrest)
        velo = spec.relative_vel(wrest*(cgm_abs.galaxy.z+1))
        # EW
        data = cgm_abs.abs_sys.ions[Zion]
        itrans = np.where(np.abs(data['LAMBDA']-wrest.value) < 1e-3)[0]
        if len(itrans) == 0:
            print('No measurement of {:g} for {:s} {:s}'.format(
                wrest,field_gal[0],field_gal[1]))
            return
        EW = data['WREST'][itrans][0]*u.AA/1e3
        sigEW = data['SIGWREST'][itrans][0]*u.AA/1e3
        if sigEW <= 0.: # Something must be wrong
            return
        if EW > 3*sigEW:
            sclr = 'k'
        else:
            sclr = 'r'
        #
        xrect = (xplt-xmnx[0])/(xmnx[1]-xmnx[0])
        yrect = (yplt-ymnx[0])/(ymnx[1]-ymnx[0])
        rect=[xrect, yrect, srect, srect] 
        axsp = xpsimp.add_subplot_axes(ax,rect)#,axisbg=axisbg)
        axsp.set_xlim(vmnx)
        axsp.set_ylim(-0.1, 1.3)
        if show_xylbl:
            axsp.set_xlabel('Relative Velocity (km/s)',fontsize=9.)
            axsp.set_ylabel('Normalized Flux',fontsize=9.)
        # Plot
        axsp.plot(velo, spec.flux, sclr, drawstyle='steps', lw=1.3)
        # Label
        #llbl = (r'$W_\lambda = $'+'{:0.2f} A \n'.format(EW.value)+
        #    '$z=${:.3f}\n'.format(cgm_abs.galaxy.z)+
        #    r'$\rho=$'+'{:d}kpc'.format(int(np.round(cgm_abs.rho.to('kpc').value))))
        if show_ewlbl:
            llbl = r'$W_\lambda = $'+'{:0.2f} A'.format(EW.value)
            axsp.text(-300., 1.1, llbl, ha='left', fontsize=7., 
                color=lclr, bbox=wbox)#multialignment='left')

    # Start the plot
    pp = PdfPages(outfil)


    for ss in range(4):
        plt.figure(figsize=(8, 5))
        plt.clf()
        gs = gridspec.GridSpec(1, 1)

        # Axes
        ax = plt.subplot(gs[0,0])
        ax.set_xlim(xmnx)
        ax.set_ylim(ymnx)
        ax.set_ylabel('sSFR')
        ax.set_xlabel('log M*')
        ax.text(9.75, -12., lbl, color=lclr, fontsize=21.)
        #ax.xaxis.set_major_locator(plt.MultipleLocator(300.))

        # 
        if ss == 0: # Single system
            field_gal = ('J0950+4831','177_27')
            sngl_abs(ax,field_gal)
        elif ss == 1: # Two systems
            cgm_list = [('J0950+4831','177_27'), ('J1245+3356', '236_36')]
            for field_gal in cgm_list:
                sngl_abs(ax,field_gal)
        elif ss > 1: # Half systems
            if ss == 2:
                iend = 20
            else:
                iend = -1
            for cgm_abs in cos_halos.cgm_abs[0:iend]:
                field_gal = (cgm_abs.field, cgm_abs.gal_id)
                sngl_abs(ax,field_gal,show_xylbl=False,srect=0.1,
                    show_ewlbl=False)
        else:
            pass

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

    # Experimental design
    if (flg_fig % 2**1) >= 2**0:
        fig_experiment()

    # Image gallery
    if (flg_fig % 2**2) >= 2**1:
        fig_images()

    # Abs gallery
    if (flg_fig % 2**3) >= 2**2:
        fig_abs_gallery()


# Command line execution
if __name__ == '__main__':

    if len(sys.argv) == 1: # Figs
        flg_fig = 0 
        #flg_fig += 2**0 # Experiment
        flg_fig += 2**1 # Image gallery
        #flg_fig += 2**2 # Abs gallery
    else:
        flg_fig = sys.argv[1]

    main(flg_fig)
