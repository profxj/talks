# Module for comparitive plots on the CGM
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
from xastropy.cgm.cos_halos import COSHalos, COSDwarfs
from xastropy.spec.lines_utils import AbsLine
from xastropy.obs import finder
from xastropy.obs import radec as xra
from xastropy.plotting import utils as xputils
from xastropy.xutils import xdebug as xdb

# Local
#sys.path.append(os.path.abspath("../Analysis/py"))
#import lls_sample as lls_s



####
#  Plot illustrating COS-Dwarfs and COS-Halos 
def fig_cos_dwarfs_images(outfil=None):

    from xastropy.obs import x_getsdssimg as xgs

    # Init COS-Dwarfs, COS-Halos 
    cos_dwarfs = COSDwarfs()
    cos_dwarfs.load_mega(skip_ions=True)
    cos_halos = COSHalos()
    cos_halos.load_mega(skip_ions=True)

    # Start the plot
    if outfil is None:
        outfil='fig_cos_dwarfs_images.pdf'
    pp = PdfPages(outfil)

    dx = 0.2
    dy = 0.40

    fig = plt.figure(figsize=(8, 5))
    fig.clf()
    # Setup for dark
    xpsimp.dark_bkgd(plt)
    gs = gridspec.GridSpec(1, 1)

    # Axes
    ax = plt.subplot(gs[0,0])
    ax.set_xlim(8.0, 11.7) 
    #ax.xaxis.set_major_locator(plt.MultipleLocator(300.))
    ax.set_ylim(-13,-8.2)
    ax.set_ylabel('sSFR')
    ax.set_xlabel('log M*')

    iend = -1
    #iend = 5
    imsize = 0.3 # arcmin
    for cgm_abs in cos_halos.cgm_abs[0:iend]:
        img, oBW = xgs.getimg(cgm_abs.galaxy.coord.ra.deg, 
            cgm_abs.galaxy.coord.dec.deg, imsize=imsize)
        # Figure out placement
        ximg = cgm_abs.galaxy.stellar_mass
        yimg = np.log10(cgm_abs.galaxy.sfr[1])-cgm_abs.galaxy.stellar_mass
        # Show
        scl=2.
        ax.imshow(img,extent=(ximg-dx/scl, ximg+dx/scl, yimg-dy/scl, yimg+dy/scl),aspect=dx/dy)

    iend = -1
    #iend = 5
    imsize = 0.8 # arcmin
    for cgm_abs in cos_dwarfs.cgm_abs[0:iend]:
        img, oBW = xgs.getimg(cgm_abs.galaxy.coord.ra.deg, 
            cgm_abs.galaxy.coord.dec.deg, imsize=imsize)
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
    plt.close()

    # Finish
    print('tlk_coshalos: Wrote {:s}'.format(outfil))
    pp.close()

####
#  Plot of all halos
def fig_all_halo_images(outfil=None):

    from xastropy.obs import x_getsdssimg as xgs

    # Init COS-Dwarfs, COS-Halos 
    cos_dwarfs = COSDwarfs()
    cos_dwarfs.load_mega(skip_ions=True)
    cos_halos = COSHalos()
    cos_halos.load_mega(skip_ions=True)

    # Start the plot
    if outfil is None:
        outfil='fig_all_halo_images.pdf'
    pp = PdfPages(outfil)

    dx = 0.2
    dy = 0.4

    fig = plt.figure(figsize=(8, 5))
    fig.clf()
    # Setup for dark
    xpsimp.dark_bkgd(plt)
    gs = gridspec.GridSpec(1, 1)

    # Axes
    ax = plt.subplot(gs[0,0])
    ax.set_xlim(10.0, 15.0) 
    #ax.xaxis.set_major_locator(plt.MultipleLocator(300.))
    ax.set_ylim(-13,-8.2)
    ax.set_ylabel('sSFR')
    ax.set_xlabel(r'$\log M_{\rm Halo}$')

    iend = -1
    #iend = 5
    imsize = 0.3 # arcmin
    for cgm_abs in cos_halos.cgm_abs[0:iend]:
        img, oBW = xgs.getimg(cgm_abs.galaxy.coord.ra.deg, 
            cgm_abs.galaxy.coord.dec.deg, imsize=imsize)
        # Figure out placement
        ximg = cgm_abs.galaxy.halo_mass
        yimg = np.log10(cgm_abs.galaxy.sfr[1])-cgm_abs.galaxy.stellar_mass
        # Show
        scl=2.
        ax.imshow(img,extent=(ximg-dx/scl, ximg+dx/scl, yimg-dy/scl, yimg+dy/scl),aspect=dx/dy)

    iend = -1
    #iend = 5
    imsize = 0.8 # arcmin
    for cgm_abs in cos_dwarfs.cgm_abs[0:iend]:
        img, oBW = xgs.getimg(cgm_abs.galaxy.coord.ra.deg, 
            cgm_abs.galaxy.coord.dec.deg, imsize=imsize)
        # Figure out placement
        ximg = cgm_abs.galaxy.halo_mass
        yimg = np.log10(cgm_abs.galaxy.sfr[1])-cgm_abs.galaxy.stellar_mass
        # Show
        scl=2.
        ax.imshow(img,extent=(ximg-dx/scl, ximg+dx/scl, yimg-dy/scl, yimg+dy/scl),aspect=dx/dy)

    xputils.set_fontsize(ax,17.)
    # Write
    fig.tight_layout(pad=0.2,h_pad=0.,w_pad=0.1)
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

    # COS-Dwarf Images
    if (flg_fig % 2**1) >= 2**0:
        fig_cos_dwarfs_images()

    # LRG finder
    if (flg_fig % 2**2) >= 2**1:
        # LRG Finder (out of order to avoid PDF issues)
        coord = xra.to_coord( (176.1888,7.24901) )
        finder.main(['LRG', coord], imsize=1.*u.arcmin, show_circ=False)

    # z~1 QPQ finder
    if (flg_fig % 2**3) >= 2**2:
        # LRG Finder (out of order to avoid PDF issues)
        coord = xra.to_coord( (26.62542, 0.25583))
        finder.main(['z1QPQ', coord], imsize=1.*u.arcmin, show_circ=False)

    # All Halos
    if (flg_fig % 2**4) >= 2**3:
        fig_all_halo_images()


# Command line execution
if __name__ == '__main__':

    if len(sys.argv) == 1: # Figs
        flg_fig = 0 
        #flg_fig += 2**0 # COS-Dwarf images
        #flg_fig += 2**1 # LRG finder
        #flg_fig += 2**2 # z~1 QPQ
        flg_fig += 2**3 # All Halos
    else:
        flg_fig = sys.argv[1]

    main(flg_fig)
