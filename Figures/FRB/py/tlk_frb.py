""" Module for the figures for talks related to FRBs
"""

# Imports
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import glob, os, sys, json
from IPython import embed

#import healpy as hp

import matplotlib as mpl

mpl.rcParams['font.family'] = 'stixgeneral'

from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import image as mpimg
from matplotlib.patches import Ellipse

from pkg_resources import resource_filename

import healpy as hp

from astropy import units
from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
from astropy.visualization.wcsaxes import SphericalCircle

from scipy.special import gammainc
from scipy.integrate import quad
from scipy.optimize import fsolve

#from casbah.gal_properties import calchalomass

from frb import hosts
flg_hosts = True
from frb.galaxies import frbgalaxy
from frb import frb
from frb.figures import utils as ffutils
from frb.figures import galaxies as ffgalaxies

# Local
sys.path.append(os.path.abspath("/home/xavier/papers/FRB/Figures/py"))
import frb_figs

# Globals
handletextpad=0.3

def _bn(n):
    f = lambda x: gammainc(2*n,x)-0.5
    return fsolve(f,2*n-1/3)
def _frac_brightness(n,r,reff):
    b_n = _bn(n)
    integrand = lambda x: np.exp(-b_n*((x/reff)**(1/n)-1))
    return quad(integrand,0,r)[0]/quad(integrand,0,np.inf)[0]


def fig_lorimer_DM(outfile='fig_lorimer_DM.png'):
    """

    """
    set_mplrc()

    lorimer = frb.FRB('FRB010724', 'J011806.0-751218.0',
                        375*units.pc/units.cm**3, z_frb=0.3)

    plt.clf()
    fig = plt.figure(figsize=(14., 10))
    f, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
    #
    frb_figs.sub_cartoon(ax1, ax2, lorimer.coord, lorimer.z, host_DM=50., halos=False,
                         ymax=lorimer.DM.value, fg_halos=None)

    # Layout and save
    plt.tight_layout(pad=0.2,h_pad=0.1,w_pad=0.1)
    plt.savefig(outfile, dpi=300)
    plt.close()
    print('Wrote {:s}'.format(outfile))

def fig_DM_maps(flag, outfile):

    if flag == 1:
        map_file = resource_filename('frb', 'data/DM/DM_ISM/fits')
        cmap=plt.cm.Greens,
        title = 'ISM'
    else:
        return
    dm_map = hp.read_map(map_file)
    print('Min/max = {}/{}'.format(np.min(dm_map), np.max(dm_map)))

    fig = plt.figure(figsize=(7., 5))
    hp.mollview(dm_map, cbar=True, xsize=800, min=10, max=275., cmap=cmap,
                title=title, unit=r'DM (pc/cm$^3$)')
    gg = plt.cm.Greys(0.8)
    hp.graticule(color=gg)
    hp.projtext(270, 2, r'$270^\circ$', lonlat=True, fontsize=14, color=gg)
    hp.projtext(0, 2, r'$0^\circ$', lonlat=True, fontsize=14, color=gg)
    hp.projtext(90, 2, r'$90^\circ$', lonlat=True, fontsize=14, color=gg)
    hp.projtext(30, 32, r'$30^\circ$', lonlat=True, fontsize=14, color=gg)
    hp.projtext(26, 62, r'$60^\circ$', lonlat=True, fontsize=14, color=gg)
    hp.projtext(30, -28, r'$-30^\circ$', lonlat=True, fontsize=14, color=gg)
    hp.projtext(30, -58, r'$-60^\circ$', lonlat=True, fontsize=14, color=gg)

    print("Writing {:s}".format(outfile))
    plt.savefig(outfile)
    plt.close()


def fig_hst27_others(outfile='fig_hst27_others.pdf'):
    """

    """
    frbfigures.utils.set_mplrc()

    # FRBs and hosts
    frb181112 = frb.FRB.by_name('FRB181112')
    HG181112 = frbgalaxy.FRBHost.by_name('181112')
    frb190102 = frb.FRB.by_name('FRB190102')
    #HG190102 = frbgalaxy.FRBHost.by_name('190102')

    # Start the plot
    fig = plt.figure(figsize=(15, 7))
    plt.clf()
    fsz = 15.
    isz = 25. # Instrument label size
    xwidth = 0.38
    ywidth = 0.9
    y0 = 0.08
    ylbl = 0.91

    # FRB 181112
    for kk, ifrb, img_file in zip([0,1], [frb181112, frb190102],
                             ['../../181112/Data/FRB181112_VLT_FORS2g.fits',
                              '../../Cosmic/Data/190102_g_coadded.fits']):
        hdu = fits.open(img_file)
        header = hdu[0].header
        vlt_g = hdu[0].data

        isize = 20.  # arcsec
        size = units.Quantity((isize, isize), units.arcsec)
        cutout = Cutout2D(vlt_g, ifrb.coord, size, wcs=WCS(header))

        if kk == 0:
            axVLT = fig.add_axes([0.10, y0, xwidth, ywidth], projection=cutout.wcs)
            vmin, vmax = 6103., 7781
            cmap = plt.get_cmap('Reds')
            frb_clr = 'black'
            #c = Ellipse(xy=(ifrb.coord.ra.value, ifrb.coord.dec.value),
            #                  width=0.5/3600., height=0.1/3600.,
            #          angle=frb181112.eellipse['theta']-90.,
            #                  facecolor=frb_clr, edgecolor=frb_clr)#, lw=1, ls='--')
            c = SphericalCircle((ifrb.coord.ra, ifrb.coord.dec),
                                0.3*units.arcsec, transform=axVLT.get_transform('icrs'),
                                edgecolor=frb_clr, facecolor=frb_clr)
        else:
            axVLT = fig.add_axes([0.58, y0, xwidth, ywidth], projection=cutout.wcs)
            vmin, vmax =-30., 500.
            cmap = plt.get_cmap('Greens')
            frb_clr = 'cyan'
            c = SphericalCircle((ifrb.coord.ra, ifrb.coord.dec),
                            0.1*units.arcsec, transform=axVLT.get_transform('icrs'),
                            edgecolor=frb_clr, facecolor=frb_clr)

        lon = axVLT.coords[0]
        lat = axVLT.coords[1]
        lon.set_ticks(exclude_overlapping=True)
        lon.set_major_formatter('hh:mm:ss')
        lon.set_ticks(number=4)
        #
        d = axVLT.imshow(cutout.data, cmap=cmap, vmin=vmin, vmax=vmax)
        plt.grid(color='gray', ls='dashed')
        axVLT.set_xlabel(r'\textbf{Right Ascension (J2000)}', fontsize=fsz)
        axVLT.set_ylabel(r'\textbf{Declination (J2000)}', fontsize=fsz, labelpad=-1.)
        #axVLT.invert_xaxis()

        axVLT.add_patch(c)

        axVLT.text(0.05, ylbl, r'\textbf{'+'{}'.format(ifrb.frb_name)+'}',
                   transform=axVLT.transAxes, fontsize=isz, ha='left', color='black')


    # Layout and save
    plt.tight_layout(pad=0.2,h_pad=0.1,w_pad=0.1)
    plt.savefig(outfile, dpi=300)
    plt.close()
    print('Wrote {:s}'.format(outfile))


def fig_hst27_pop(outfile='fig_hst27_pop.pdf'):
    """
    BPT and SFR vs. Mstar
    """
    ffutils.set_mplrc()

    # FRBs and hosts
    HG121102 = frbgalaxy.FRBHost.by_name('121102')
    HG180924 = frbgalaxy.FRBHost.by_name('180924')
    HG181112 = frbgalaxy.FRBHost.by_name('181112')

    # Edit 1811122 because of Telluric
    HG181112.neb_lines['[NII] 6584'] = HG181112.neb_lines['[NII] 6584'] / 4
    #HG181112.neb_lines['Halpha_err'] = HG181112.neb_lines['Halpha_err']

    galaxies = [HG121102, HG180924, HG181112]
    clrs = ['black', 'blue', 'red']
    markers = ['s', 'o', 'o']

    # Start the plot
    fig = plt.figure(figsize=(15, 7))
    plt.clf()
    xwidth = 0.38
    ywidth = 0.9
    y0 = 0.08

    # BPT
    axBPT = fig.add_axes([0.10, y0, xwidth, ywidth])
    ffgalaxies.sub_bpt(axBPT, galaxies, clrs, markers, SDSS_clr='Greys')

    # SFR vs. Mstar
    axSFMS = fig.add_axes([0.57, y0, xwidth, ywidth])
    frb_figs.sub_sfms(axSFMS, galaxies, clrs, markers)

    # Layout and save
    plt.tight_layout(pad=0.2,h_pad=0.1,w_pad=0.1)
    plt.savefig(outfile, dpi=300)
    plt.close()
    print('Wrote {:s}'.format(outfile))

def log_me(val, err):
    xerr = np.array([[np.log10(val) - np.log10(val - err)],
                     [-np.log10(val) + np.log10(val + err)]])
    return np.log10(val), xerr

def set_fontsize(ax,fsz):
    '''
    Parameters
    ----------
    ax : Matplotlib ax class
    fsz : float
      Font size
    '''
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(fsz)


def set_mplrc():
    mpl.rcParams['mathtext.default'] = 'it'
    mpl.rcParams['font.size'] = 12
    mpl.rc('font',family='Times New Roman')
    mpl.rcParams['text.latex.preamble'] = [r'\boldmath']
    mpl.rc('text', usetex=True)


#### ########################## #########################
def main(flg_fig):

    if flg_fig == 'all':
        flg_fig = np.sum( np.array( [2**ii for ii in range(25)] ))
    else:
        flg_fig = int(flg_fig)

    # Lorimer DM cartoon
    if flg_fig & (2**0):
        #fig_hst27_180924()
        fig_lorimer_DM()

    # Other images
    if flg_fig & (2**1):
        fig_hst27_others() #'fig_hst27_others.png')

    # Other images
    if flg_fig & (2**2):
        fig_hst27_pop() #'fig_hst27_pop.png')



# Command line execution
if __name__ == '__main__':

    if len(sys.argv) == 1:
        flg_fig = 0
        #flg_fig += 2**0   # Lorimer DM
        flg_fig += 2**2   # DM maps
    else:
        flg_fig = sys.argv[1]

    main(flg_fig)
