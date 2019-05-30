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


def fig_lorimer_DM(outfile='fig_lorimer_DM.png', z_frb=0.3):
    """
    DM Cartoon for the Lorimer burst

    """
    set_mplrc()

    lorimer = frb.FRB('FRB010724', 'J011806.0-751218.0',
                        375*units.pc/units.cm**3, z_frb=z_frb)

    plt.clf()
    fig = plt.figure(figsize=(14., 10))
    f, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
    #
    frb_figs.sub_cartoon(ax1, ax2, lorimer.coord, lorimer.z, host_DM=50., halos=False,
                         FRB_DM=lorimer.DM.value, fg_halos=None, yscl=0.97)

    # Layout and save
    plt.tight_layout(pad=0.2,h_pad=0.1,w_pad=0.1)
    plt.savefig(outfile, dpi=300)
    plt.close()
    print('Wrote {:s}'.format(outfile))


def fig_repeater_DM(outfile='fig_repeater_DM.png'):
    """

    """
    set_mplrc()

    repeater = frb.FRB.by_name('FRB121102')

    plt.clf()
    fig = plt.figure(figsize=(14., 10))
    f, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
    #
    frb_figs.sub_cartoon(ax1, ax2, repeater.coord, repeater.z, host_DM=150., halos=False,
                         ymax=repeater.DM.value, fg_halos=None)

    # Layout and save
    plt.tight_layout(pad=0.2,h_pad=0.1,w_pad=0.1)
    plt.savefig(outfile, dpi=300)
    plt.close()
    print('Wrote {:s}'.format(outfile))

def fig_DM_maps(flag, outfile):

    if flag == 1:
        map_file = resource_filename('frb', 'data/DM/hp_DM_ISM.fits')
        cmap=plt.get_cmap('Greens')
        title = 'ISM'
    elif flag == 2:
        map_file = resource_filename('frb', 'data/DM/hp_DM_LG.fits')
        cmap=plt.get_cmap('Blues')
        mnmx = 20., 300.
    else:
        return
    dm_map = hp.read_map(map_file)
    print('Min/max = {}/{}'.format(np.min(dm_map), np.max(dm_map)))

    fig = plt.figure(figsize=(7., 5))
    hp.mollview(dm_map, cbar=True, xsize=800, min=mnmx[0], max=mnmx[1], cmap=cmap,
                unit=r'DM (pc/cm$^3$)', norm='log', title=None)
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
    plt.savefig(outfile, dpi=300)
    plt.close()



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
        fig_lorimer_DM()

    # Other images
    if flg_fig & (2**1):
        fig_repeater_DM()

    # DM maps
    if flg_fig & (2**2):
        #fig_DM_maps(1, 'fig_DM_map_ISM.png')
        fig_DM_maps(2, 'fig_DM_map_LG.png')


# Command line execution
if __name__ == '__main__':

    if len(sys.argv) == 1:
        flg_fig = 0
        #flg_fig += 2**0   # Lorimer DM
        flg_fig += 2**1   # Repeater DM
        #flg_fig += 2**2   # DM maps
    else:
        flg_fig = sys.argv[1]

    main(flg_fig)
