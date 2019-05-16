# Module for the Spectra figures for talks
# Imports
from __future__ import print_function

import numpy as np
import glob, os, sys, imp
import pdb

import matplotlib as mpl
mpl.rcParams['font.family'] = 'stixgeneral'
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec

from astropy.io import ascii
from astropy import units as u
from astropy.io import fits

from linetools.spectra.xspectrum1d import XSpectrum1D
from linetools.spectra import io as lsio

from specdb.specdb import IgmSpec

# Local
sys.path.append(os.path.abspath("../py"))
from tlk_fig_utils import set_fontsize, set_spines



####
#  f(N) at z~2.5; Data + Model
def fig_fn(outfil=None, data_list=None):

    from xastropy.igm.fN import model as xifm
    from xastropy.igm.fN import data as xifd

    # Init
    if outfil == None:
        outfil='fig_fN.pdf'

    # Start the plot
    if outfil != None:
        pp = PdfPages(outfil)

    fn_file = xa_path+'/igm/fN/fn_constraints_z2.5_vanilla.fits'
    k13r13_file = xa_path+'/igm/fN/fn_constraints_K13R13_vanilla.fits'
    n12_file = xa_path+'/igm/fN/fn_constraints_N12_vanilla.fits'
    all_fN_cs = xifd.fn_data_from_fits([k13r13_file, fn_file, n12_file])

    fN_model = xifm.default_model()


    
    # Remove K12
    #data_list = ['K13R13','OPB07', 'N12']
    #outfil = 'tmp.png'
    if data_list is None:
        fN_cs = [fN_c for fN_c in all_fN_cs if ((fN_c.ref != 'K02') & (fN_c.ref != 'PW09'))]
    else:
        fN_cs = [fN_c for fN_c in all_fN_cs if fN_c.ref in data_list]
    fN_dtype = [fc.fN_dtype for fc in fN_cs]

    # Plot
    ymnx = [-26, -10.]
    xmnx = [12., 22.]

    plt.figure(figsize=(8, 5))
    plt.clf()
    gs = gridspec.GridSpec(1, 1)
    main = plt.subplot(gs[:, :])

    # f(N) data
    main.set_ylabel(r'$\log \, f(N_{\rm HI})$')
    main.set_xlabel(r'$\log \, N_{\rm HI}$')
    main.set_ylim(ymnx)
    main.set_xlim(xmnx)

    syms = ['o', 's', 'D']
    isyms = 0
    for fN_c in fN_cs: 
        if fN_c.fN_dtype == 'fN':
            # Length
            ip = range(fN_c.data['NPT'])
            #xdb.set_trace()
            val = np.where(fN_c.data['FN'][ip] > -90)[0]
            #xdb.set_trace()
            if len(val) > 0:
                #xdb.set_trace()
                ipv = np.array(ip)[val]
                xval = np.median(fN_c.data['BINS'][:,ipv],0)
                xerror = [ fN_c.data['BINS'][1,ipv]-xval, xval-fN_c.data['BINS'][0,ipv] ]
                yerror = [ fN_c.data['SIG_FN'][1,ipv], fN_c.data['SIG_FN'][0,ipv] ] # Inverted!
                main.errorbar(xval, fN_c.data['FN'][ipv], xerr=xerror, 
                    yerr=yerror, fmt=syms[isyms], 
                    label=fN_c.ref,capthick=2, color='gray')
                isyms += 1

    # Paint on the Media
    lsz = 16.
    main.fill_between( [xmnx[0],14.5], ymnx[0], ymnx[1],
        color='red', alpha=0.3) 
    main.text(13.1, -18., 'IGM\n(85%)', color='red', size=lsz, ha='center')
    main.fill_between( [14.5,19.], ymnx[0], ymnx[1],
        color='blue', alpha=0.3) 
    main.text(16.5, -22., 'CGM\n(10%)', color='blue', size=lsz, ha='center')
    main.fill_between( [19., xmnx[1]], ymnx[0], ymnx[1],
        color='green', alpha=0.3) 
    main.text(20.5, -18., 'ISM\n(5%)', color='green', size=lsz, ha='center')

    # Label
    main.text(np.sum(xmnx)/2., -12., r'$z \approx 2.5$'+'\n(P+14)', 
        color='black', size=lsz, ha='center', va='top')

    # Legend
    main.legend(loc='lower left', numpoints=1)


    # Model?
    xplt = 12.01 + 0.01*np.arange(1100)
    yplt = fN_model.eval(xplt, 2.4)
    main.plot(xplt,yplt,'-',color='black')

    # Fonts
    xputils.set_fontsize(main,16.)

    # Write
    plt.tight_layout(pad=0.2,h_pad=0.,w_pad=0.1)
    pp.savefig()
    plt.close()

    # Finish
    print('tlk_igm.fig_fn: Wrote {:s}'.format(outfil))
    pp.close()

def fig_lya_forest(zoom_in=False, redshift=False):

    outfile = 'fig_lya_forest.png'
    xmnx = (3950., 4363)

    if zoom_in:
        xmnx = (4130., 4183)
        outfile = 'fig_lya_forest_zoom.png'

    if redshift:
        outfile = 'fig_lya_forest_z.png'
        xmnx = np.array(xmnx)/1215.6701 - 1.

    # Read spectrum
    spec_file = os.getenv('DROPBOX_DIR')+'Keck/HIRES/RedData/Q1759+75/Q1759+75T_f.fits'
    xspec = XSpectrum1D.from_file(spec_file)


    plt.figure(figsize=(9, 5))
    plt.clf()
    gs = gridspec.GridSpec(1, 1)
    ax = plt.subplot(gs[:, :])

    # Plot
    if redshift:
        xplt = xspec.wavelength.value/1215.6701 - 1.
    else:
        xplt = xspec.wavelength.value
    ax.plot(xplt, xspec.flux, 'k')
    ax.plot([0.,5000], [0., 0.], 'g--')
    ax.set_xlim(xmnx)
    ax.set_ylim(-0.05, 1.1)
    ax.set_ylabel('Normalized Flux')
    if redshift:
        ax.set_xlabel('Redshift (z)')
    else:
        ax.set_xlabel('Wavelength (A)')

    set_fontsize(ax, 15.)

    # Write
    plt.tight_layout(pad=0.2,h_pad=0.,w_pad=0.1)
    plt.savefig(outfile, dpi=500)
    plt.close()
    print("Wrote: {:s}".format(outfile))

def fig_lowz_hiz(outfil='fig_loz_hiz.png'):
    """ Show varying IGM transmission
    """
    #hdlls_path = '/u/xavier/paper/LLS/Optical/Data/DR1/Spectra/'
    esi_path = '/u/xavier/Keck/ESI/RedData/'
    hst_path = '/u/xavier/HST/Cycle23/z1IGM/Archive/PG1206+459/'
    #
    igmsp = IgmSpec()
    idicts = [
        dict(filename='Data/3C273_STIS_E140M_F.fits'),
        dict(filename=hst_path+'PG1206+459_E230M_f.fits'),
        dict(coord='J212329.50-005052.9', group=['HD-LLS_DR1']),
        dict(coord='J113621.00+005021.0', group=['HD-LLS_DR1']),
        dict(coord='J113246.5+120901.6', group=['ESI_DLA']),
        dict(filename=esi_path+'J1148+5251/SDSSJ1148+5251_stack.fits'),  # z=6
        ]
    lbls = [
        'HST/STIS: 3C273',
        'HST/STIS: PG1206+459',
        'Keck/HIRES: J2123-0050',  # 2.26
        'Magellan/MIKE: J1136+0050', # 3.43
        'Keck/ESI: J1132+1209', # 5.17
        'Keck/ESI: J1148+5251', # 6.4
        ]
    zems = [0.17, 1.16, 2.26, 3.43, 5.17, 6.4]
    xrest = np.array([1080, 1200.])
    ymnx = [-0.1, 1.1]

    # Cut down
    idicts = [idicts[0], idicts[3]]
    lbls = [lbls[0], lbls[3]]
    zems = [zems[0], zems[3]]

    lw = 1.
    csz = 19.

    # Start the plot
    #fig = plt.figure(figsize=(5.0, 8.0))
    fig = plt.figure(figsize=(8.0, 5.0))

    plt.clf()
    gs = gridspec.GridSpec(len(lbls),1)

    # Loop
    for qq, lbl in enumerate(lbls):

        # Grab data
        idict = idicts[qq]
        if 'coord' in idict.keys():
            qdict = {}
            for key in idict.keys():
                if key not in ['coord','group']:
                    qdict[key] = idict[key]
            spec, meta = igmsp.spectra_from_coord(idict['coord'], tol=5.*u.arcsec,
                                                  groups=idict['group'])#, query_dict=qdict)
            if meta is None:
                print("Bad coord?")
                pdb.set_trace()
            elif len(meta) > 1:
                pdb.set_trace()
        else:
            spec = lsio.readspec(idict['filename'])

        if lbl == 'HST/STIS: 3C273':
            hdu = fits.open('Data/3C273_STIS_E140M_c.fits')
            conti_3c273 = hdu[0].data
            spec.co = conti_3c273
            spec.normed = True

        # Spectrum
        ax = plt.subplot(gs[qq])
        ax.set_xlim(xrest*(1+zems[qq])/1215.67 - 1)
        ax.set_ylim(ymnx)
        if qq == 3:
            ax.set_ylabel('Normalized Flux')
        if qq == len(lbls)-1:
            ax.set_xlabel(r'Redshift of Ly$\alpha$')


        ax.plot(spec.wavelength.value/1215.6701 - 1, spec.flux, 'k', linewidth=lw)

        # Label
        #ax.text(0.05, 0.95, lbl+' zem={:0.1f}'.format(zems[qq]), color='blue',
        #    transform=ax.transAxes, size=csz, ha='left', bbox={'facecolor':'white'})
        #
        set_fontsize(ax, 12.)

    # Layout and save
    #plt.subplots_adjust(hspace=0)
    plt.tight_layout(pad=0.2,h_pad=0.0,w_pad=0.4)
    plt.savefig(outfil, dpi=600)
    plt.close()
    # Finish
    print('Writing {:s}'.format(outfil))

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

    # f(N)
    if (flg_fig & 2**0):
        fig_fn()

    # Lya forest
    if (flg_fig & 2**1):
        fig_lya_forest()
        fig_lya_forest(zoom_in=True)
        fig_lya_forest(redshift=True)

    # Low-z vs. high-z
    if flg_fig & 2**2:
        fig_lowz_hiz()


# Command line execution
if __name__ == '__main__':

    if len(sys.argv) == 1: # Figs
        flg_fig = 0 
        #flg_fig += 2**0 # fN
        flg_fig += 2**1 # Lya forest
        flg_fig += 2**2 # lowz vs. high z
    else:
        flg_fig = sys.argv[1]

    main(flg_fig)
