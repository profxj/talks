# Module for the Spectra figures for talks
# Imports
from __future__ import print_function

import numpy as np
import glob, os, sys, imp

import matplotlib as mpl
mpl.rcParams['font.family'] = 'stixgeneral'
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec

from astropy.io import ascii
from astropy.relativity import velocities as arv
from astropy import units as u
from astropy import constants as const

from xastropy.spec.lines_utils import AbsLine
from xastropy.igm.abs_sys.lls_utils import LLS_System, LLS_Survey
from xastropy.igm.abs_sys import abssys_utils as abssys
from xastropy import spec as xsp
from xastropy.plotting import simple as xplots
from xastropy.plotting import utils as xputils
from xastropy.xutils import xdebug as xdb

xa_path = imp.find_module('xastropy')[1]

# Local
#sys.path.append(os.path.abspath("../Analysis/py"))
#import lls_sample as lls_s


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
        fig_fn()


# Command line execution
if __name__ == '__main__':

    if len(sys.argv) == 1: # Figs
        flg_fig = 0 
        flg_fig += 2**0 # fN
        #flg_fig += 2**1 # COG
    else:
        flg_fig = sys.argv[1]

    main(flg_fig)
