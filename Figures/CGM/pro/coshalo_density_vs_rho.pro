
pro coshalo_density_vs_rho, ionz, ioni, galprop, psfile = psfile, white = white,  anyion = anyion, eqw = eqw, xion = xion, ionratio = ionratio, XPS=xps, MEGASTRUCT=megastruct

;;;;; Program to plot any column density or ratio of column densities
;;;;; vs. any galaxy property using the megastructures. 
;;;;; INPUTS: 
;;;;; ionz -- the atomic number of the element, an integer
;;;;; ioni -- the ionic state, an integer
;;;;; e.g., If you want plots of Si III, then ionz = 14 and ioni = 3

;;;;; GALPROP MUST BE ONE OF THE FOLLOWING (i.e. the quantity on the X-axis)
;;;;; galprop is a string of any one of the following: 'smass' for
;;;;; stellar mass (actually a log); 'rhokpc' for impact parameter in kpc; 'sfr' for
;;;;; the star formation rate in msun per yr; 'metal' for galaxy
;;;;; metallicity, this will be the 12 + Log (O/H) abundance; and
;;;;; 'nhi' for Log HI column density in cm^-2, and 'nion' for a
;;;;; column density. If this is set, then you must specify xion =
;;;;; [x_ionz, x_ioni] otherwise, things break badly. 
;;;;;
;;;;; IF YOU WANT A COLUMN ON THE X-AXIS, then set xion = [x_ionz,
;;;;; xioni] with galprop = 'nion'
;;;;;
;;;;; IF you'd like an ionratio, then ionratio should be set to
;;;;; the two element vector of the ion you'd like to have in
;;;;; the demoninator of the ion ratio plot, i.e. ionratio = [14,3]
;;;;; and ionz = 14 and ioni = 2 for SiII/SiIII plots. 

;;;;; Keyword WHITE is set for white background plot. Otherwise, will
;;;;; make plots with dark background. 

;;;;; Set psfile to the name of the plot you're creating. If
;;;;; you do not, the plot name is set automatically to be something
;;;;; sensible. 

;;;;; Added Keyword EQW -- setting this keyword will give you
;;;;;                      equivalent widths on the Y-axis for
;;;;;                      a specific transition. Still set the ionz
;;;;;                      and ioni. Currently, SiII will be 1260,
;;;;;                      SiIII is 1206, CII is 1036, CIII is 977,
;;;;;                      NII is 1083 and NIII is 989 and HI is
;;;;;                      1215 MgII is 2796. 

;;;;; Added Keyword ANYION -- INACTIVE KEYWORD. 
;;;;;                         this will give you the maximum column
;;;;;                         density of any low ion (CII, CIII,
;;;;;                         SiIII, SiII, NII, NIII) if it is
;;;;;                         detected. In effect, a detection
;;;;;                         statistic for the lower ionization
;;;;;                         states. Non-detections are plotted as
;;;;;                         smaller, open symbols. The ionz and ioni
;;;;;                         must still be input, but they are
;;;;;                         completely irrelevant. perhaps this bug
;;;;;                         could be worked out, since it is cludgy
;;;;;                         but as of 4/27/2012 it's not worth it. 

;;;;; Keyword XION -- This is set the same way as the inputs and
;;;;;                 ionratio, by atomic mass and ionization
;;;;;                 state. If this is set to [14,2], for instance,
;;;;;                 then the SiII column density will be on the
;;;;;                 xaxis.


 xps = 1



;;;;;;; SET PLOT NAME ;;;;;;;;
 psfile = 'coshalo_density_vs_rho.ps'

;;;;;;; THE DEFAULTS SHOULD BE SENSIBLE ;;;;;

ldir = getenv('DROPBOX_DIR')+'/COS-Halos/lowions/'
tdir = getenv('DROPBOX_DIR')+'/COS-Halos/Targets/'

cldy_str=xmrdfits('coshaloscloudysol.fits',1)

;pushd, ldir
   


;;red/blue cut
red = where(cldy_str.sfruplim eq 'yes', nred, complement = blue, ncomplement = nblue)

;;;;;;;;;;;;;;;;;; GET XY RANGES, AXIS LABELS AND
;;;;;;;;;;;;;;;;;; QUANTITIES;;;;;;;;;;;;;;;;;;


;; impact parameter
xrange = [0.0, 200.0]
xlabel = greek('rho', /FORCE_PS)+'!X [kpc]'
xquant =  cldy_str.rhofinal
xuplim = -1
xlowlim = -1
nxup = 0
nxlow = 0


;; here we are plotting straight up column densities. 
yquant = cldy_str.hden
yrange = [-5, -1]
ylabel = 'log n!dH!N (cm!u-3!N)'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;              PLOTTING HERE         ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;; PLOT ALL GALAXIES ;;;;;;;;;;;;;;;;
if not keyword_set(XPS) then ct_psopen, psfile, xs = 10, ys = 8 else $
	x_psopen, psfile, /maxs
;ct_load_colours
 clr = getcolor(/load)

;;Set up colour indices. These correspond to colours defined by jw_setclrs
if keyword_set(white) then begin
    xclr = jw_setclrs(/white)
   bluecol = xclr[1]
   redcol = xclr[2]
   fgcol = 0
   bgcol = 255
endif else begin
   xclr = jw_setclrs(/dark)
   bluecol = clr.cyan
   redcol = clr.tomato
   fgcol = 255
   bgcol = 0
endelse 

xmrg = [8,2]
ymrg = [4,1]

ct_psfill_black, col = bgcol
plot, [0], [0], /nodata, yr = yrange, ys = 1, xr = xrange, xs = 1, $
      ytit = ylabel, xtit = xlabel, /noerase, col = fgcol, $
      chars = 2.5, xmar=xmrg, ymar=ymrg


;; Plot the detections. Don't just plot everything, since we
;; want to only outline non-detections
   det = lindgen(n_elements(red))
   oplot, xquant[red[det]], yquant[red[det]], psym = symcat(14, col = fgcol), syms = 3.0
   oplot, xquant[red[det]], yquant[red[det]], psym = symcat(14, col = redcol), syms = 2.0

   det2 = lindgen(n_elements(blue))
   oplot, xquant[blue[det2]], yquant[blue[det2]], psym = symcat(15, col = fgcol), syms = 2.3
   oplot, xquant[blue[det2]], yquant[blue[det2]], psym = symcat(15, col = bluecol), syms = 1.5

;;Panel label for HST proposal
x_psclose

;popd
end
