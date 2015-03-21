
pro fig_sivovi_vs_nvovi, ionz, ioni, galprop, psfile = psfile, white = white,  $ 
       anyion = anyion, eqw = eqw, xion = xion, ionratio = ionratio, pxionratio = pxionratio, $
       topxionz = topxionz, topxioni = topxioni, XPS=xps, xdev=xdev,  tagid = tagid, dwarf = dwarf, inyrange = inyrange, $
                      inxrange = inxrange, _extra=_extra 

;;; HACKED version of jw_plotempirical to run nv ratios. 
;;; run with the following command: 

if not keyword_set(ionz) then ionz = 14
if not keyword_set(ioni) then ioni = 4
if not keyword_set(galprop) then galprop = 'nionratio'
if not keyword_set(ionratio) then ionratio = [8,6]
if not keyword_set(pxionratio) then pxionratio = [8,6]
if not keyword_set(xion) then xion = [8,6]
if not keyword_set(topxionz) then topxionz = 7
if not keyword_set(topxioni) then topxioni = 5
WHITE = 1
if not keyword_set(inyrange) then inyrange = [-2.0, 0.]
if not keyword_set(inxrange) then inxrange = [-1.5, 0.]

;;;; 
;;;;;jw_plotnvoviratios, 14, 4, 'nionratio', ionratio=[8,6], pxionratio = [8,6], xion = [8,6], topxionz = 7, topxioni = 5, /white, inyrange = [-3.5, 0.5], inxrange = [-2.5, 0.2]
;;;;

;;;;;;;;;;;;;;;;;;;;;;;;
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

;;;; GALPROP CAN NOW BE 'qsoaz' or 'mhalo' or 'rvir' They
;;;;are what you'd  think. ADDED 1/11/13
;;;;;
;;;;; IF YOU WANT A COLUMN ON THE X-AXIS, then set xion = [x_ionz,
;;;;; xioni] with galprop = 'nion'
;;;;;
;;;;;
;;;;; IF you'd like an ionratio, then ionratio should be set to
;;;;; the two element vector of the ion you'd like to have in
;;;;; the demoninator of the ion ratio plot, i.e. ionratio = [14,3]
;;;;; and ionz = 14 and ioni = 2 for SiII/SiIII plots. 

;;;;; IF you'd like an ionratio on the xaxis, then you are a
;;;;; little crazy, but set galprop to 'nionratio', and then xionratio
;;;;; = [14,3] or whatever, and xionz = 14 and xioni = 2. This is
;;;;; fairly complex...

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

;;;;; Keyword XDEV -- Will not create PS files, will instead 
;;;;;  write plots to the x device - useful for overplotting other stuff. 
;;;;;  JT 050212 

;;;;; Keyword TAGID -- Will place a special symbol (right now a thick circle) 
;;;;;  			over points with the given subscripts. Tagid should 
;;;;;  			have values that subscript for the 51 galaxies in the 
;;;;; 			masterstructure 

;;;;; Keyword dwarfs -- works on the dwarf megastructure! JW added 5/24/2012

; if (N_params() LT 3) then begin 
;    print,'Syntax - ' + $
;             'jw_plotempirical,ionz,ioni,galproperty, psfile = '', /WHITE, /ANYION, IONRATIO = [z,i], /XPS (v1.1)'
;    return
;  endif 

 if (~keyword_set(xion) AND galprop eq 'nion') THEN BEGIN
    print, 'If you want a column density on the xaxis, please specify which ion by xion = [x_ioni, x_ionz]'
    return
 ENDIF

;;;;;;; SET PLOT NAME ;;;;;;;;
;;;;;;; THE DEFAULTS SHOULD BE SENSIBLE ;;;;;

 program = 'COS-Halos' 
 plotprefix = 'halos_'
   if (keyword_set(dwarf)) then program = 'COS-Dwarfs' 

   if (keyword_set(dwarf)) then plotprefix = 'dwarfs_' 


if (~keyword_set(psfile) AND ~keyword_set(ionratio)) then psfile = plotprefix+strtrim(string(ionz, format = '(i2)'), 2)+$
                                      '_'+strtrim(string(ioni, format = '(i1)'), 2)+'_vs_'+$
                                      strtrim(string(galprop, format = '(a8)'), 2)+'.eps'
if (keyword_set(ionratio) AND ~keyword_set(psfile)) then psfile = plotprefix+$
   strtrim(string(ionz, format = '(i2)'), 2)+$
   '_'+strtrim(string(ioni, format = '(i1)'), 2)+$
   '-'+strtrim(string(ionratio[0], format = '(i2)'), 2)+$
   '_'+strtrim(string(ionratio[1], format = '(i1)'), 2)+$ $
   '_vs_'+strtrim(string(galprop, format = '(a8)'), 2)+'.eps'

if (~keyword_set(ionratio) AND ~keyword_set(xion) AND keyword_set(eqw)) then begin
   psfile = plotprefix+strtrim(string(ionz, format = '(i2)'), 2)+$
  '_'+strtrim(string(ioni, format = '(i1)'), 2)+'EW'+'_vs_'+$
   strtrim(string(galprop, format = '(a8)'), 2)+'.eps'
endif

if (~keyword_set(ionratio) AND keyword_set(xion)) then psfile = plotprefix+strtrim(string(ionz, format = '(i2)'), 2)+$
  '_'+strtrim(string(ioni, format = '(i1)'), 2)+'_vs_'+$
   strtrim(string(galprop, format = '(a8)'), 2)+'_'+$
   strtrim(string(xion[0], format = '(i2)'), 2)+'_'+$
   strtrim(string(xion[1], format = '(i2)'), 2)+'.eps'

if (keyword_set(ionratio) AND keyword_set(xion)) then psfile = plotprefix+$
   strtrim(string(ionz, format = '(i2)'), 2)+$
   '_'+strtrim(string(ioni, format = '(i1)'), 2)+$
   '-'+strtrim(string(ionratio[0], format = '(i2)'), 2)+$
   '_'+strtrim(string(ionratio[1], format = '(i1)'), 2)+$ 
   '_vs_'+strtrim(string(galprop, format = '(a8)'), 2)+'_'+$
   strtrim(string(xion[0], format = '(i2)'), 2)+'_'+$
   strtrim(string(xion[1], format = '(i2)'), 2)+'.eps'

if ( ~keyword_set(ionratio) AND keyword_set(xion) AND keyword_set(eqw)) then psfile = plotprefix+strtrim(string(ionz, format = '(i2)'), 2)+$
  '_'+strtrim(string(ioni, format = '(i1)'), 2)+'EW'+'_vs_'+$
   strtrim(string(galprop, format = '(a8)'), 2)+'_'+$
   strtrim(string(xion[0], format = '(i2)'), 2)+'_'+$
   strtrim(string(xion[1], format = '(i2)'), 2)+'.eps'

;;;; All for the new ratio utility on xaxis


if (~keyword_set(ionratio) AND keyword_set(pxionratio)) then psfile = plotprefix+strtrim(string(ionz, format = '(i2)'), 2)+$
  '_'+strtrim(string(ioni, format = '(i1)'), 2)+'_vs_'+$
   strtrim(string(galprop, format = '(a8)'), 2)+'_'+$
   strtrim(string(topxionz, format = '(i2)'), 2)+'_'+$
   strtrim(string(topxioni, format = '(i2)'), 2)+'_over_'+$
    strtrim(string(xion[0], format = '(i2)'), 2)+'_'+$
   strtrim(string(xion[1], format = '(i2)'), 2)+'.eps'

if (keyword_set(ionratio) AND keyword_set(pxionratio)) then psfile = plotprefix+$
   strtrim(string(ionz, format = '(i2)'), 2)+$
   '_'+strtrim(string(ioni, format = '(i1)'), 2)+$
   '-'+strtrim(string(ionratio[0], format = '(i2)'), 2)+$
   '_'+strtrim(string(ionratio[1], format = '(i1)'), 2)+$ 
   '_vs_'+strtrim(string(galprop, format = '(a8)'), 2)+'_'+$
   strtrim(string(topxionz, format = '(i2)'), 2)+'_'+$
   strtrim(string(topxioni, format = '(i2)'), 2)+'_over_'+$
   strtrim(string(xion[0], format = '(i2)'), 2)+'_'+$
   strtrim(string(xion[1], format = '(i2)'), 2)+'.eps'

if ( ~keyword_set(ionratio) AND keyword_set(pxionratio) AND keyword_set(eqw)) then psfile = plotprefix+strtrim(string(ionz, format = '(i2)'), 2)+$
  '_'+strtrim(string(ioni, format = '(i1)'), 2)+'EW'+'_vs_'+$
   strtrim(string(galprop, format = '(a8)'), 2)+'_'+$
   strtrim(string(topxionz, format = '(i2)'), 2)+'_'+$
   strtrim(string(topxioni, format = '(i2)'), 2)+'_over_'+$
   strtrim(string(xion[0], format = '(i2)'), 2)+'_'+$
   strtrim(string(xion[1], format = '(i2)'), 2)+'.eps'
;;;;;;;;;;;;;;;


if keyword_set(anyion) then begin
   strput, psfile, '_any', rstrpos(psfile, '.eps')
   psfile = psfile + '.eps'
   
endif

;;;;;; END SET PLOT NAME ;;;;;;;;

psfile = 'fig_sivovi_vs_nvovi.ps'

ldir = getenv('DROPBOX_DIR')+'/'+program+'/lowions/'
tdir = getenv('DROPBOX_DIR')+'/'+program+'/Targets/'

restore, ldir+'/cosmetals_megastructure.sav' 

;pushd, ldir
   


;;red/blue cut
;red = where(megastruct.galaxy.sfr_uplim eq 'yes', nred, complement = blue, ncomplement = nblue)

red = where(megastruct.ion[8,6].nflg EQ 3, nred, complement = blue, ncomplement = nblue)

;; For dwarfs. If there are no red galaxies, as there will not be, set
;; red = 0. Too bad. This will have to be until the megastructure
;; really gets updated. 

If (red[0] ge 0) THEN BEGIN
print, 'Okay, this is fine. There are some red galaxies.'
ENDIF ELSE BEGIN
red = 0 
print, 'Warning. Majorly cludgy!! Basically you have no red galaxies, so you just made one up as index 0 to save time. Do this properly, PLEASE!'
color = megastruct.galaxy.sdssu - megastruct.galaxy.sdssr
red = where(color gt 1.8, nred, complement = blue, ncomplement = nblue)
ENDELSE

;;; Now truly define red/blue cut for dwarfs. Need to get rid of pesky
;;; AGN-dominated ones. 

IF (program eq 'COS-Dwarfs') THEN BEGIN
sflimit = 0.06
red = where(megastruct.galaxy.sfr le sflimit, nred, complement = blue, ncomplement = nblue )
megastruct[red].galaxy.sfr_uplim = 'yes'
ENDIF ELSE BEGIN
sflimit = 0.12
ENDELSE


;;;;;;;;;;;;;;;;;; GET XY RANGES, AXIS LABELS AND
;;;;;;;;;;;;;;;;;; QUANTITIES;;;;;;;;;;;;;;;;;;

;; First, get the quantities. 

IF (~keyword_set(xion)) THEN BEGIN

;;; X-axis is the galaxy property ;;;

;;;;; galprop is a string of any one of the following: 'smass' for
;;;;; stellar mass (actually a log); 'rhokpc' for impact parameter in kpc; 'sfr' for
;;;;; the star formation rate in msun per yr; 'metal' for galaxy
;;;;; metallicity, this will be the 12 + Log (O/H) abundance; 

IF (galprop NE 'smass' AND galprop NE 'rhokpc' AND galprop NE 'sfr' AND galprop NE 'metal' AND galprop NE 'nhi' AND galprop NE 'nion' $
   AND galprop NE 'rvir' AND galprop NE 'qsoaz' AND galprop NE 'mhalo') THEN BEGIN
print, 'galprop must be one of the following: smass, rhokpc, sfr, metal, or nhi or now, rvir, qsoaz, or mhalo'
return
ENDIF

IF (galprop EQ 'smass') THEN BEGIN

;; log stellar mass
xrange = [9.0, 11.5]
xlabel = 'Log M!d*!n [M!dsun!n]'
xquant =  megastruct.logmfinal
xuplim = -1
xlowlim = -1
xlim3 = -1
nxlim3 = 0
nxup = 0
nxlow = 0

ENDIF

IF (galprop EQ 'mhalo') THEN BEGIN

;; log halo mass
xrange = [10.5, 13.0]
xlabel = 'Log M!dHalo!n [M!dsun!n]'
xquant =  megastruct.logmhalo
xuplim = -1
xlowlim = -1
xlim3 = -1
nxlim3 = 0
nxup = 0
nxlow = 0

ENDIF

IF (galprop EQ 'rhokpc') THEN BEGIN

;; impact parameter
xrange = [0.0, 200.0]
xlabel = greek('rho', /FORCE_PS)+'!X [kpc]'
xquant =  megastruct.rhofinal
xuplim = -1
xlowlim = -1
xlim3 = -1
nxlim3 = 0
nxup = 0
nxlow = 0
ENDIF

IF (galprop EQ 'rvir') THEN BEGIN

;; r/rvirial radius of the galaxy
xrange = [0.0, 1.0]
xlabel = greek('rho', /FORCE_PS)+'!X / R!dvir!n'
xquant =  (megastruct.rhofinal / megastruct.rvir)
xuplim = -1
xlowlim = -1
xlim3 = -1
nxlim3 = 0
nxup = 0
nxlow = 0
ENDIF

IF (galprop EQ 'qsoaz') THEN BEGIN

;; r/rvirial radius of the galaxy
xrange = [0.0, 90.0]
xlabel = greek('Phi', /FORCE_PS)+'!X [deg]'
xquant =  megastruct.qsoaz
xuplim = -1
xlowlim = -1
xlim3 = -1
nxlim3 = 0
nxup = 0
nxlow = 0
ENDIF


IF (galprop EQ 'sfr') THEN BEGIN

;; Log SFR 
xrange = [-1.7, 1.4]
IF (program eq 'COS-Dwarfs') THEN BEGIN
xrange = [-2.9, 1.4]
ENDIF

xlabel = 'Log SFR [M!dsun!n yr!u-1!n]'
xquant =  alog10(megastruct.galaxy.sfr)
IF  (program eq 'COS-Dwarfs') THEN BEGIN
xuplim = where(megastruct.galaxy.sfr le sflimit, nxup)
ENDIF ELSE BEGIN
xuplim = where(megastruct.galaxy.sfr_uplim eq 'yes', nxup)
ENDELSE

xlowlim = -1
nxlow = 0
xlim3 = -1
nxlim3 = 0

ENDIF

IF (galprop EQ 'metal') THEN BEGIN

;; Galaxy Abundance
xrange = [-0.5, 0.5]
xlabel = 'Log Z [Z/Z!dsun!n]'
totalo = 10.^(megastruct.galaxy.abun - 12.0)
totalsun = 10.^(8.69 - 12.0) ;solar oxygen abundance. 
xquant =  alog10(totalo/totalsun)
xuplim = -1
xlowlim = -1
nxup = 0
nxlow = 0
xlim3 = -1
nxlim3 = 0
ENDIF

IF (galprop EQ 'nhi') THEN BEGIN

;; HI Column
xrange = [12.0, 19.0]
xlabel = 'Log N!dHI!n [cm!u-2!n]' 
xquant =  megastruct.ion[1,1].lognion
xuplim = where(megastruct.ion[1,1].nflg eq 3, nxup)
xlowlim = where(megastruct.ion[1,1].nflg eq 2, nxlow)
xlim3 = -1
nxlim3 = 0
IF NOT keyword_set(inxrange) THEN BEGIN
xrange = [xmin, xmax]
ENDIF ELSE BEGIN
xrange = inxrange
ENDELSE
ENDIF

ENDIF ELSE BEGIN

;; Now, we'd like the X axis quantities to be able to be an ion
;; column. Exactly as I define the Y axis. No ratios. That is way too
;; complicated. 

;; First, just ionic column densities. 

IF (galprop EQ 'nion') THEN BEGIN
;; here we are plotting straight up column densities. 
xionz = xion[0]
xioni = xion[1]
xquant = megastruct.ion[xionz, xioni].lognion
getabnd, xionnm, xionz, abnd, flag = 1
IF xioni eq 1 then strxioni = 'I'
IF xioni eq 2 then strxioni = 'II'
IF xioni eq 3 then strxioni = 'III'
IF xioni eq 4 then strxioni = 'IV'
IF xioni eq 5 then strxioni = 'V'
IF xioni eq 6 then strxioni = 'VI'
xlabel = 'Log N!d'+xionnm+strxioni+'!n [cm!u-2!n]'
goodcol = where(xquant gt 0.)
xmin = min(xquant[goodcol]) - 0.5
xmax = max(xquant[goodcol]) + 0.5
xrange = [xmin, xmax]

;; Plot the detections. Don't just plot everything, since we
;; want to only outline non-detections
;; but also account for the "detected blends" that are flagged as
;; upper limits. 
blendsredx = where((megastruct[red].ion[xionz,xioni].nflg eq 3 AND megastruct[red].ion[xionz,xioni].trans[0].flg eq 3) OR $
                 (megastruct[red].ion[xionz,xioni].nflg eq 3 AND megastruct[red].ion[xionz,xioni].trans[1].flg eq 3) OR $
                 (megastruct[red].ion[xionz,xioni].nflg eq 3 AND megastruct[red].ion[xionz,xioni].trans[2].flg eq 3), nblendredx)
blendsbluex =  where((megastruct[blue].ion[xionz,xioni].nflg eq 3 AND megastruct[blue].ion[xionz,xioni].trans[0].flg eq 3) OR $
                 (megastruct[blue].ion[xionz,xioni].nflg eq 3 AND megastruct[blue].ion[xionz,xioni].trans[1].flg eq 3) OR $
                 (megastruct[blue].ion[xionz,xioni].nflg eq 3 AND megastruct[blue].ion[xionz,xioni].trans[2].flg eq 3), nblendbluex)



detx = where(megastruct[red].ion[xionz, xioni].nflg eq 1 or megastruct[red].ion[xionz, xioni].nflg eq 2, ndetx)
det2x = where(megastruct[blue].ion[xionz, xioni].nflg eq 1 or megastruct[blue].ion[xionz, xioni].nflg eq 2, ndet2x)

;;Add arrows for the Upper limits (non-detections)
xuplim = where(megastruct.ion[xionz, xioni].nflg eq 3, nxup)
xlowlim = where(megastruct.ion[xionz, xioni].nflg eq 2, nxlow)

;;;; For consistency, these have to be defined. 
blendsredratx = -1
blendsblueratx = -1
nblendredratx = 0
nblendblueratx = 0

lim3 = -1
nlim3 = 0 ; for consistency with stupid where. 
limplotx = -1 
nlimplotx = 0
xlim3 = -1
nxlim3 = 0
ENDIF 
;;;;;;; EDIT THIS HERE> 
IF (galprop EQ 'nionratio') THEN BEGIN
iontop = megastruct.ion[topxionz,topxioni].lognion
topflg = megastruct.ion[topxionz,topxioni].nflg
ionbot =  megastruct.ion[xion[0],xion[1]].lognion
botflg = megastruct.ion[xion[0],xion[1]].nflg
xquant = iontop - ionbot
topy = megastruct.ion[ionz,ioni].nflg
;goodcol = where(iontop NE 0. AND ionbot NE 0. AND megastruct.ion[xion[0],xion[1]].nflg NE 3 $
    ;            AND megastruct.ion[8,6].nflg NE 3)

goodcol = where(iontop NE 0. AND ionbot NE 0. $
             AND megastruct.ion[8,6].nflg NE 3)
goodval = where(topflg EQ 1 AND topy EQ 1 and megastruct.ion[8,6].nflg NE 3)

xmin = min(xquant[goodcol]) - 0.2
xmax = max(xquant[goodcol]) + 0.2
IF NOT keyword_set(inxrange) THEN BEGIN
xrange = [xmin, xmax]
ENDIF ELSE BEGIN
xrange = inxrange
ENDELSE

;;; labels ;;;

getabnd, xionnm, topxionz, abnd, flag = 1
IF topxioni eq 1 then strioni = 'I'
IF topxioni eq 2 then strioni = 'II'
IF topxioni eq 3 then strioni = 'III'
IF topxioni eq 4 then strioni = 'IV'
IF topxioni eq 5 then strioni = 'V'
IF topxioni eq 6 then strioni = 'VI'
getabnd, xionnmb, xion[0], abnd, flag = 1
IF xion[1] eq 1 then strionib = 'I'
IF xion[1] eq 2 then strionib = 'II'
IF xion[1] eq 3 then strionib = 'III'
IF xion[1] eq 4 then strionib = 'IV'
IF xion[1] eq 5 then strionib = 'V'
IF xion[1] eq 6 then strionib = 'VI'

xlabel = 'Log N!d'+xionnm+strioni+'!n/N!d'+xionnmb+strionib+'!n'




;;;; flagging is a little weird. 

;; Plot the detections. Don't just plot everything, since we
;; want to only outline non-detections

blendsredx = where((megastruct[red].ion[topxionz,topxioni].nflg eq 3 AND megastruct[red].ion[topxionz,topxioni].trans[0].flg eq 3) OR $
                 (megastruct[red].ion[topxionz,topxioni].nflg eq 3 AND megastruct[red].ion[topxionz,topxioni].trans[1].flg eq 3) OR $
                 (megastruct[red].ion[topxionz,topxioni].nflg eq 3 AND megastruct[red].ion[topxionz,topxioni].trans[2].flg eq 3), nblendredx)
blendsbluex =  where((megastruct[blue].ion[topxionz,topxioni].nflg eq 3 AND megastruct[blue].ion[topxionz,topxioni].trans[0].flg eq 3) OR $
                 (megastruct[blue].ion[topxionz,topxioni].nflg eq 3 AND megastruct[blue].ion[topxionz,topxioni].trans[1].flg eq 3) OR $
                 (megastruct[blue].ion[topxionz,topxioni].nflg eq 3 AND megastruct[blue].ion[topxionz,topxioni].trans[2].flg eq 3), nblendbluex)

blendsredratx = where((megastruct[red].ion[xion[0], xion[1]].nflg eq 3 AND $
                       megastruct[red].ion[xion[0], xion[1]].trans[0].flg eq 3) OR $
                 (megastruct[red].ion[xion[0], xion[1]].nflg eq 3 AND megastruct[red].ion[xion[0], xion[1]].trans[1].flg eq 3) OR $
                 (megastruct[red].ion[xion[0], xion[1]].nflg eq 3 AND megastruct[red].ion[xion[0], xion[1]].trans[2].flg eq 3), $
                      nblendredratx)
blendsblueratx =  where((megastruct[blue].ion[xion[0], xion[1]].nflg eq 3 AND $
                         megastruct[blue].ion[xion[0], xion[1]].trans[0].flg eq 3) OR $
                 (megastruct[blue].ion[xion[0], xion[1]].nflg eq 3 AND megastruct[blue].ion[xion[0], xion[1]].trans[1].flg eq 3) OR $
                 (megastruct[blue].ion[xion[0], xion[1]].nflg eq 3 AND megastruct[blue].ion[xion[0], xion[1]].trans[2].flg eq 3), $
                        nblendblueratx)


detx = where((megastruct[red].ion[topxionz, topxioni].nflg eq 1  AND megastruct[red].ion[xion[0], xion[1]].nflg eq 1) $
            OR (megastruct[red].ion[topxionz, topxioni].nflg eq 1  AND megastruct[red].ion[xion[0], xion[1]].nflg eq 2) $
            OR  (megastruct[red].ion[topxionz, topxioni].nflg eq 2  AND megastruct[red].ion[xion[0], xion[1]].nflg eq 1) $
            OR (megastruct[red].ion[topxionz, topxioni].nflg eq 2  AND megastruct[red].ion[xion[0], xion[1]].nflg eq 2), ndetx)

det2x = where((megastruct[blue].ion[topxionz, topxioni].nflg eq 1  AND megastruct[blue].ion[xion[0], xion[1]].nflg eq 1) $
            OR (megastruct[blue].ion[topxionz, topxioni].nflg eq 1  AND megastruct[blue].ion[xion[0], xion[1]].nflg eq 2) $
            OR  (megastruct[blue].ion[topxionz, topxioni].nflg eq 2  AND megastruct[blue].ion[xion[0], xion[1]].nflg eq 1) $
            OR (megastruct[blue].ion[topxionz, topxioni].nflg eq 2  AND megastruct[blue].ion[xion[0], xion[1]].nflg eq 2), ndet2x)

;;Add arrows for the Upper limits (non-detections) AND Lower limits
;;(saturated)

;;;;lower limits
xlowlim = where((megastruct.ion[topxionz, topxioni].nflg eq 2  AND megastruct.ion[xion[0], xion[1]].nflg eq 1) $
            OR  (megastruct.ion[topxionz, topxioni].nflg eq 1  AND megastruct.ion[xion[0], xion[1]].nflg eq 3), nxlow)
           

;;upper limits
xuplim = where((megastruct.ion[topxionz, topxioni].nflg eq 3  AND megastruct.ion[xion[0], xion[1]].nflg eq 1) $
            OR (megastruct.ion[topxionz, topxioni].nflg eq 1  AND megastruct.ion[xion[0], xion[1]].nflg eq 2) $
            OR  (megastruct.ion[topxionz, topxioni].nflg eq 3  AND megastruct.ion[xion[0], xion[1]].nflg eq 2), nxup)

;;; ambiguous weirdo limits


;;; things that are upper and lower limits at the same time ;;;
;;; aka lower limit over lower limit
;;; or non-detection over non-detection (should be clear, since they
;;; are open points. 

;;;;;;;;---------->>>> CUT THIS PART> 
;xlim3 = where((megastruct.ion[topxionz, topxioni].nflg eq 2  AND megastruct.ion[xion[0], xion[1]].nflg eq 2) $
  ;        OR (megastruct.ion[topxionz, topxioni].nflg eq 3  AND megastruct.ion[xion[0], xion[1]].nflg eq 3), nxlim3)

xlim3 = where((megastruct.ion[topxionz, topxioni].nflg eq 2  AND megastruct.ion[xion[0], xion[1]].nflg eq 2) , nxlim3)

;; The open symbols get plotted if any of the quantities is an upper
;; limit

;limplotx = where((megastruct.ion[topxionz, topxioni].nflg eq 3  AND megastruct.ion[xion[0], xion[1]].nflg eq 3) $
 ;             OR  (megastruct.ion[topxionz, topxioni].nflg eq 1  AND megastruct.ion[xion[0], xion[1]].nflg eq 3) $
    ;            OR (megastruct.ion[topxionz, topxioni].nflg eq 2  AND megastruct.ion[xion[0], xion[1]].nflg eq 3) $
      ;          OR (megastruct.ion[topxionz, topxioni].nflg eq 3  AND megastruct.ion[xion[0], xion[1]].nflg eq 1) $
       ;         OR  (megastruct.ion[topxionz, topxioni].nflg eq 3  AND megastruct.ion[xion[0], xion[1]].nflg eq 2),nlimplotx)      

limplotx = where( (megastruct.ion[topxionz, topxioni].nflg eq 1  AND megastruct.ion[xion[0], xion[1]].nflg eq 3) $
                OR (megastruct.ion[topxionz, topxioni].nflg eq 3  AND megastruct.ion[xion[0], xion[1]].nflg eq 1) $
                OR  (megastruct.ion[topxionz, topxioni].nflg eq 3  AND megastruct.ion[xion[0], xion[1]].nflg eq 2),nlimplotx)      

ENDIF



IF (galprop NE 'nion' AND galprop NE 'nionratio') THEN BEGIN

print, 'You have set xion without specifying nion or nionratio, which is wrong.'
STOP
ENDIF


ENDELSE






;;; Y-axis is the Ion Column or Ratio (Always);;;

IF ~keyword_set(ionratio) then begin

;; here we are plotting straight up column densities. 

yquant = megastruct.ion[ionz,ioni].lognion
getabnd, ionnm, ionz, abnd, flag = 1
IF ioni eq 1 then strioni = 'I'
IF ioni eq 2 then strioni = 'II'
IF ioni eq 3 then strioni = 'III'
IF ioni eq 4 then strioni = 'IV'
IF ioni eq 5 then strioni = 'V'
IF ioni eq 6 then strioni = 'VI'
ylabel = 'Log N!d'+ionnm+strioni+'!n [cm!u-2!n]'
goodcol = where(yquant gt 0.)
ymin = min(yquant[goodcol]) - 0.5
ymax = max(yquant[goodcol]) + 0.5
IF NOT keyword_set(inyrange) THEN BEGIN
yrange = [ymin, ymax]
ENDIF ELSE BEGIN
yrange = inyrange
ENDELSE


;; Plot the detections. Don't just plot everything, since we
;; want to only outline non-detections
;; but also account for the "detected blends" that are flagged as
;; upper limits. 
blendsred = where((megastruct[red].ion[ionz,ioni].nflg eq 3 AND megastruct[red].ion[ionz,ioni].trans[0].flg eq 3) OR $
                 (megastruct[red].ion[ionz,ioni].nflg eq 3 AND megastruct[red].ion[ionz,ioni].trans[1].flg eq 3) OR $
                 (megastruct[red].ion[ionz,ioni].nflg eq 3 AND megastruct[red].ion[ionz,ioni].trans[2].flg eq 3), nblendred)
blendsblue =  where((megastruct[blue].ion[ionz,ioni].nflg eq 3 AND megastruct[blue].ion[ionz,ioni].trans[0].flg eq 3) OR $
                 (megastruct[blue].ion[ionz,ioni].nflg eq 3 AND megastruct[blue].ion[ionz,ioni].trans[1].flg eq 3) OR $
                 (megastruct[blue].ion[ionz,ioni].nflg eq 3 AND megastruct[blue].ion[ionz,ioni].trans[2].flg eq 3), nblendblue)

;;;; these have to be defined. Are only valid paramters for the ratio
;;;; option, however. 

blendsredrat = -1
blendsbluerat = -1
nblendredrat = 0
nblendbluerat = 0


det = where(megastruct[red].ion[ionz, ioni].nflg eq 1 or megastruct[red].ion[ionz, ioni].nflg eq 2, ndet)
det2 = where(megastruct[blue].ion[ionz, ioni].nflg eq 1 or megastruct[blue].ion[ionz, ioni].nflg eq 2, ndet2)

;;Add arrows for the Upper limits (non-detections)
lim = where(megastruct.ion[ionz, ioni].nflg eq 3, nlim)
lim2 = where(megastruct.ion[ionz, ioni].nflg eq 2, nlim2)
lim3 = -1
nlim3 = 0 ; for consistency with stupid where. 
limplot = -1 
nlimplot = 0

endif ELSE BEGIN

;;; here we are plotting an ionic ratio of some kind. 

iontop = megastruct.ion[ionz,ioni].lognion
ionbot =  megastruct.ion[ionratio[0],ionratio[1]].lognion
yquant = iontop - ionbot
;goodcol = where(iontop NE 0. AND ionbot NE 0. AND  megastruct.ion[ionratio[0],ionratio[1]].nflg NE 3 $
           ;    AND megastruct.galaxy.sfr_uplim EQ 'no')

goodcol = where(iontop NE 0. AND ionbot NE 0. AND  megastruct.ion[8,6].nflg NE 3)
ymin = min(yquant[goodcol]) - 0.2
ymax = max(yquant[goodcol]) + 0.2
IF NOT keyword_set(inyrange) THEN BEGIN
yrange = [ymin, ymax]
ENDIF ELSE BEGIN
yrange = inyrange
ENDELSE

;;; labels ;;;

getabnd, ionnm, ionz, abnd, flag = 1
IF ioni eq 1 then strioni = 'I'
IF ioni eq 2 then strioni = 'II'
IF ioni eq 3 then strioni = 'III'
IF ioni eq 4 then strioni = 'IV'
IF ioni eq 5 then strioni = 'V'
IF ioni eq 6 then strioni = 'VI'
getabnd, ionnmb, ionratio[0], abnd, flag = 1
IF ionratio[1] eq 1 then strionib = 'I'
IF ionratio[1] eq 2 then strionib = 'II'
IF ionratio[1] eq 3 then strionib = 'III'
IF ionratio[1] eq 4 then strionib = 'IV'
IF ionratio[1] eq 5 then strionib = 'V'
IF ionratio[1] eq 6 then strionib = 'VI'

ylabel = 'Log N!d'+ionnm+strioni+'!n/N!d'+ionnmb+strionib+'!n'




;;;; flagging is a little weird. 

;; Plot the detections. Don't just plot everything, since we
;; want to only outline non-detections

blendsred = where((megastruct[red].ion[ionz,ioni].nflg eq 3 AND megastruct[red].ion[ionz,ioni].trans[0].flg eq 3) OR $
                 (megastruct[red].ion[ionz,ioni].nflg eq 3 AND megastruct[red].ion[ionz,ioni].trans[1].flg eq 3) OR $
                 (megastruct[red].ion[ionz,ioni].nflg eq 3 AND megastruct[red].ion[ionz,ioni].trans[2].flg eq 3), nblendred)
blendsblue =  where((megastruct[blue].ion[ionz,ioni].nflg eq 3 AND megastruct[blue].ion[ionz,ioni].trans[0].flg eq 3) OR $
                 (megastruct[blue].ion[ionz,ioni].nflg eq 3 AND megastruct[blue].ion[ionz,ioni].trans[1].flg eq 3) OR $
                 (megastruct[blue].ion[ionz,ioni].nflg eq 3 AND megastruct[blue].ion[ionz,ioni].trans[2].flg eq 3), nblendblue)

blendsredrat = where((megastruct[red].ion[ionratio[0], ionratio[1]].nflg eq 3 AND megastruct[red].ion[ionratio[0], ionratio[1]].trans[0].flg eq 3) OR $
                 (megastruct[red].ion[ionratio[0], ionratio[1]].nflg eq 3 AND megastruct[red].ion[ionratio[0], ionratio[1]].trans[1].flg eq 3) OR $
                 (megastruct[red].ion[ionratio[0], ionratio[1]].nflg eq 3 AND megastruct[red].ion[ionratio[0], ionratio[1]].trans[2].flg eq 3), nblendredrat)
blendsbluerat =  where((megastruct[blue].ion[ionratio[0], ionratio[1]].nflg eq 3 AND megastruct[blue].ion[ionratio[0], ionratio[1]].trans[0].flg eq 3) OR $
                 (megastruct[blue].ion[ionratio[0], ionratio[1]].nflg eq 3 AND megastruct[blue].ion[ionratio[0], ionratio[1]].trans[1].flg eq 3) OR $
                 (megastruct[blue].ion[ionratio[0], ionratio[1]].nflg eq 3 AND megastruct[blue].ion[ionratio[0], ionratio[1]].trans[2].flg eq 3), nblendbluerat)


det = where((megastruct[red].ion[ionz, ioni].nflg eq 1  AND megastruct[red].ion[ionratio[0], ionratio[1]].nflg eq 1) $
            OR (megastruct[red].ion[ionz, ioni].nflg eq 1  AND megastruct[red].ion[ionratio[0], ionratio[1]].nflg eq 2) $
            OR  (megastruct[red].ion[ionz, ioni].nflg eq 2  AND megastruct[red].ion[ionratio[0], ionratio[1]].nflg eq 1) $
            OR (megastruct[red].ion[ionz, ioni].nflg eq 2  AND megastruct[red].ion[ionratio[0], ionratio[1]].nflg eq 2), ndet)

det2 = where((megastruct[blue].ion[ionz, ioni].nflg eq 1  AND megastruct[blue].ion[ionratio[0], ionratio[1]].nflg eq 1) $
            OR (megastruct[blue].ion[ionz, ioni].nflg eq 1  AND megastruct[blue].ion[ionratio[0], ionratio[1]].nflg eq 2) $
            OR  (megastruct[blue].ion[ionz, ioni].nflg eq 2  AND megastruct[blue].ion[ionratio[0], ionratio[1]].nflg eq 1) $
            OR (megastruct[blue].ion[ionz, ioni].nflg eq 2  AND megastruct[blue].ion[ionratio[0], ionratio[1]].nflg eq 2), ndet2)

;;Add arrows for the Upper limits (non-detections) AND Lower limits
;;(saturated)

;;;;lower limits
lim2 = where((megastruct.ion[ionz, ioni].nflg eq 2  AND megastruct.ion[ionratio[0], ionratio[1]].nflg eq 1) $
            OR  (megastruct.ion[ionz, ioni].nflg eq 1  AND megastruct.ion[ionratio[0], ionratio[1]].nflg eq 3), nlim2)
           

;;upper limits
lim = where((megastruct.ion[ionz, ioni].nflg eq 3  AND megastruct.ion[ionratio[0], ionratio[1]].nflg eq 1) $
            OR (megastruct.ion[ionz, ioni].nflg eq 1  AND megastruct.ion[ionratio[0], ionratio[1]].nflg eq 2) $
            OR  (megastruct.ion[ionz, ioni].nflg eq 3  AND megastruct.ion[ionratio[0], ionratio[1]].nflg eq 2), nlim)

;;; ambiguous weirdo limits


;;; things that are upper and lower limits at the same time ;;;
;;; aka lower limit over lower limit
;;; or non-detection over non-detection (should be clear, since they
;;; are open points. 

;lim3 = where((megastruct.ion[ionz, ioni].nflg eq 2  AND megastruct.ion[ionratio[0], ionratio[1]].nflg eq 2) $
 ;           OR (megastruct.ion[ionz, ioni].nflg eq 3  AND megastruct.ion[ionratio[0], ionratio[1]].nflg eq 3), nlim3)

lim3 = where((megastruct.ion[ionz, ioni].nflg eq 2  AND megastruct.ion[ionratio[0], ionratio[1]].nflg eq 2), nlim3)

;; The open symbols get plotted if any of the quantities is a lower
;; limit

;limplot = where((megastruct.ion[ionz, ioni].nflg eq 3  AND megastruct.ion[ionratio[0], ionratio[1]].nflg eq 3) $
  ;              OR  (megastruct.ion[ionz, ioni].nflg eq 1  AND megastruct.ion[ionratio[0], ionratio[1]].nflg eq 3) $
   ;             OR (megastruct.ion[ionz, ioni].nflg eq 2  AND megastruct.ion[ionratio[0], ionratio[1]].nflg eq 3) $
    ;            OR (megastruct.ion[ionz, ioni].nflg eq 3  AND megastruct.ion[ionratio[0], ionratio[1]].nflg eq 1) $
    ;            OR  (megastruct.ion[ionz, ioni].nflg eq 3  AND megastruct.ion[ionratio[0], ionratio[1]].nflg eq 2),nlimplot)      

limplot = where((megastruct.ion[ionz, ioni].nflg eq 1  AND megastruct.ion[ionratio[0], ionratio[1]].nflg eq 3) $
                OR (megastruct.ion[ionz, ioni].nflg eq 3  AND megastruct.ion[ionratio[0], ionratio[1]].nflg eq 1) $
                OR  (megastruct.ion[ionz, ioni].nflg eq 3  AND megastruct.ion[ionratio[0], ionratio[1]].nflg eq 2),nlimplot)      

ENDELSE

;;;;;;;;;Anyion, just for Y axis at this point. 

;;; currently not an active area Will stop. May return to this...

IF keyword_set(anyion) then begin
print, 'WARNING: GO BACK! THIS OPTION IS NOT SUPPORTED!'
STOP
;; we must not have ionratio flagged. Make sure of this. 

IF keyword_set(ionratio) then begin
print, 'You cannot have anyion and ionratio both set here. It makes no sense.'
STOP 
ENDIF


;;; first, we must re-define ionz and ioni. 

nsys = n_elements(megastruct.galaxy.galid)
dummyionz = indgen(nsys) * 0
dummyioni = indgen(nsys) * 0

FOR izz=0, nsys -1 do begin
goodlines = where((megastruct[ii].ion.trans.flg EQ 1.0 OR $
                  megastruct[ii].ion.trans.flg EQ 3.0 OR $
                  megastruct[ii].ion.trans.flg EQ 9.0 OR $
                 megastruct[ii].ion.trans.flg EQ 11.0) AND $
                 megastruct[ii].ion.trans.name NE 'OVI' AND $
                  megastruct[ii].ion.trans.name NE 'HI' AND $
                  megastruct[ii].ion.trans.name NE 'PV' AND $
                 megastruct[ii].ion.trans.name NE 'PII' AND $
                 megastruct[ii].ion.trans.name NE 'AlII' AND $
                 megastruct[ii].ion.trans.name NE 'SIV' AND $
                 megastruct[ii].ion.trans.name NE 'SII' AND $
                  megastruct[ii].ion.trans.name NE 'SVI' AND $
                 megastruct[ii].ion.trans.name NE 'NV') ;;; AAAAAAHHHHHH!!!!

idets = where((megastruct[ii].ion.trans.flg)[goodlines] EQ 1.0, nidets)
iblends =  where((megastruct[ii].ion.trans.flg)[goodlines] EQ 3.0, niblends)
isat =  where((megastruct[ii].ion.trans.flg)[goodlines] EQ 9.0, nisat)
isatblend = where((megastruct[ii].ion.trans.flg)[goodlines] EQ 11.0, nisatblend)

;; give preference to the detections, then the minor blends, then
;; saturated, then minor saturated. 

IF (nidets EQ 1 AND idets[0] ge 0) then begin
;; the ion structure is [27,16] in size. this is immutable. hard-wired
;; in is bad form, but I don't care! 

dummyionz[izz] = (goodlines[idets]/16) MOD 27
dummyioni[izz] = (goodlines[idets]/16) / 27

ENDIF

IF (nidets GT 1 AND idets[0] ge 0) then begin
;; the ion structure is [27,16] in size. this is immutable. hard-wired
;; in is bad form, but I don't care! 

;; in the case of several good detections, we use the one with the
;; minimum column. 

mincol = min(megastruct[ii].ion[( (goodlines[idets]/16) MOD 27), ((goodlines[idets]/16) / 27) ].lognion, mini)
dummyionz[izz] = (goodlines[idets[mini]]/16) MOD 27
dummyioni[izz] = (goodlines[idets[mini]]/16) / 27

ENDIF



endfor

;; here we are plotting straight up column densities. 

yquant = megastruct.ion[ionz,ioni].lognion
getabnd, ionnm, ionz, abnd, flag = 1
IF ioni eq 1 then strioni = 'I'
IF ioni eq 2 then strioni = 'II'
IF ioni eq 3 then strioni = 'III'
IF ioni eq 4 then strioni = 'IV'
IF ioni eq 5 then strioni = 'V'
IF ioni eq 6 then strioni = 'VI'
ylabel = 'Log N!d'+ionnm+strioni+'!n [cm!u-2!n]'
goodcol = where(yquant gt 0.)
ymin = min(yquant[goodcol]) - 0.5
ymax = max(yquant[goodcol]) + 0.5
yrange = [ymin, ymax]

;; Plot the detections. Don't just plot everything, since we
;; want to only outline non-detections
;; but also account for the "detected blends" that are flagged as
;; upper limits. 
blendsred = where((megastruct[red].ion[ionz,ioni].nflg eq 3 AND megastruct[red].ion[ionz,ioni].trans[0].flg eq 3) OR $
                 (megastruct[red].ion[ionz,ioni].nflg eq 3 AND megastruct[red].ion[ionz,ioni].trans[1].flg eq 3) OR $
                 (megastruct[red].ion[ionz,ioni].nflg eq 3 AND megastruct[red].ion[ionz,ioni].trans[2].flg eq 3), nblendred)
blendsblue =  where((megastruct[blue].ion[ionz,ioni].nflg eq 3 AND megastruct[blue].ion[ionz,ioni].trans[0].flg eq 3) OR $
                 (megastruct[blue].ion[ionz,ioni].nflg eq 3 AND megastruct[blue].ion[ionz,ioni].trans[1].flg eq 3) OR $
                 (megastruct[blue].ion[ionz,ioni].nflg eq 3 AND megastruct[blue].ion[ionz,ioni].trans[2].flg eq 3), nblendblue)

;;;; these have to be defined. Are only valid paramters for the ratio
;;;; option, however. 

blendsredrat = -1
blendsbluerat = -1
nblendredrat = 0
nblendbluerat = 0


det = where(megastruct[red].ion[ionz, ioni].nflg eq 1 or megastruct[red].ion[ionz, ioni].nflg eq 2, ndet)
det2 = where(megastruct[blue].ion[ionz, ioni].nflg eq 1 or megastruct[blue].ion[ionz, ioni].nflg eq 2, ndet2)

;;Add arrows for the Upper limits (non-detections)
lim = where(megastruct.ion[ionz, ioni].nflg eq 3, nlim)
lim2 = where(megastruct.ion[ionz, ioni].nflg eq 2, nlim2)
lim3 = -1
nlim3 = 0 ; for consistency with stupid where. 
limplot = -1 
nlimplot = 0

endif 

;;;;; Now, Y-axis is all moot if we just want an EW
IF keyword_set(eqw) THEN BEGIN
;;; As of right now, this is taking the lambda of the most common ion
;;; of a given species. 
ngals  = n_elements(megastruct.galaxy.galid)
thistrans = indgen(ngals)
thisblend = sindgen(ngals)
thislam = 1750. ; should give null results if we're not in a supported transition
IF (ionz eq 14 AND ioni eq 2) THEN BEGIN
thislam = 1260.
ENDIF
IF (ionz eq 6 AND ioni eq 2) THEN BEGIN
thislam = 1036.
;thislam = 1334.
ENDIF
IF (ionz eq 7 AND ioni eq 2) THEN BEGIN
thislam = 1083.
ENDIF
IF (ionz eq 12 AND ioni eq 2) THEN BEGIN
thislam = 2796.
ENDIF
IF (ionz eq 14 AND ioni eq 3) THEN BEGIN
thislam = 1206.
ENDIF
IF (ionz eq 6 AND ioni eq 3) THEN BEGIN
thislam = 977.
ENDIF
IF (ionz eq 7 AND ioni eq 3) THEN BEGIN
thislam = 989.
ENDIF
IF (ionz eq 14 AND ioni eq 4) THEN BEGIN
thislam = 1393.
ENDIF
IF (ionz eq 6 AND ioni eq 4) THEN BEGIN
thislam = 1548.
ENDIF
IF (ionz eq 8 AND ioni eq 6) THEN BEGIN
thislam = 1031.
ENDIF
IF (ionz eq 1 AND ioni eq 1) THEN BEGIN
thislam = 1215.
ENDIF

FOR i = 0, ngals - 1 do begin
thisblend[i] = 'noblend'
thistrans[i] = where(abs(megastruct[i].ion[ionz, ioni].trans.lambda - thislam) lt 1.0)
IF (thistrans[i] lt 0) THEN BEGIN
thistrans[i] = 11 ; the values here are zeroed out for all known ions, I think. 
ENDIF
IF ((megastruct[i].ion[ionz, ioni].trans[thistrans[i]].flg EQ 5.) OR $
     (megastruct[i].ion[ionz, ioni].trans[thistrans[i]].flg EQ 7.)) THEN BEGIN
yquant[i] = 2.0 * megastruct[i].ion[ionz,ioni].trans[thistrans[i]].sigwrest ; 2 sigma upper limits
ENDIF ELSE BEGIN

yquant[i] = megastruct[i].ion[ionz,ioni].trans[thistrans[i]].wrest ;will give infinity for many undefined. sigh.
ENDELSE
IF ((megastruct[i].ion[ionz, ioni].trans[thistrans[i]].flg EQ 3.)) THEN BEGIN
thisblend[i] = 'blend'
ENDIF

endfor

 

; In Angstrom
yquant = yquant / 1000.
ylabel = 'EW !drest,'+ionnm+strioni+','+strtrim(string(thislam, format = '(i4)'), 2)+'!n [A]'
goodcol = where(yquant gt 0.)
ymin = 0.0060 ;min(yquant[goodcol])

ymax = max(yquant[goodcol]) + 1.0
yrange = [ymin, ymax]


;; Plot the detections. Don't just plot everything, since we
;; want to only outline non-detections
;; but also account for the "detected blends" that are flagged as
;; upper limits. 
blendsred = -1
blendsblue = -1

;;;; these have to be defined. Are only valid paramters for the ratio
;;;; option, however. 

blendsredrat = -1
blendsbluerat = -1
nblendredrat = 0
nblendbluerat = 0


det = where(megastruct[red].ion[ionz, ioni].nflg eq 1 or megastruct[red].ion[ionz, ioni].nflg eq 2, ndet)
det2 = where(megastruct[blue].ion[ionz, ioni].nflg eq 1 or megastruct[blue].ion[ionz, ioni].nflg eq 2, ndet2)

;;Add arrows for the Upper limits (non-detections)
lim = where((megastruct.ion[ionz, ioni].nflg eq 3 AND thisblend NE 'blend'), nlim) ; will have not blends in it. whatever. I think this is proper. 
lim2 = -1 ; irrelevant
lim3 = -1
nlim3 = 0 ; for consistency with stupid where. 
limplot = -1 
nlimplot = 0
ENDIF

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                    ;;
;;              PLOTTING HERE         ;;
;;                                    ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;; PLOT ALL GALAXIES ;;;;;;;;;;;;;;;;
if not keyword_set(xdev) then begin ; if we have not been told to use the x device 
  if not keyword_set(XPS) then ct_psopen, psfile, xs = 10, ys = 8 else $
        	x_psopen, psfile, /maxs
endif 
;ct_load_colours
 clr = getcolor(/load)
  device, /encapsulated

 thickratio = 1.
 if (!d.name eq 'X') then thickratio = 4.

;;Set up colour indices. These correspond to colours defined by jw_setclrs
if keyword_set(white) then begin
    xclr = jw_setclrs(/white)
   bluecol = xclr[21]
   redcol = xclr[22]
   fgcol = 0
   bgcol = 255
endif else begin
   xclr = jw_setclrs(/dark)
   bluecol = xclr[21]
   redcol = xclr[22]
   fgcol = 255
   bgcol = 0
endelse 

ct_psfill_black, col = bgcol

IF ~keyword_set(EQW) THEN BEGIN
plot, [0], [0], /nodata, yr = yrange, ys = 1, xr = xrange, xs = 1, $
      ytit = ylabel, xtit = xlabel, /noerase, col = fgcol, $
      chars = 2.0, _extra=_extra 
ENDIF ELSE BEGIN
plot, [0], [0], /nodata, yr = yrange, ys = 1, xr = xrange, xs = 1, $
      ytit = ylabel, xtit = xlabel, /noerase, col = fgcol, $
      chars = 2.0, /ylog, _extra=_extra 
ENDELSE

;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; DATA

;; Plot the detections. Don't just plot everything, since we
;; want to only outline non-detections

oplot, xquant[goodcol], yquant[goodcol], psym = symcat(16, col = cgcolor('navy')), syms = 1.5


 

IF (lim[0] GE 0) THEN BEGIN
for ii = 0, nlim-1 do begin
   ;;First plot a blue outlined square or red-outlined triangle as the
   ;;symbol
  ; IF ~keyword_set(ionratio) THEN BEGIN
 ;  if megastruct[lim[ii]].galaxy.sfr_uplim eq "yes" then begin
      ;;Red
 ;     oplot, [xquant[lim[ii]]], [yquant[lim[ii]]], $
  ;           psym = 4, col = redcol, syms = 3.0, thick = 8 / thickratio 
      
  ; endif else begin
  ;    ;;blue
  ;    oplot, [xquant[lim[ii]]], [yquant[lim[ii]]], $
   ;          psym = 6, col = bluecol, syms = 2.3, thick = 8 / thickratio 
   ;  
 ;  endelse 
;ENDIF

   ;;Now plot the down-facing arrow
   plotsym, 1, 2.5, thick = 10 / thickratio 
   IF (megastruct[lim[ii]].galaxy.sfr_uplim EQ 'no') THEN BEGIN
   oplot, [xquant[lim[ii]]], [yquant[lim[ii]]], psym = 8, col = fgcol
ENDIF

endfor 
;;Add arrows for the lower limits (saturated)
ENDIF

IF (lim2[0] GE 0) THEN BEGIN
for ii = 0, nlim2-1 do begin
   plotsym, 2, 2, thick = 10 / thickratio 
   IF (megastruct[lim2[ii]].galaxy.sfr_uplim EQ 'no') THEN BEGIN
   oplot, [xquant[lim2[ii]]], [yquant[lim2[ii]]], psym = 8, col = fgcol
ENDIF

endfor 
ENDIF

;;;; Galaxy Property Limits
IF (xuplim[0] GE 0) THEN BEGIN
for ii = 0, nxup-1 do begin
   plotsym, 6, 2.5, thick = 10 / thickratio 
   IF (megastruct[xuplim[ii]].galaxy.sfr_uplim EQ 'no') THEN BEGIN
   oplot, [xquant[xuplim[ii]]], [yquant[xuplim[ii]]], psym = 8, col = fgcol
   ENDIF
endfor 
ENDIF

IF (xlowlim[0] GE 0) THEN BEGIN
for ii = 0, nxlow-1 do begin
   plotsym, 7, 2, thick = 10 / thickratio 
   IF (megastruct[xlowlim[ii]].galaxy.sfr_uplim EQ 'no') THEN BEGIN
   oplot, [xquant[xlowlim[ii]]], [yquant[xlowlim[ii]]], psym = 8, col = fgcol
   ENDIF
endfor 
ENDIF


;;; weirdos. upper and lower limits. 
IF (lim3[0] GE 0) THEN BEGIN
print, 'AAAAH'
for ii = 0, nlim3-1 do begin
   plotsym, 2, 2, thick = 10 / thickratio 
   IF (megastruct[lim3[ii]].galaxy.sfr_uplim EQ 'no') THEN BEGIN
   oplot, [xquant[lim3[ii]]], [yquant[lim3[ii]]], psym = 8, col = fgcol

    plotsym, 1, 2, thick = 10 / thickratio 
  
   oplot, [xquant[lim3[ii]]], [yquant[lim3[ii]]], psym = 8, col = fgcol
ENDIF

endfor 
ENDIF

;;; weirdos. upper and lower limits. 
IF (xlim3[0] GE 0) THEN BEGIN
print, 'AAAAHX'
for ii = 0, nxlim3-1 do begin
   IF (megastruct[xlim3[ii]].galaxy.sfr_uplim EQ 'no') THEN BEGIN
   plotsym, 6, 2, thick = 10 / thickratio 
   oplot, [xquant[xlim3[ii]]], [yquant[xlim3[ii]]], psym = 8, col = fgcol
  plotsym, 7, 2, thick = 10 / thickratio 
   oplot, [xquant[xlim3[ii]]], [yquant[xlim3[ii]]], psym = 8, col = fgcol
ENDIF

endfor 
ENDIF

;oplot, xquant[goodcol], yquant[goodcol], psym = symcat(16, col = cgcolor('cornflower blue')), syms = 1.1
oplot, xquant[goodcol], yquant[goodcol], psym = symcat(16, col = cgcolor('red')), syms = 1.1
oplot, xquant[goodval], yquant[goodval], psym = symcat(16, col = cgcolor('cornflower blue')), syms = 1.1
oplot, xquant[goodcol], yquant[goodcol], psym = symcat(7, col = cgcolor('black')), syms = 0.7

;; label the model lines

;xyouts, -0.2, -0.25, 'Photoionization, Log U', orientation = 55.0, charsize = 0.8
;xyouts, -0.35, -2.3, 'Collisional Ionization, T', orientation = 27.0, charsize = 0.8
if (keyword_set(tagid)) then begin 
   plotsym, 0, 3, thick=3 
   oplot, [xquant[tagid]], [yquant[tagid]], psym = 8, col =fgcol   
endif  


;;Panel label for HST proposal
;xyouts, 0.20, 0.87, 'D', chars = 3.0, /normal, font = 1

;ct_psclose, /noslug
if not keyword_set(xdev) then begin 
  if not keyword_set(XPS) then ct_psclose, /noslug else x_psclose
   spawn, 'idlepstopdf.sh '+psfile
endif 

;popd
end
