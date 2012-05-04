
pro coshalo_lowions, ionz, ioni, galprop, psfile = psfile, white = white,  anyion = anyion, eqw = eqw, xion = xion, ionratio = ionratio, MEGASTRUCT=megastruct

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


;;;;;;; SET PLOT NAME ;;;;;;;;
;;;;;;; THE DEFAULTS SHOULD BE SENSIBLE ;;;;;

 psfile = 'coshalo_lowions.ps'

;;;;;; END SET PLOT NAME ;;;;;;;;

ldir = getenv('DROPBOX_DIR')+'/COS-Halos/lowions/'
tdir = getenv('DROPBOX_DIR')+'/COS-Halos/Targets/'

if not keyword_set(MEGASTRUCT) then restore, ldir+'/cosmetals_megastructure.sav' 

if not keyword_set(LSZ) then lsz = 2.0
   yrange = [0.003, 3.]

;pushd, ldir
   


;;red/blue cut
red = where(megastruct.galaxy.sfr_uplim eq 'yes', nred, complement = blue, ncomplement = nblue)

;;;;;;;;;;;;;;;;;; GET XY RANGES, AXIS LABELS AND
;;;;;;;;;;;;;;;;;; QUANTITIES;;;;;;;;;;;;;;;;;;

;; First, get the quantities. 

   
;;; X-axis is the galaxy property ;;;
   galprop = 'rhokpc'

;;;;; galprop is a string of any one of the following: 'smass' for
;;;;; stellar mass (actually a log); 'rhokpc' for impact parameter in kpc; 'sfr' for
;;;;; the star formation rate in msun per yr; 'metal' for galaxy
;;;;; metallicity, this will be the 12 + Log (O/H) abundance; 

   ;; impact parameter
   xrange = [0.0, 200.0]
   xlabel = greek('rho', /FORCE_PS)+'!X [kpc]'
   xquant =  megastruct.rhofinal
   xuplim = -1
   xlowlim = -1
   nxup = 0
   nxlow = 0

;;;;; Now, Y-axis is all moot if we just want an EW

;Those are flagged in the transition structures, and are our bitwise flags :
;1 = use
;0 = not use
;
;2 = blend
;4 = non-detection
;8 = saturated
;
;Thus, 11 = 1 + 8 + 2 = saturated blend that we want to use
;And 9 = saturated line
;5 = non-detection
;1 = good.
;3 = usable blend
;etc., etc.
;
;In the case of EW, when the saturated flag is triggered, I do not call the line a lower limit, since EW is EW... but the flag is still there. 

   EQW = 1
   ngals  = n_elements(megastruct.galaxy.galid)
   thistrans = indgen(ngals)
   thisblend = sindgen(ngals)
   take_tran = fltarr(ngals)
   yquant = fltarr(ngals)
   flg_det = intarr(ngals)
   thislam = [1260.4, 1334.5, 2796.3]  ;; SiII, CII, MgII
   aionz = [14, 6, 12]
   ioni = 2
   nlam = n_elements(thislam)

   FOR i = 0, ngals - 1 do begin
      thisblend[i] = 'noblend'

      for jj=0L,nlam-1 do begin
         ionz = aionz[jj]
         flg_low = 0
         yquant[i] = 99.
         thistrans[i] = where(abs(megastruct[i].ion[ionz, ioni].trans.lambda - thislam[jj]) lt 0.2)

         IF (thistrans[i] lt 0) THEN BEGIN
            thistrans[i] = 11   ; the values here are zeroed out for all known ions, I think. 
         ENDIF

         IF ((megastruct[i].ion[ionz, ioni].trans[thistrans[i]].flg EQ 5.) OR $
             (megastruct[i].ion[ionz, ioni].trans[thistrans[i]].flg EQ 7.)) THEN BEGIN
            if flg_low LE 0 then begin
               yquant[i] = yquant[i] < (2.0 * megastruct[i].ion[ionz,ioni].trans[thistrans[i]].sigwrest) ; 2 sigma upper limits
               flg_det[i] = 3
            endif
         ENDIF ELSE BEGIN
            ;; Value
            if thistrans[i] NE 11 then begin
               if flg_low EQ 1 then begin
                  if megastruct[i].ion[ionz,ioni].trans[thistrans[i]].wrest  GT yquant[i] then take = 1 else take = 0
               endif else take = 1
               if TAKE EQ 1 then begin
                  yquant[i] = megastruct[i].ion[ionz,ioni].trans[thistrans[i]].wrest ;will give infinity for many undefined. sigh.
                  IF ((megastruct[i].ion[ionz, ioni].trans[thistrans[i]].flg EQ 3.)) THEN BEGIN
                     thisblend[i] = 'blend'
                  ENDIF
                  take_tran[i] = thislam[jj]
               endif
               print, i, jj, take_tran[i], xquant[i], yquant[i] 
               flg_low = 1
               flg_det[i] = 1
            endif
         ENDELSE

      endfor
      
   endfor

 

   ;; In Angstrom
   yquant = yquant / 1000.
   ;ylabel = 'EW !drest,'+ionnm+strioni+','+strtrim(string(thislam, format = '(i4)'), 2)+'!n [A]'
   ylabel = 'W!drest!N (Low Ions) [A]'
   goodcol = where(yquant gt 0.)
   ymin = 0.0060                ;min(yquant[goodcol])

   ymax = max(yquant[goodcol]) + 1.0

   
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
   
   
   det = where(flg_det[red] EQ 1, ndet) ;megastruct[red].ion[ionz, ioni].nflg eq 1 or megastruct[red].ion[ionz, ioni].nflg eq 2, ndet)
   det2 = where(flg_det[blue] EQ 1, ndet2) ;megastruct[blue].ion[ionz, ioni].nflg eq 1 or megastruct[blue].ion[ionz, ioni].nflg eq 2, ndet2)
   
;;Add arrows for the Upper limits (non-detections)
   lim = where(flg_det EQ 3 AND thisblend NE 'blend', nlim) ; will have not blends in it. whatever. I think this is proper. 

   lim2 = -1                                                                          ; irrelevant
   lim3 = -1
   nlim3 = 1                    ; for consistency with stupid where. 
   limplot = -1 
   nlimplot = 1

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                    ;;
;;              PLOTTING HERE         ;;
;;                                    ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;; PLOT ALL GALAXIES ;;;;;;;;;;;;;;;;
;ct_psopen, psfile, xs = 10, ys = 8
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


;; Loop !

for ss=0,2 do begin

   ct_psfill_black, col = bgcol
   plot, [0], [0], /nodata, yr = yrange, ys = 1, xr = xrange, xs = 1, $
         ytit = ylabel, xtit = xlabel, col = fgcol, $
         chars = 2.0, /ylog, xmarg=[8,2], ymarg=[4,1], /noeras
   
;; Plot the detections. Don't just plot everything, since we
;; want to only outline non-detections
   IF (det[0] GE 0) THEN BEGIN
      oplot, xquant[red[det]], yquant[red[det]], psym = symcat(14, col = fgcol), syms = 3.0
      oplot, xquant[red[det]], yquant[red[det]], psym = symcat(14, col = redcol), syms = 2.0
   ENDIF
   
   IF (det2[0] GE 0) THEN BEGIN
      oplot, xquant[blue[det2]], yquant[blue[det2]], psym = symcat(15, col = fgcol), syms = 2.3
      oplot, xquant[blue[det2]], yquant[blue[det2]], psym = symcat(15, col = bluecol), syms = 1.5
   ENDIF
   

   IF (lim[0] GE 0) THEN BEGIN
      for ii = 0, nlim-1 do begin
         ;;First plot a blue outlined square or red-outlined triangle as the
         ;;symbol
         IF ~keyword_set(ionratio) THEN BEGIN
            if megastruct[lim[ii]].galaxy.sfr_uplim eq "yes" then begin
               ;;Red
               oplot, [xquant[lim[ii]]], [yquant[lim[ii]]], $
                      psym = 4, col = redcol, syms = 3.0, thick = 8
            
            endif else begin
               ;;blue
               oplot, [xquant[lim[ii]]], [yquant[lim[ii]]], $
                      psym = 6, col = bluecol, syms = 2.3, thick = 8
            endelse 
         ENDIF
         
         ;;Now plot the down-facing arrow
         plotsym, 1, 4, thick = 10
         oplot, [xquant[lim[ii]]], [yquant[lim[ii]]], psym = 8, col = fgcol
      endfor 
;;Add arrows for the lower limits (saturated)
   ENDIF

   IF (limplot[0] GE 0) THEN BEGIN
      for ii = 0, nlimplot-1 do begin
         ;;First plot a blue outlined square or red-outlined triangle as the
         ;;symbol
         
         if megastruct[limplot[ii]].galaxy.sfr_uplim eq "yes" then begin
            ;;Red
            oplot, [xquant[limplot[ii]]], [yquant[limplot[ii]]], $
                   psym = 4, col = redcol, syms = 3.0, thick = 8
            
         endif else begin
            ;;blue
            oplot, [xquant[limplot[ii]]], [yquant[limplot[ii]]], $
                   psym = 6, col = bluecol, syms = 2.3, thick = 8
         endelse 
      endfor 
   ENDIF

   IF (lim2[0] GE 0) THEN BEGIN
      for ii = 0, nlim2-1 do begin
         plotsym, 2, 4, thick = 10
         oplot, [xquant[lim2[ii]]], [yquant[lim2[ii]]], psym = 8, col = fgcol
      endfor 
   ENDIF
   
;;;; Galaxy Property Limits
   IF (xuplim[0] GE 0) THEN BEGIN
      for ii = 0, nxup-1 do begin
         plotsym, 6, 4, thick = 10
         oplot, [xquant[xuplim[ii]]], [yquant[xuplim[ii]]], psym = 8, col = fgcol
      endfor 
   ENDIF
   
   IF (xlowlim[0] GE 0) THEN BEGIN
      for ii = 0, nxlow-1 do begin
         plotsym, 7, 4, thick = 10
         oplot, [xquant[xlowlim[ii]]], [yquant[xlowlim[ii]]], psym = 8, col = fgcol
      endfor 
   ENDIF

;;; weirdos. upper and lower limits. 
   IF (lim3[0] GE 0) THEN BEGIN
      for ii = 0, nlim3-1 do begin
         plotsym, 2, 4, thick = 10
         oplot, [xquant[lim3[ii]]], [yquant[lim3[ii]]], psym = 8, col = fgcol
         plotsym, 1, 4, thick = 10
         oplot, [xquant[lim3[ii]]], [yquant[lim3[ii]]], psym = 8, col = fgcol
      endfor 
   ENDIF
   
   
   ;; Covering fraction
   if ss GE 1 then begin
      rho_cut = 70.
      ;; Split them
      oplot, replicate(rho_cut, 2), [1e-10,1e10], color=clr.green, linesty=2, thick=7
      oplot, replicate(160., 2), [1e-10,1e10], color=clr.green, linesty=2, thick=7

      ;; Stats
      bin1 = where(xquant LT rho_cut, nbin1, complement=bin2, ncomplement=nbin2)
      
      for jj=0,1 do begin
         case jj of 
            0: begin
               binv = bin1
               xxy = rho_cut/2.
               nbin = nbin1
            end
            1: begin
               binv = bin2
               nbin = nbin2
               xxy = mean([rho_cut,160.])
            end
            else: stop
         endcase
         det_10 = where(yquant[binv] GT 0.1 AND flg_det[binv], npos)
         cov_10 = float(npos)/nbin
         det_03 = where(yquant[binv] GT 0.03 AND flg_det[binv], npos)
         cov_03 = float(npos)/nbin
         xyouts, xxy, yrange[1]/1.5, 'C!df!N = '+string(cov_10,format='(f4.2)')+$
                 ' ('+string(cov_03,format='(f4.2)')+')', color=clr.yellow, charsiz=lsz, $
                 align=0.5
      endfor
   endif

   ;; Chen+10b model (Using h_100 instead of h_72)
   if ss GT 1 then begin
      rho = findgen(200.) + 1
      a0 = 1.5
      a1 = -1.8
      EW = exp(a0 + a1 * alog10(rho * (72./100)))
      
      oplot, rho, EW, color=clr.magenta, linesty=2, thick=7
   endif

   plot, [0], [0], /nodata

endfor
;stop



;;Panel label for HST proposal
;xyouts, 0.20, 0.87, 'D', chars = 3.0, /normal, font = 1

;ct_psclose, /noslug
x_psclose
;spawn, 'idlepstopdf.sh '+psfile

;popd
end
