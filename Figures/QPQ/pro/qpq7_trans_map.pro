;+
;procedure che fa il ps di una simulation e aggiunge color bar
;range--> range of colorbar
;zoom is the the factor to zoom in a region 
;eps for encapsulated
;ps is the plate scale
;print save output in ps
;log for log scale
;zscale --> set zscale, if not set defualt if min max
;text allows label
;tlocation --> where to put the label
;colormap  allows to change the colormap
;grey     -->set to colormap 0
;reverse  reverse color map
;format  set the format for the colorbar ('F4.3')
;title -- display title
;coltitle --> title for the color bar
;colpos --> position for the color bar
;rotate --> 0-7 angle parameter for  rotation
;smooth --> the width of the smmothing function
;fits--> if set, im1 is passed as an array
;nocolbar  --> disable colorbar 
;x/ytitle   --> put a title in the bottom/left of the plot
;white     --> if set, the text is printed in black and background is
;              in white
;left      --> move the label to the left
;black     --> force black background
;circ      --> circle or radius pix, centered
;
;-
;; EW maps for QPQ6 HI Lya

; qpq7_trans_map, 1334.5323d, range=[0., 3], SIG_FINE=0.05
PRO qpq7_trans_map, wrest, RANGE=range, ZOOM=zoom, EPS=eps,$
                   PS=ps, print=print, log=log, zscale=zscale,  units=units, text=text, tlocation=tlocation,$
                   colormap=colormap,grey=grey,  reverse=reverse, format=format, title=title,$
                   coltitle=coltitle, rotate=rotate, PSFILE=psfile, $
                   smooth=smooth, fits=fits, colpos=colpos, nocolbar=nocolbar, xtitle=xtitle, ytitle=ytitle,$
                   white=white, left=left, black=black, _extra=extra, annotatecolor=annotatecolor, circ=circ, $
                 SIG_CLM=sig_clm, SIG_FINE=sig_fine


  ;; Initialize
  getion, wrest, ioni, elm, Z=ionz, NM=strioni
  transc = elm+strioni+strtrim(fix(wrest),2)
  if not keyword_set( PSFILE ) then psfile = 'qpq7_'+transc+'_map.ps'

  if not keyword_set(RANGE) then range = [0., 10]

  if not keyword_set(LOW_BOOST) then low_boost = 1.  ;; Boost saturated
  if not keyword_set(UPPER_BOOST) then upper_boost = 1.  ;; Drop lower limit

  if not keyword_set(SIG_CLM) then sig_clm = 0.05 ; Ang
  if not keyword_set(SIG_FINE) then sig_fine = 0.05 ; dex

  if not keyword_set(SAMP_LIM) then samp_lim = 6L  ;; Minimum number for stats

  ;; Dimensions
  grid_sz = 200. ; kpc
  if not keyword_set(PIX_SZ) then pix_sz = 5. ; kpc

  ;; QPQ7
  qpq_fil = '~/Dropbox/QSOPairs/qpq7_pairs.fits'
  qpq_strct = xmrdfits(qpq_fil, 1)

  if not keyword_set(FOREST) then begin ;; Cut out Lya forest
     good = where( (qpq_strct.z_fg+1)*wrest GT ((qpq_strct.z_bg+1)*1215.6701 + 20.), ngd)
     qpq_strct = qpq_strct[good]
  endif

  lowR = where(qpq_strct.R_phys LE 200., nlow)
  gd_strct = qpq_strct[lowR]

  ;; Push into this programs format
  rho_obs = gd_strct.R_phys ; kpc

  ;; Values
  EW_obs = gd_strct[lowR].EWLya

  gd = where( abs(gd_strct.metal_wrest-wrest) LT 1e-3 AND $
             gd_strct.metal_sigEW LT 0.3 AND $
              (gd_strct.flg_metal_eye MOD 256) GE 128 AND $
              gd_strct.flg_metal_EW EQ 1, ngd)
  sz = size(gd_strct.metal_wrest,/dimen)
  gdi = gd / sz[0]
  QPQ7_EW = (gd_strct.metal_EW)[gd] > 0.013
  flg_obs = replicate(1, ngd)

  ;; Upper Limits
  ulim = where( abs(gd_strct.metal_wrest-wrest) LT 1e-3 AND $
                gd_strct.metal_sigEW LT 0.3 AND $
                (gd_strct.flg_metal_eye MOD 256) GE 128 AND $
                gd_strct.flg_metal_EW EQ 3, nulim)
  if nulim GT 0 then begin
     ui = ulim / sz[0]
     lim_QPQ7_EW = (gd_strct.metal_EW)[ulim] 
     QPQ7_sigEW = (gd_strct.metal_sigEW)[ulim] 
     QPQ7_lim_EW = lim_QPQ7_EW > (2.* QPQ7_sigEW)
     QPQ7_EW = [QPQ7_EW, QPQ7_lim_EW]
     flg_obs = [flg_obs, replicate(2, nulim)]
  endif

  ;; set color map
  if ~keyword_set(colormap) then colormap=33
  if keyword_set(grey) then colormap=0
  if ~keyword_set(format) then format='(F5.2)'
  if keyword_set(white) then colormap=39
  if ~keyword_set(annotatecolor) then annotatecolor='black'

  if ~keyword_set(colpos) then colpos=[0.88,0.05,0.92,0.92]
  if ~keyword_set(reverse) then ctload, colormap
  if keyword_set(reverse) then ctload, colormap, /reverse


  ;; Generate the "Image"
  npix = round(grid_sz/pix_sz) * 2 + 1
  low_fits = fltarr(npix,npix)
  xval = pix_sz*findgen(npix) - grid_sz + 0.5 * pix_sz
  yval = pix_sz*findgen(npix) - grid_sz + 0.5 * pix_sz

  xgrid = xval # replicate(1., npix)
  ygrid = replicate(1., npix) # yval
  rgrid = sqrt(xgrid^2 + ygrid^2)

  if not keyword_set(RBINS) then $  ; kpc
     rbins = [ [0., 100], $
               [100, 200]] 
  nbins = (size(rbins,/dimen))[1]

  seed = -32188L
  for qq=0L,nbins-1 do begin
     ;; Grid
     gdr = where( rgrid GE rbins[0,qq] and rgrid LT rbins[1,qq], ngdr )
     if ngdr EQ 0 then continue

     ;; Observation
     gdo = where( rho_obs GE rbins[0,qq] and rho_obs LT rbins[1,qq], ngdo )

     if ngdo LT SAMP_LIM then stop ;; Too few data
     print, 'N samp = ', ngdo
     sub_colm = QPQ7_EW[gdo]
     sub_flg = flg_obs[gdo]

     ;; Random drawing
     rani = fix( randomu(seed, ngdr) * ngdo ) 
     ran_clm = sub_colm[ rani ]
     ran_flg = sub_flg[ rani ]

     ;; Lower limits
     low = where(ran_flg EQ 3, nlow)
     if nlow GT 0 then ran_clm[low] = ran_clm[low] + LOW_BOOST * randomu(seed, nlow)

     ;; Upper limits
     upper = where(ran_flg EQ 2, nupper)
     if nupper GT 0 then ran_clm[upper] = ran_clm[upper] - UPPER_BOOST * randomu(seed, nupper)

     ;; Add an error
     ran_clm = ran_clm + SIG_CLM * randomn(seed, ngdr)

     ;; Fill in
     low_fits[gdr] = ran_clm
  endfor
     
  ;; Remap to a 1kpc grid and add scatter
  new_bin = round(pix_sz)
  fits = rebin( low_fits, npix*new_bin, npix*new_bin)

  ;; Add scatter
  nfine = n_elements(fits)
  fits = fits + SIG_FINE * randomn(seed, nfine)

  ;; Smooth?
  
  ;;set zoom
  IF KEYWORD_SET(ZOOM) THEN BEGIN
      fsz=size(fits)
      xsize=fsz[1]-1
      ysize=fsz[2]-1
      DeltaX=(xsize-xsize/zoom)*0.5
      DeltaY=(ysize-ysize/zoom)*0.5
      fits=fits[DeltaX:(xsize-DeltaX),DeltaY:(ysize-DeltaY)]
      splog, 'Final size ', (size(fits))[1], (size(fits))[2]
      splog, 'Final size ', ps*(size(fits))[1], ps*(size(fits))[2]
      
  ENDIF
    
  IF keyword_set(PSFILE) THEN x_psopen, psfile,/maxs, /encapsulated,/color
    
  ;;set scales
  IF KEYWORD_SET(range) THEN scale=range ELSE scale=[MIN(fits),MAX(fits)]
  IF KEYWORD_SET(zscale) THEN scale=ZSCALE_RANGE(fits)
 
  
  ;;if white set, saturate the background white
  IF KEYWORD_SET(white) THEN BEGIN
      ;;;fix white into red
      ;TVLCT, [[255], [0], [0]], 255
      ;;;put black into white
      ;TVLCT, [[255], [255], [255]], 0
      low=where(fits lt scale[0],nl)
      if(nl gt 0) then fits[low]=abs(scale[1]*10)
   ENDIF
  

  if KEYWORD_SET(black) THEN TVLCT, [[0], [0], [0]], 0
      

  !x.style=1
  !y.style=1
  
  IF KEYWORD_SET(log) THEN BEGIN
      IMDISP, ALOG10(fits), RANGE=ALOG10(scale), /axis, XTickformat='(A1)', YTickformat='(A1)',$
        XTICKLEN=1D-5,YTICKLEN=1D-5, title=title, xtitle=xtitle, ytitle=ytitle, $
        margin=0.025, charsize=2.5, position=[0.04,0.04,0.94,0.94], margin=0.
      if ~keyword_set(nocolbar) then fsc_colorbar, range=ALOG10(scale), FORMAT=format,$
        /vertical, title=coltitle, POSITION=colpos, annotatecolor=annotatecolor, _extra=extra
  ENDIF ELSE BEGIN
      IMDISP, fits, RANGE=scale, /axis, XTickformat='(A1)', YTickformat='(A1)',$
        XTICKLEN=1D-5,YTICKLEN=1D-5, title=title, xtitle=xtitle , ytitle=ytitle,$
        charsize=2.5, position=[0.01,0.02,0.9,0.92], margin=0.
      if ~keyword_set(nocolbar) then fsc_colorbar, range=scale,FORMAT=format, /vertical,$
        title=coltitle, POSITION=colpos, charsize=2., annotatecolor=annotatecolor, _extra=extra
  ENDELSE
  
 if keyword_set(left) then start=0.2 else start=0.8

 if keyword_set(white) then textcol=fsc_color("black") else textcol=fsc_color("white")


  IF KEYWORD_SET(units) THEN BEGIN
     xmax=N_ELEMENTS(fits[0,*])-1
     ymax=N_ELEMENTS(fits[0,*])-1
     oplot, [start*xmax-1./PS,start*xmax], [0.1*ymax,0.1*ymax], THICK=10, color=textcol
     lab='$Form.LABEL' 
     lab=units
     xyouts, (0.9*start*xmax-1./PS), (0.05*ymax), lab, CHARSIZE=1.5,CHARTHICK=3, color=textcol
  ENDIF
  
  
  if keyword_Set(circ) then x_oplotcirc, circ/ps, x0=((size(fits))[1])*0.5, y0=((size(fits))[2])*0.5, $
                                         color=textcol

  IF KEYWORD_SET(text) THEN BEGIN
     ;loadct, 0
     xyouts, tlocation[0],tlocation[1], text, CHARSIZE=1.5,CHARTHICK=2.5 , color=textcol 
  ENDIF
  
  IF keyword_set(PSFILE) THEN x_psclose




  ;;RESTORE COLOR BAR
  ctload, 39


END
