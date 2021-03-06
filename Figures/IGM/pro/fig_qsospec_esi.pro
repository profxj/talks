;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro fig_qsospec_esi

  compile_opt strictarr

  ;; Get structure if necessary
  if not keyword_set( PSFILE ) then psfile = 'fig_qsospec_esi.ps'
  if not keyword_set( CSZ ) then csz = 2.3
  if not keyword_set( LSZ ) then lsz = 1.9
  if not keyword_set( BLSZ ) then blsz = 2.0

  if not keyword_set( MAXCHI ) then MAXCHI = 1.3
  if not keyword_set( WVMIN ) then wvmin = 3900. 
  if not keyword_set(OMEGA_M) then omega_m = 0.3
  if not keyword_set(H0) then H0 = 72.
  if not keyword_set(ZQSO) then zqso = 3.8
  if not keyword_set(FWHM) then fwhm = 2.

  xrng1=[800., 1500.]

  ;; Read in Data
  zem = 3.561
  scale = 1e15
  flux = x_readspec('~/Keck/ESI/RedData/Q2223+20/Q2223+20_xF.fits', wav=wave)
  flux = flux * scale
  conti = xmrdfits('~/Keck/ESI/RedData/Q2223+20/Q2223+20_c.fits')
  conti = conti*scale
  rwave = wave / (1+zem)
  gd = where(rwave GT xrng1[0] and rwave LT xrng1[1])

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; PLOT
  x_psopen, psfile, /maxs
  clr = getcolor(/load)
  !p.multi=[0,1,1]
  
  xmrg = [8,2]
  ymrg = [4.0,4]
  ;xrng2=xrng1*(1+zem)
  xrng2=[4200., 6600]

  thk = 5
;  yrng=[-0.02, 3.3]
  yrng=[-0.02, 2.]

  for kk=0,3 do begin

     ;; Start plot
     plot, [0], [0], color=clr.black, background=clr.white, $
           charsize=csz,$
           xmargin=xmrg, ymargin=ymrg, xtitle='Wavelength (Angstroms)', $
           ytitle='Brightness', yrange=yrng, thick=4, $
           xrange=xrng2, ystyle=1, xstyle=1, psym=1, /nodata ;, /ylog
     
     
     if kk LT 2 then $
        oplot, wave[gd], smooth(flux[gd],3), color=clr.black, thick=5, psym=10
     if kk GE 1 then $
        oplot, wave[gd], conti[gd], color=clr.red, thick=15

     ;; Silhouette
     silh = (conti-(smooth(flux,3)>0.)) > 0.
     if kk GE 3 then $
        oplot, wave[gd], silh, color=clr.darkgray, thick=5

     ;; Labels
     xlbl = 5700.
     if kk LT 2 then $
        xyouts, xlbl, 1.7, 'Keck/ESI Spectrum', color=clr.black, charsi=lsz, align=0.
     if kk GT 0 then $
        xyouts, xlbl, 1.55, 'Intrinsic QSO Spectrum', color=clr.red, charsi=lsz, align=0.
     if kk EQ 3 then $
        xyouts, xlbl, 1.4, 'Silhouette (Cosmic Web)', color=clr.darkgray, charsi=lsz, align=0.

  endfor

  ;; Close Ps
  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]

  print, 'fig_sdss:  All done!'
       
  return
end
      
      
