;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro fig_q1937

  compile_opt strictarr

  ;; Get structure if necessary
  if not keyword_set( PSFILE ) then psfile = 'fig_q1937.ps'
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
  fx_hires = xmrdfits('~/Keck/HIRES/RedData/Q1937-1009/q1937-hires.fits', 0, hhead) 
  npix = n_elements(fx_hires)
  wv_hires = 
  fx_lris = x_readspec('~/Keck/HIRES/RedData/Q1937-1009/q1937.fits', inf=2, wav=wv_lris)
  gd = where(rwave GT xrng1[0] and rwave LT xrng1[1])

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; PLOT
  x_psopen, psfile, /maxs
  clr = getcolor(/load)
  !p.multi=[0,1,1]
  
  
  ;; LRIS first
  xmrg = [8,2]
  ymrg = [15.0,1]
  xrng=[4000., 6500]

  thk = 5
  yrng=[0., 25.]

  ;; Start plot
  plot, [0], [0], color=clr.black, background=clr.white, $
        charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, xtitle='Wavelength (Angstroms)', $
        ytitle='Intensity', yrange=yrng, thick=4, $
        xrange=xrng2, ystyle=1, xstyle=1, psym=1, /nodata ;, /ylog
     
     
  oplot, wv_lris, fx_lris, color=clr.black, thick=5, psym=10
  ;; Close Ps
  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]

  print, 'fig_sdss:  All done!'
       
  return
end
      
      
