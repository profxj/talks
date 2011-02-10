;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro fig_q1937

  compile_opt strictarr

  ;; Get structure if necessary
  if not keyword_set( PSFILE ) then psfile = 'fig_q1937.ps'
  if not keyword_set( CSZ ) then csz = 1.9
  if not keyword_set( LSZ ) then lsz = 1.9
  if not keyword_set( BLSZ ) then blsz = 2.0

  if not keyword_set( MAXCHI ) then MAXCHI = 1.3
  if not keyword_set( WVMIN ) then wvmin = 3900. 
  if not keyword_set(OMEGA_M) then omega_m = 0.3
  if not keyword_set(H0) then H0 = 72.
  if not keyword_set(ZQSO) then zqso = 3.8
  if not keyword_set(FWHM) then fwhm = 2.


  ;; Read in Data
  zem = 3.561
  fx_hires = xmrdfits('~/Keck/HIRES/RedData/Q1937-1009/q1937-hires.fits', 0, hhead) 
  npix = n_elements(fx_hires)
  wv_hires = 10.d^(sxpar(hhead,'CRVAL1') + dindgen(sxpar(hhead,'NAXIS1'))*sxpar(hhead,'CDELT1'))
  fx_lris = x_readspec('~/Keck/HIRES/RedData/Q1937-1009/q1937.fits', inf=2, wav=wv_lris)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; PLOT
  x_psopen, psfile, /maxs
  clr = getcolor(/load)
  !p.multi=[0,1,1]
  
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; LRIS first
  xmrg = [8,2]
  ymrg = [19.0,1]
  xrng=[4000., 6500]

  thk = 5
  yrng=[0., 25.]

  ;; Start plot
  plot, [0], [0], color=clr.black, background=clr.white, $
        charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, xtitle='Wavelength (Angstroms)', $
        ytitle='Intensity', yrange=yrng, thick=4, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata ;, /ylog
     
     
  oplot, wv_lris, fx_lris, color=clr.black, thick=4, psym=10
  xyouts, 4150., 20., 'Keck/LRIS Spectrum of Q1937-1009', color=clr.darkgreen, charsi=lsz, align=0.
  xyouts, 4150., 17., 'Tytler et al. (1996)', color=clr.darkgreen, charsi=lsz, align=0.

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; HIRES next
  xmrg = [8,2]
  ymrg = [4.0,1]
  xrng=[5550., 5567]

  thk = 5
  yrng=[-0.05, 0.7]

  ;; Start plot
  plot, [0], [0], color=clr.black, background=clr.white, $
        charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, xtitle='Wavelength (Angstroms)', $
        ytitle='Intensity', yrange=yrng, thick=6, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata ;, /ylog
     
     
  oplot, wv_hires, fx_hires, color=clr.black, thick=7, psym=10

  xyouts, 5551., 0.65, 'Keck/HIRES Spectrum of Q1937-1009', color=clr.red, charsi=lsz, align=0.
  xyouts, 5551., 0.62, 'Tytler et al. (1996)', color=clr.red, charsi=lsz, align=0.


  ;; Close Ps
  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]

  print, 'fig_sdss:  All done!'
       
  return
end
      
      
