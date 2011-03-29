;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro fig_solar_s2n

  compile_opt strictarr

  ;; Get structure if necessary
  if not keyword_set( PSFILE ) then psfile = 'fig_solar_s2n.ps'
  if not keyword_set( CSZ ) then csz = 1.8
  if not keyword_set( LSZ ) then lsz = 1.9
  if not keyword_set( BLSZ ) then blsz = 2.0

  if not keyword_set( MAXCHI ) then MAXCHI = 1.3
  if not keyword_set( WVMIN ) then wvmin = 3900. 
  if not keyword_set(OMEGA_M) then omega_m = 0.3
  if not keyword_set(H0) then H0 = 72.
  if not keyword_set(ZQSO) then zqso = 3.8
  if not keyword_set(FWHM) then fwhm = 2.

  xrng=[4350., 4372.]

  ;; Read in Data
  flux = xmrdfits('~/Stars/Solar/SolkurB_f.fits',0)
  wave = xmrdfits('~/Stars/Solar/SolkurB_l.fits',0)
  wv = wave[*,0]
  fx = flux[*,0]
  npix = n_elements(wv)

  ;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Degrade

  ;; Smooth
  dwv = wv[1] - wv[0]
  fwhm_pix = 12  ;; Corresponds to ~10 km/s FWHM
  nsmooth = fwhm_pix/(2.*sqrt(2*alog(2)))
  kernel = gauss_kernel(nsmooth)
  smth = convol(fx, kernel,/edge_wrap)
  
  ;; Rebin
  nlow = npix/4
  low_fx = congrid(smth, nlow)
  xval = congrid(wv, nlow)
  yval = x_addnoise(low_fx, 30.)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; PLOT
  x_psopen, psfile, /maxs
  clr = getcolor(/load)
  !p.multi=[0,1,1]
  
  xmrg = [8,2]
  ymrg = [15.0,4]

  thk = 5
;  yrng=[-0.02, 3.3]
  yrng=[0.0, 1.2]
  
  ;; Start plot
  plot, [0], [0], color=clr.black, background=clr.white, $
        charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, xtitle='Wavelength (Angstroms)', $
        ytitle='Normalized Intensity', yrange=yrng, thick=4, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata ;, /ylog
     
;  oplot, wave[*,0], flux[*,0], color=clr.black, thick=3, psym=10
  oplot, xval, yval, color=clr.black, thick=3, psym=10

  xyouts, 4351., 1.1, 'Solar Spectrum at R=30,000, S/N=30', color=clr.red, $
          charsi=lsz, align=0.

  ;; Close Ps
  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]

  print, 'fig_sdss:  All done!'
       
  return
end
      
      
