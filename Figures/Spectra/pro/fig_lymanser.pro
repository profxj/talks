pro fig_lymanser

  if not keyword_set(PSFILE) then psfile = 'fig_lymanser.ps'

  ;; Read in spectrum
  fx = xmrdfits( getenv('XIDL_DIR')+'/SDSS/LLS/nhi16_19b30.fits',10)
  npnh = n_elements(fx) 
  wv_mod = 10^(2.7d + dindgen(npnh)*1e-4)

  ;; Plot
  if keyword_set( PSFILE ) then x_psopen, psfile, /maxs
  !p.multi=[0,1,1]
  clr = getcolor(/load)

  xrng = [900., 1230.]
  yrng = [-0.1, 1.1]
  csz = 2.

  plot, [0],  [0],  color=clr.lightgray, $
    background=clr.black, charsize=csz,$
    xmargin=[7,1.5], ymargin=[5,7], xtitle='Wavelength (Ang)', $
    ytitle='Normalized Flux', /nodata, xthick=7, ythick=7, xstyle=1, ystyle=1, $
    yr=yrng, xr=xrng

  oplot, wv_mod, fx, color=clr.lightgray, thick=7

  ;; Label
  lsz = 2.
  xyouts, 1200., 0.0, 'Ly!9a!X', color=clr.yellow, charsiz=lsz, alignm=0.5
  xyouts, 1025., 0.0, 'Ly!9b!X', color=clr.yellow, charsiz=lsz, alignm=0.5
  xyouts, 972., 0.0, 'Ly!9g!X', color=clr.yellow, charsiz=lsz, alignm=0.5
  xyouts, 949., 0.0, 'Ly!9d!X', color=clr.yellow, charsiz=lsz, alignm=0.5

  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]

  return
end
