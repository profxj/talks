;; 
pro fig_lowhigh_spec, psfile

  lsz = 1.5
  csz = 2.5

 ;; 
  if not keyword_set(PSFILE) then psfile = 'fig_lowhigh_spec.ps'


  x_psopen, psfile, /portrait
  !p.multi=[0,1,3]
  clr = getcolor(/load)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; GRB 050401
  grb050401 = '/u/xavier/GRB/data/050401/VLT/FORS2/norm.fits'
  grb = x_readspec(grb050401, wav=wave)
  xrng = [4000., 8990]
  yrng = [-0.1, 1.3]

  ;; Plot
  plot, [0],  [0],  color=clr.lightgray, $
    background=clr.black, charsize=csz,$
    xmargin=[7,2.5], ymargin=[3.5,0.5], xtitle='Wavelength (Ang)', $
    ytitle='Normalized Flux', /nodata, xthick=5, ythick=5, xstyle=1, ystyle=1, $
    yr=yrng, xr=xrng

  oplot, wave, grb, color=clr.lightgray, thick=5

  ;; Label
  xyouts, xrng[1]-1200., 0.15, 'GRB 050401', color=clr.yellow, charsiz=lsz, align=0.
  xyouts, xrng[1]-1200., 0.0, 'Watson+ 06', color=clr.yellow, charsiz=lsz, align=0.

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; GRB 050401
  grb051111 = '/u/xavier/Keck/HIRES/RedData/GRB051111/GRB051111_f.fits'
  grb = x_readspec(grb051111, wav=wave)

  xrng = [6010., 6090]
  yrng = [-0.1, 1.2]

  ;; Plot
  plot, [0],  [0],  color=clr.lightgray, $
    background=clr.black, charsize=csz,$
    xmargin=[7,2.5], ymargin=[3.5,0.5], xtitle='Wavelength (Ang)', $
    ytitle='Normalized Flux', /nodata, xthick=5, ythick=5, xstyle=1, ystyle=1, $
    yr=yrng, xr=xrng

  oplot, wave, grb, color=clr.lightgray, thick=5

  ;; Label
  xyouts, xrng[0]+5., 0.15, 'GRB 051111', color=clr.yellow, charsiz=lsz, align=0.
  xyouts, xrng[0]+5., 0.0, 'Pro+ 06', color=clr.yellow, charsiz=lsz, align=0.

  x_psclose
  !p.multi=[0,1,1]
;  x_splot, wave, smooth(grb,3), /bloc
;  stop
  

  return
end
