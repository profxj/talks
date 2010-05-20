;; 
pro fig_tauvslambda, psfile, MAX=max, LOR=LOR, VOIGT=VOIGT, DELTA=delta

  lsz = 2.0
  csz = 1.5

 ;; 
  c = x_constants()
  nplt = 100000L
  if not keyword_set(PSFILE) then psfile = 'fig_tauvslambda.ps'
  x_psopen, psfile, /portrait
  !p.multi=[0,1,2]
  clr = getcolor(/load)

  xrng = [-1000, 1000]
  yrng = [-8., 1.]

  ;; Velocity
  vplt = findgen(nplt)*(xrng[1])/(nplt-1)
  vplt = [-1*vplt,vplt]
  vplt = vplt[sort(vplt)]

  ;; Voigt
  b = 30.
  gamma = 6.265d8
  nu0 = c.c / (1215.67 * 1e-8)
  dnu = vplt*1e5 / c.c * nu0 
  nuD = nu0 * b*1e5 / c.c
  a = gamma / (4*!dpi*nuD)
  u = dnu / nuD
  tau_nu = voigt(a,u) / nuD / sqrt(!dpi)
  tau = tau_nu / max(tau_nu)

  ;; Optical depth
  plot, [0],  [0],  color=clr.lightgray, $
    background=clr.black, charsize=csz,$
    xmargin=[7,2.5], ymargin=[3.5,0.5], xtitle='Relative Velocity (km/s)', $
    ytitle='log !9t!X(v)', /nodata, xthick=5, ythick=5, xstyle=1, ystyle=1, $
    yr=yrng, xr=xrng

  oplot, vplt, alog10(tau), color=clr.lightgray

  ;; Label
  xyouts, xrng[0]+100., 0., '!9t!X!d0!N=1', color=clr.lightgray, charsiz=lsz

  ;; Observable
  yspaces = replicate(' ',30) 
  yrng = [0., 1.05]
  plot, [0],  [0],  color=clr.lightgray, $
    background=clr.black, charsize=csz,$
    xmargin=[7,2.5], ymargin=[3.5,0.5], xtitle='Relative Velocity (km/s)', $
    ytitle='Intensity', /nodata, xthick=5, ythick=5, xstyle=1, ystyle=1, $
    yr=yrng, xr=xrng, ytickn=yspaces

  oplot, xrng, [1.,1.], color=clr.orange, linesty=1
  oplot, vplt, exp(-1*tau), color=clr.lightgray

  ;; Label
  xyouts, xrng[0]-150., 0.95, 'I!d!9n!X!u0!N', color=clr.lightgray, charsiz=lsz


  x_psclose
  !p.multi=[0,1,1]

  return
end
