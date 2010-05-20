;; 
; lineprof, 'line_delta.ps', /DELTA
; lineprof, 'line_lorentz.ps', /DELTA, /LOR
; lineprof, 'line_max.ps', /DELTA, /LOR, /MAX
; lineprof, 'line_voigt.ps', /DELTA, /LOR, /MAX, /VOIGT
pro lineprof, psfile, MAX=max, LOR=LOR, VOIGT=VOIGT, DELTA=delta


  lsz = 2.3

 ;; 
  c = x_constants()
  csz = 2.1
  nplt = 100000L
  if not keyword_set(PSFILE) then psfile = 'lineprof.ps'
  x_psopen, psfile, /maxs
  !p.multi=[0,1,1]
  clr = getcolor(/load)

  xrng = [-1000, 1000]
  yrng = [-15., 0.1]

  vplt = findgen(nplt)*(xrng[1])/(nplt-1)
  vplt = [-1*vplt,vplt]
  vplt = vplt[sort(vplt)]

  ;; Overall
  plot, [0],  [0],  color=clr.lightgray, $
    background=clr.black, charsize=csz,$
    xmargin=[7,1.5], ymargin=[3.5,0.5], xtitle='Relative Velocity (km/s)', $
    ytitle='log !9f!X(v)', /nodata, xthick=7, ythick=5, xstyle=1, ystyle=1, $
    yr=yrng, xr=xrng

  ;; Delta
  b = 1e-3
  tau = exp(-vplt^2/b^2)

  if keyword_set(DELTA) then oplot, vplt, alog10(tau), color=clr.green

  ;; Maxwell
  b = 30.
  if keyword_set(MAX) then begin
      tau = exp(-vplt^2/b^2) / (b*1e5 * sqrt(!dpi))
      oplot, vplt, alog10(tau), color=clr.cyan
  endif

  ;; Lorentzian
  gamma = 6.265d8
  nu0 = c.c / (1215.67 * 1e-8)
  dnu = vplt*1e5 / c.c * nu0 
  tau_nu = (gamma/4/!dpi^2) / ( dnu^2 + gamma^2/(4*!dpi)^2 )
  tau = tau_nu * nu0 / c.c 

  if keyword_set(LOR) then oplot, vplt, alog10(tau), color=clr.yellow

  ;; Voigt
  nuD = nu0 * b*1e5 / c.c
  a = gamma / (4*!dpi*nuD)
  u = dnu / nuD
  tau_nu = voigt(a,u) / nuD / sqrt(!dpi)
  tau = tau_nu * nu0 / c.c 

  if keyword_set(VOIGT) then oplot, vplt, alog10(tau), color=clr.lightgray

  ;; Label
  if keyword_set(DELTA) then $
    xyouts, xrng[0]+100., -1, 'Delta function', color=clr.green, charsiz=lsz
  if keyword_set(MAX) then $
    xyouts, xrng[0]+100., -3., 'Doppler (b=30 km/s)', color=clr.cyan, charsiz=lsz
  if keyword_set(LOR) then $
    xyouts, xrng[0]+100., -2., 'Natural (!9g!x= 6x10!u8!N)', color=clr.yellow, charsiz=lsz
  if keyword_set(VOIGT) then $
    xyouts, xrng[0]+100., -4., 'Voigt', color=clr.lightgray, charsiz=lsz

  x_psclose
  !p.multi=[0,1,1]

  return
end
