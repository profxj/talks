pro fig_lbg_vr, RREAL=rreal, STRCT=strct

  if not keyword_set( PSFILE ) then psfile = 'fig_lbg_vr.ps'
  if not keyword_set(CSZ) then csz = 1.3
  if not keyword_set(lSZ) then lsz = 1.5
  if not keyword_set(XNCOLORS) then xncolors=200L

  if not keyword_set(dv) then dv = 1.
  if not keyword_set(NPTS1) then npts1 = 10000L                  ; Log steps
  if not keyword_set(NPTS2) then npts2 = 10000L                   ; Log steps

  ;; LBG stuff
  if not keyword_set(r_min) then r_min = 1.0 ; kpc
  if not keyword_set(r_max) then r_max = 100. ; kpc
  if not keyword_set(gamma) then gamma = 0.5 ; Covering fraction
  if not keyword_set(fcmax) then fcmax = 0.6 ; Covering fraction
  if not keyword_set(alpha) then alpha = 1.3 ; Covering fraction
  if not keyword_set(vmax) then vmax = 800. ; Covering fraction
  if not keyword_set(WREST) then wrest = 2796.352d


  ;; Begin
  c = x_constants()
  rcut = 1.2
  rval_lo = r_min * 10.^(alog10(rcut) * findgen(npts1) / npts1) ; kpc
  rval_hi = (r_min*rcut) * 10.^(alog10(r_max/r_min/rcut) * findgen(npts2) / (npts2-1)) ; kpc
  rval = [rval_lo, rval_hi]
  npts = npts1+npts2
  dr = rval - shift(rval,1)  ; kpc
  dr[0] = dr[1] 

  ;; Wavelength array
  npix = 2000L
  wav = 10.^(alog10(2790.) + dindgen(npix)*1.449d-6)
  vel = (wav-wrest)/wrest * c.c/1e5
  dvel = median(vel-shift(vel,1))  ; Should be 1 km/s


  ;;;;;;;;;;;;;;;;
  ;; LBG
  A_lbg = vmax^2 * (1-alpha)
  print, 'A_LBG = ', A_LBG
  v_lbg = -1. * sqrt(A_lbg/(1-alpha)) * sqrt(r_min^(1-alpha) - rval^(1-alpha))
  fc_lbg = fcmax * (rval/r_min)^(-1*gamma)
  I_lbg = 1 - fc_lbg
;  I_v = 1. - fcmax * (1 - (1-alpha)*v_lbg^2/A_lbg)^(-1*gamma/(1-alpha))

  ;;;;;;;;;;;;;;;;
  ;; Model (i)  [Sobolev]

  ;; Key stuff
  tau_r = -1*alog(I_lbg)
  dvdr = sqrt(A_lbg/(1-alpha)) * 0.5 / sqrt(1-rval^(1-alpha)) * (alpha-1) * rval^(-1*alpha)
  dvdr[0] = 2*dvdr[1] ; Kludge

  ;; Density
  getfnam, wrest, fval, nam
  kappa_l = !pi*c.e^2/c.me/c.c * fval
  n_Mg = tau_r * (dvdr*1e5/c.kpc) / (wrest*1e-8) / kappa_l 
  METAL = -0.3
  DUST = 0.1
  b_val = 5.
  nr_i =  n_Mg /  (10.^(7.53-12.+METAL) * DUST)  ;; Hydrogen
;  x_splot, v_lbg, n_H, /block
;  stop
  ;; dlogn/dlogr
  lgn = alog(nr_i)
  dlgn = lgn - shift(lgn,1)
  dlgn[0] = dlgn[1]
  lgr = alog(rval-1)
  dlgr = lgr - shift(lgr,1)
  dlgr[0] = dlgr[1]
  dlgnr = dlgn/dlgr
  ;x_splot, rval-1, dlgn/dlgr, /bloc
  lowr = where(rval LT 1.1)
  low_dlnr = median(dlgnr[lowr])
  print, 'dln[n] / dln[r]:  Low -- ',  low_dlnr
  highr = where(rval GT 4)
  hi_dlnr = median(dlgnr[highr])
  print, 'dln[n] / dln[r]:  High -- ',  hi_dlnr

  ;; Optical depth
  mgii = x_setline(wrest)
  lines = replicate(mgii, npts)
  lines.b = b_val
  lines.N = alog10(dr * c.kpc * nr_i * 10.^(7.53-12.+METAL) * DUST)
  lines.zabs = v_lbg/3e5
;  printcol, rval[0:10], dr, v_lbg, lines.N 
;  stop

;  fx_i = x_voigt(wav, lines, /nosmooth, TAU=tau)
;  stop

;  x_splot, rval, v_lbg, /blo
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Plot
  x_psopen, psfile, /maxs
  !p.multi = [0,1,1]
  clr = getcolor(/load)
  lclr = clr.white

  xmrg = [8,7]
  ymrg = [3.7,0.5]

  ;; Radial stuff
  xrng=[1, 100]
;  yrng=[1e-6, 1.]
  yrng=[0., 800]
  csz = 2.2
  plot, [0], [0], color=lclr, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, $
        ytitle='v!dr!N (km/s)', $
        xtitle='r (kpc)', yrange=yrng, thick=4, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata, /xlog, /noeras;, $

  oplot, rval, abs(v_lbg), color=lclr

  lsz = 2.2
;  xyouts, 10., 50., '!9t!X!dS!N(r) ~ n(r) * (dv/dr)!u-1!N', color=lclr, charsi=lsz
  xx=30.
  xyouts, xx, 500., 'Very little scattering', color=clr.cyan, charsi=lsz, align=0.5
  xyouts, xx, 450., 'at large radii', color=clr.cyan, charsi=lsz, align=0.5


  if keyword_set( PSFILE ) then x_psclose
  !p.multi = [0,1,1]

  return
end
