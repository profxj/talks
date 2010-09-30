pro fig_lbg_spec

  if not keyword_set( PSFILE ) then psfile = 'fig_lbg_spec.ps'
  if not keyword_set(PAD_FRAC) then pad_frac = 0.1
  if not keyword_set(CSZ) then csz = 1.3
  if not keyword_set(CSZ2) then csz2 = 2.3
  if not keyword_set(lSZ) then lsz = 1.5
  if not keyword_set(lSZ2) then lsz2 = 1.3
  if not keyword_set(XNCOLORS) then xncolors=200L
  if not keyword_set(YSTP) then ystp = 0.07

  xlbl = 0.05
  xlbl2 = 0.05
  ylbl = 0.90
  thk = 9

  close, /all

  ;;; BEGIN PLOTS
  x_psopen, psfile, /portrait
  !p.multi = [0,1,2]
  clr = getcolor(/load)
  clrs = x_setclrs()

  xmrg = [8,4]
  ymrg = [3.0,3.5]

  ;; Plot FeII
  yrng=[0., 1.8]
  xrng=[2580, 2618]
  xcut = 2605.
  off = 1.

  for ss=0,1 do begin
     ;; Plot FeII
     case ss of 
        0: begin
           yrng=[0., 2.0] 
           ysty = 9
           wvmnx = [xrng[0], xcut-off]
        end
        1: begin
              yrng=[0.90,1.9]
              ysty = 5
              wvmnx = [xcut+off,xrng[1]]
           end
        else: stop
     endcase

     plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
           xmargin=xmrg, ymargin=ymrg, ytitle='Relative Flux', $
           xtitle='Wavelength (Ang)', yrange=yrng, thick=4, $
           xrange=xrng, ystyle=ysty, xstyle=9, psym=1, /nodata, /noerase, $
           xthi=thk, ythi=thk
     if ss EQ 1 then axis, yaxis=1, charsiz=csz, ysty=1, xrang=yrng, ytickint=0.2, xthi=thk, $
                           ythi=thk
  
     ;; Sobolev
     Fe_fil = getenv('WFTY')+'/Analysis/LBG/Sobolev/Output/spec_FeII_lbg_sobolev.dat'
     readcol, Fe_fil, wv, fx, noscatt_fx, /silen
     nrm = median(fx[where(wv GT 2605)])
     pix = where(wv GT wvmnx[0] and wv LT wvmnx[1])
     oplot, wv[pix], fx[pix]/nrm, color=clr.red, psym=10, thick=thk

     ;; f_c method
     Fe_fil = getenv('WFTY')+'/Analysis/LBG/Covering/Output/spec_FeII_lbg_covering.dat'
     readcol, Fe_fil, wv, fx, noscatt_fx, /silen
     nrm = median(fx[where(wv GT 2605)])
     pix = where(wv GT wvmnx[0] and wv LT wvmnx[1])
     oplot, wv[pix], fx[pix]/nrm, color=clr.black, psym=10, thick=thk
     nrm = median(noscatt_fx[where(wv GT 2605)])
     if ss EQ 0 then $
        oplot, wv[pix], noscatt_fx[pix]/nrm, color=clr.darkgray, thick=thk 
  endfor
  xrng2 = (xrng/2600.173 - 1)*3e5
  axis, xaxis=1, charsiz=csz, xsty=1, xrang=xrng2, xtitl='Velocity Relative to FeII 2600 (km/s)'

  oplot, replicate(2586.650,2), yrng, color=clr.orange, linesty=2, thick=5
  oplot, replicate(2600.173,2), yrng, color=clr.orange, linesty=2, thick=5
  oplot, replicate(2612.6542,2), yrng, color=clr.orange, linesty=2, thick=5
  xyouts, xrng[0]+xlbl*(xrng[1]-xrng[0]), yrng[0]+(yrng[1]-yrng[0])*ylbl, $
          'FeII', color=clr.black, charsiz=lsz
;     oplot, xrng, [1., 1.], color=clr.red, linestyle=1, thick=1
  x_curvefill, [xcut-off,xcut+off], [0., 0.], [10., 10], color=clr.tan
  
  ;;;;;;;;;;;;;;;;;;;;;;;;
  ;;; MgII

  !p.multi = [1,1,2]
  ;; Plot MgII
  yrng=[0., 2.5]
  xrng=[2786., 2809.8]
  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, ytitle='Relative Flux', $
        xtitle='Wavelength (Ang)', yrange=yrng, thick=4, $
        xrange=xrng, ystyle=1, xstyle=9, psym=1, /nodata, $
        xthi=thk, ythi=thk

  ;; Sobolev
  Mg_fil = getenv('WFTY')+'/Analysis/LBG/Sobolev/Output/spec_MgII_lbg_sobolev.dat'
  readcol, Mg_fil, wv, fx, /silen
  nrm = median(fx[where(wv GT 2815)])
  oplot, wv, fx/nrm, color=clr.red, psym=10, thick=thk

  ;; Covering
  Mg_fil = getenv('WFTY')+'/Analysis/LBG/Covering/Output/spec_MgII_lbg_covering.dat'
  readcol, Mg_fil, wv, fx, noscatt_fx, /silen
  nrm = median(fx[where(wv GT 2815)])
  oplot, wv, fx/nrm, color=clr.black, psym=10, thick=thk

  nrm = median(noscatt_fx[where(wv GT 2815)])
  oplot, wv, noscatt_fx/nrm, color=clr.darkgray, thick=thk

  xrng2 = (xrng/2796.352 - 1)*3e5
  axis, xaxis=1, charsiz=csz, xsty=1, xrang=xrng2, xtitl='Velocity Relative to MgII 2796 (km/s)', $
        xthi=thk, ythi=thk

  oplot, replicate(2796.352,2), yrng, color=clr.orange, linesty=2, thick=5
  oplot, replicate(2803.531,2), yrng, color=clr.orange, linesty=2, thick=5
  xyouts, xrng[0]+xlbl*(xrng[1]-xrng[0]), yrng[1]*ylbl, $
          'MgII', color=clr.black, charsiz=lsz
;  oplot, xrng, [1., 1.], color=clr.red, linestyle=1, thick=1

  if keyword_set( PSFILE ) then x_psclose
  !p.multi = [0,1,1]
  close, /all

  return

end
