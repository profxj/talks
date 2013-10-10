;;  Plots GRB 050730 for comparison
pro fig_grb050730, REST=rest, ZOOM=zoom, GRB=grb

  if not keyword_set(CSZ) then csz = 2.1
  if not keyword_set(lSZ) then lsz = 1.8
  if not keyword_set(lSZ2) then lsz2 = 1.5

  ;; Data
  fxb = x_readspec('~/GRB/data/050730/MIKE/FSpec/GRB050730a_b_F.fits', wav=wvb) 
  fxr = x_readspec('~/GRB/data/050730/MIKE/FSpec/GRB050730a_r_F.fits', wav=wvr) 
  gdb = where(wvb LT 5010)
  gdr = where(wvr GT 5010)
  fx = [fxb[gdb], fxr[gdr]]
  wv = [wvb[gdb], wvr[gdr]]
  fxp = where(abs(wv-6301.3) LT 0.5)
  fx[fxp] = 0.57
  fxp = where(abs(wv-6171.8) LT 0.5)
  fx[fxp] = 0.28
  fxp = where(abs(wv-6364.9) LT 0.5)
  fx[fxp] = 0.67

  ;; Emission Plot
  psfil = 'fig_grb050730.ps'

  ;;
  x_psopen, psfil, /maxs
  !p.multi = [0,1,1]
  clr = getcolor(/load)
  pclr = clr.white
  ;pclr = clr.lightgray

  xmrg = [7,7]
  ymrg = [11,1]
  yrng=[-0.02, 1.4]
  xrng=[900., 1350]

  ytit = 'Relative Intensity'
  zgrb = 3.96855
  xtit = 'Rest Wavelength (Ang)'
  wv = wv / (1+zgrb)

  ;; Plot
  plot, [0], [0], color=pclr, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, ytitle=ytit, $
        xtitle=xtit, yrange=yrng, thick=5, $;ytickinter=5.0, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata ;, xtickint=1.
  
  oplot, xrng, [0., 0], color=clr.green, linest=1, thick=3
  oplot, wv, smooth(fx,3), color=pclr, psym=10, thick=2

  lclr = clr.cyan
  xplt = 1100.
  xyouts, xplt, 1.3, 'GRB 057030  (Magellan/MIKE)', color=lclr, charsi=lsz, align=0.
  xyouts, xplt+100, 1.20, 'z=3.969', color=lclr, charsi=lsz, align=0.

;  xyouts, 8300, 13.0, 'Fumagalli et al. (2013)', color=lclr, charsi=lsz2, align=0.
        
  ;; End
  x_psclose
  !p.multi = [0,1,1]


  return

end
