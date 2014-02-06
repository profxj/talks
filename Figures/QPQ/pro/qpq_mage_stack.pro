;; Plots HE0940 in a few ways
pro qpq_mage_stack, REST=rest, ZOOM=zoom

  if not keyword_set(CSZ) then csz = 2.1
  if not keyword_set(lSZ) then lsz = 1.8
  if not keyword_set(lSZ2) then lsz2 = 1.5

  ;; Data
  stack = xmrdfits('../MFP/llsz3mage_stack.fits',1)
  gd = where(finite(stack.flux_mean))
  fx = stack.flux_mean[gd]
  wv = stack.wave[gd]

  ;; Emission Plot
  psfil = 'qpq_mage_stack.ps'

  ;;
  x_psopen, psfil, /maxs
  !p.multi = [0,1,1]
  clr = getcolor(/load)
  ;pclr = clr.white
  fclr = clr.lightgray
  pclr = clr.white

  xmrg = [7,7]
  ymrg = [11,1]
  yrng=[0.0, 5]
  ;xrng = [3390, 5150.]/(1+3.076) ;; To match HE0940
  xrng = [1155., 1250.] ;; Ang
  xplt = 4000.
  REST=1

  ytit = 'Relative Intensity'
  zem = 3.088
  xtit = 'Rest Wavelength (Ang)'
  xplt = xplt / (1+zem)

  ;; Plot
  plot, [0], [0], color=fclr, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, ytitle=ytit, $
        xtitle=xtit, yrange=yrng, thick=5, ytickinter=1.0, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata ;, xtickint=1.
  
  oplot, wv, fx, color=pclr, psym=10, thick=2

  lclr = clr.cyan
  xyouts, 1160., 4.3, 'Stacked Quasars (Magellan/MagE)', color=lclr, charsi=lsz, align=0.
  ;xyouts, xplt+20, 1.65, 'z=3.088', color=lclr, charsi=lsz, align=0.

;  xyouts, 8300, 13.0, 'Fumagalli et al. (2013)', color=lclr, charsi=lsz2, align=0.
        
  ;; End
  x_psclose
  !p.multi = [0,1,1]


  return

end
