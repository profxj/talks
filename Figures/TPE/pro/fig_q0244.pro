;; Plots HE0940 in a few ways
pro fig_q0244, REST=rest, ZOOM=zoom

  if not keyword_set(CSZ) then csz = 2.1
  if not keyword_set(lSZ) then lsz = 1.8
  if not keyword_set(lSZ2) then lsz2 = 1.5

  ;; Data
  fx = x_readspec('~/Dropbox/llsz3mage/spectra/Q0244-302_F.fits', $
                     wav=wv) 
  ;; Emission Plot
  psfil = 'fig_q0244.ps'

  ;;
  x_psopen, psfil, /maxs
  !p.multi = [0,1,1]
  clr = getcolor(/load)
  pclr = clr.white
  ;pclr = clr.lightgray

  xmrg = [7,7]
  ymrg = [11,1]
  yrng=[0.0, 2]
  xrng = [3390, 5150.]/(1+3.076) ;; To match HE0940
  xplt = 4000.
  REST=1

  ytit = 'Relative Intensity'
  zem = 3.088
  xtit = 'Rest Wavelength (Ang)'
  wv = wv / (1+zem)
  xplt = xplt / (1+zem)

  ;; Plot
  plot, [0], [0], color=pclr, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, ytitle=ytit, $
        xtitle=xtit, yrange=yrng, thick=5, ytickinter=5.0, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata ;, xtickint=1.
  
  oplot, wv, fx, color=pclr, psym=10, thick=2

  lclr = clr.cyan
  xyouts, xplt, 1.8, 'Q0244-302  (Magellan/MagE)', color=lclr, charsi=lsz, align=0.
  xyouts, xplt+20, 1.65, 'z=3.088', color=lclr, charsi=lsz, align=0.

;  xyouts, 8300, 13.0, 'Fumagalli et al. (2013)', color=lclr, charsi=lsz2, align=0.
        
  ;; End
  x_psclose
  !p.multi = [0,1,1]


  return

end
