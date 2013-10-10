;; Plots HE0940 in a few ways
pro lpe_he0940, REST=rest, ZOOM=zoom

  if not keyword_set(CSZ) then csz = 2.1
  if not keyword_set(lSZ) then lsz = 1.8
  if not keyword_set(lSZ2) then lsz2 = 1.5

  ;; Data
  fx = x_readspec('~/Dropbox/llsz3mage/spectra/HE0940-1050_F.fits', $
                     wav=wv) 
  ;; Emission Plot
  psfil = 'lpe_he0940.ps'

  ;;
  x_psopen, psfil, /maxs
  !p.multi = [0,1,1]
  clr = getcolor(/load)
  pclr = clr.white
  ;pclr = clr.lightgray

  xmrg = [7,7]
  ymrg = [9,3]
  yrng=[0.0, 16]
  xrng=[4740., 5070]
  xplt = 4980.

  ytit = 'Relative Intensity'
  xtit = 'Observed Wavelength (Ang)'

  ;; Plot
  plot, [0], [0], color=pclr, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, ytitle=ytit, $
        xtitle=xtit, yrange=yrng, thick=5, ytickinter=5.0, $
        xrange=xrng, ystyle=1, xstyle=9, psym=1, /nodata ;, xtickint=1.
  
  oplot, wv, fx, color=pclr, psym=10, thick=2
  oplot, xrng, replicate(yrng[1]*0.999,2), color=pclr, thick=7

  lclr = clr.cyan
  xyouts, xplt, 3.0, 'HE0940-1050', color=lclr, charsi=lsz, align=0.
  xyouts, xplt+10, 1.7, 'z=3.076', color=lclr, charsi=lsz, align=0.

  ;; Distance on upper axis
  zem = 3.076
  wv_lbl = [4800., 4850., 4900,4950]
  zlya = wv_lbl/1215.6701 - 1.d
  nlbl = n_elements(zlya)

  for kk=0L,nlbl-1 do begin
     dist = abs(cosm_dist(zem,/w05map,/physical) - cosm_dist(zlya[kk],/w05map,/physical))
     xyouts, wv_lbl[kk], yrng[1]+0.3, strtrim(string(round(dist),format='(i2)'),2), $
             color=pclr, charsi=csz,  align=0.5
     oplot, replicate(wv_lbl[kk],2), yrng[1]+[0,-0.5], color=pclr, thick=5
  endfor
  xyouts, mean(xrng), yrng[1]+1.5, 'Distance from HE0940-1050 (pMpc)', $
          color=pclr, charsi=csz, align=0.5

;  xyouts, 8300, 13.0, 'Fumagalli et al. (2013)', color=lclr, charsi=lsz2, align=0.
        
  ;; End
  x_psclose
  !p.multi = [0,1,1]


  return

end
