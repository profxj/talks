;; Shows Best estimates of D/H and thereby Omega b
pro fig_deut_omegab

  if not keyword_set(CSZ) then csz = 1.7
  if not keyword_set(lSZ) then lsz = 1.8
  if not keyword_set(lSZ2) then lsz2 = 1.5

  ;; Emission Plot
  x_psopen, 'fig_deut_omegab.ps', /maxs
  !p.multi = [0,1,1]
  clr = getcolor(/load)

  ;; Eta formulation from Steidgman 2007
  eta = 2 + 18*findgen(10000L)/9999.
  DH = 1d-5 * 2.68 * (6/eta)^(1.6)
  omegabh2 = eta / 274.

  xmrg = [7,27]
  ymrg = [10,6]
  yrng=[1e-5, 1e-4]
  xrng=[0.01, 0.04]

  ytit = 'D/H'
  xtit = '!9W!X!db!N h!u2!N'

  ;; Plot
  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, ytitle=ytit, $
        xtitle=xtit, yrange=yrng, thick=5, xtickinter=0.01, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata, /ylog
  
  ;; Value
  DHp = -4.585 ; 0.02
  x_curvefill, xrng, replicate(10.^(DHp+0.02),2), replicate(10.^(DHp-0.02),2), $
               color=clr.pink

  oplot, omegabh2, DH, color=clr.blue

  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, ytitle=ytit, $
        xtitle=xtit, yrange=yrng, thick=5, xtickinter=0.01, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata, /ylog, /noerase

;  xyouts, 8300, 1.30, 'Rafelski et al. (2012)', color=clr.darkgreen, charsi=lsz2, align=0.
        
  x_psclose
  !p.multi = [0,1,1]

  ;stop

  return

end
