;; Shows Best estimates of D/H and thereby Omega b
pro fig_q1946_lyaf

  if not keyword_set(CSZ) then csz = 2.3
  if not keyword_set(lSZ) then lsz = 2.0
  if not keyword_set(lSZ2) then lsz2 = 1.5

  ;; Read
  q1946 = x_readspec('~/Keck/HIRES/RedData/Q1946+76/Q1946+76UV_f.fits', /struct)


  ;; Emission Plot
  x_psopen, 'fig_q1946_lyaf.ps', /maxs
  !p.multi = [0,1,3]
  clr = getcolor(/load)

  ;; Eta formulation from Steidgman 2007
  xmrg = [9,1]
  ymrg = [3,0.1]


  xrng=[3650, 5025]
  yrng=[0., 1.05]

  ytit = 'Normalized Flux'
  xtit = 'Wavelength (Ang)'

  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Q1946
  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, ytitle=ytit, $
        xtitle=xtit, yrange=yrng, thick=5, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata
     
  ;; Data
  oplot, q1946.wv, q1946.fx, psym=10, color=clr.black, thick=3

     
  x_psclose
  !p.multi = [0,1,1]

  return

end
