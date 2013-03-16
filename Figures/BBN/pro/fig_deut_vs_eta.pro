;; This builds in Lya at z=2, 3, and 4
pro fig_deut_vs_eta

  if not keyword_set(CSZ) then csz = 2.1
  if not keyword_set(lSZ) then lsz = 1.8
  if not keyword_set(lSZ2) then lsz2 = 1.5

  ;; Emission Plot
  x_psopen, 'fig_deut_vs_eta.ps', /maxs
  !p.multi = [0,1,1]
  clr = getcolor(/load)

  xmrg = [7,17]
  ymrg = [10,6]
  yrng=[1e-6, 1e-2]
  xrng=[0.5, 30]

  ytit = 'D/H'
  xtit = '!9h!X!d10!N'

  ;; Plot
  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, ytitle=ytit, $
        xtitle=xtit, yrange=yrng, thick=5, ytickinter=1.0, $
        xrange=xrng, ystyle=1, xstyle=9, psym=1, /nodata, /xlog, /ylog
  
  ;; Eta formulation from Steidgman 2007
  eta = xrng[0] + (xrng[1]-xrng[0])*findgen(10000L)/9999.
  DH = 1d-5 * 2.68 * (6/eta)^(1.6)

  oplot, eta, DH, color=clr.blue

  ;; Omega_b axis
  omegabh2 = eta / 274.
  xrng2 = xrng / 274.
  axis, xaxis = 1, color = clr.black, charsi = csz, xrang = xrng2, $
        xsty = 1, xtit = '!9W!X!db!N h!u2!N', /save

;  xyouts, 8300, 1.30, 'Rafelski et al. (2012)', color=clr.darkgreen, charsi=lsz2, align=0.
        
  ;; End
  x_psclose
  !p.multi = [0,1,1]


  return

end
