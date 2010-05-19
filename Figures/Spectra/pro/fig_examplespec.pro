pro fig_examplespec

  if not keyword_set(CSZ) then csz = 1.4
  if not keyword_set(lSZ) then lsz = 1.5

  ;; Star
  !p.multi = [0,1,1]

  if not keyword_set( starFIL ) then $
    starFIL = '/u/xavier/HST/STIS/HD195965/HD195965_A.fits'
  star = x_readspec(starfil, wav=wave, inf=2)

  x_psopen, 'fig_examplespec_star.ps', /maxs
  clr = getcolor(/load)
  yrng=[-0.05, 1.1]
  xrng=[1234., 1390]
  plot, [0], [0], color=clr.lightgray, background=clr.black, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, ytitle='Normalized Flux', $
        xtitle='Wavelength', yrange=yrng, thick=4, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata, $
        pos=[0.1, 0.75, 0.95, 0.95]

  oplot, wave, smooth(star,3), color=clr.lightgray, thick=3
  oplot, xrng, [0., 0.], color=clr.yellow, lines=2, thick=3

  x_psclose
  !p.multi = [0,1,1]

  ;; QSO spectrum
  if not keyword_set( qsoFIL ) then $
    qsoFIL = '/u/xavier/Keck/ESI/RedData/FJ0812+32/FJ0812+32_f.fits'
  qso = x_readspec(qsofil, wav=wave)

  x_psopen, 'fig_examplespec_qso.ps', /maxs
  clr = getcolor(/load)
  yrng=[-0.05, 1.1]
  xrng=[4100, 8950.]
  plot, [0], [0], color=clr.lightgray, background=clr.black, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, ytitle='Normalized Flux', $
        xtitle='Wavelength', yrange=yrng, thick=4, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata, $
        pos=[0.1, 0.75, 0.95, 0.95]

  oplot, wave, smooth(qso,3), color=clr.lightgray, thick=3
  oplot, xrng, [0., 0.], color=clr.yellow, lines=2, thick=3

  x_psclose
  !p.multi = [0,1,1]
  return

end
