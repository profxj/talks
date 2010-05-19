pro fig_silhouette

  if not keyword_set( SPECFIL ) then $
    SPECFIL = '/u/xavier/HST/STIS/HD195965/HD195965_A.fits'
  if not keyword_set(CSZ) then csz = 2.1
  if not keyword_set(lSZ) then lsz = 1.5

  spec = x_readspec(specfil, wav=wave, inf=2)
  smth_spec = smooth(spec,5)
  emiss = 2. - smth_spec


  ;; Emission Plot
  xspaces = replicate(' ',30) 
  x_psopen, 'fig_silhoutte_emiss.ps', /maxs
  !p.multi = [0,1,1]
  clr = getcolor(/load)
  xmrg = [8,1]
  ymrg = [4.0,1]
  yrng=[0., 2.3]
  xrng=[1300., 1320]
  plot, [0], [0], color=clr.lightgray, background=clr.black, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, ytitle='Relative Flux', $
        xtitle='Wavelength', yrange=yrng, thick=4, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata, $
        xtickn=xspaces

  oplot, wave, emiss, color=clr.lightgray, thick=5

  x_psclose
  !p.multi = [0,1,1]

  ;; Absorption Plot
  xspaces = replicate(' ',30) 
  x_psopen, 'fig_silhoutte_abs.ps', /maxs
  !p.multi = [0,1,1]
  clr = getcolor(/load)
  xmrg = [8,1]
  ymrg = [4.0,1]
  yrng=[-0.05, 1.1]
  xrng=[1300., 1320]
  plot, [0], [0], color=clr.lightgray, background=clr.black, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, ytitle='Relative Flux', $
        xtitle='Wavelength', yrange=yrng, thick=4, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata, $
        xtickn=xspaces

  oplot, wave, smth_spec, color=clr.lightgray, thick=5
  oplot, xrng, [0., 0.], color=clr.yellow, lines=2, thick=3

  x_psclose

  return

end
