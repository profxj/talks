pro fig_mgii, RREAL=rreal

  if not keyword_set( PSFILE ) then psfile = 'fig_mgii.ps'
  if not keyword_set( GRID_FILE ) then $
    Mg_fil = getenv('WFTY')+'/Analysis/Fiducial/Output/spec_MgII_fiducial.dat'
  if not keyword_set(PAD_FRAC) then pad_frac = 0.1
  if not keyword_set(CSZ) then csz = 1.9
  if not keyword_set(lSZ) then lsz = 1.9

  if not keyword_set(MGII_REST) then mgii_rest = 2796.35d

  ;; Read and normalize
  readcol, Mg_fil, wv, fx, noscatt_fx, /silen
  nrm = median(fx[where(wv GT 2815)])
  fx = fx/nrm
  nrm2 = median(noscatt_fx[where(wv GT 2815)])
  noscatt_fx = noscatt_fx/nrm2

  thk = 11
  ;;; BEGIN PLOTS
  x_psopen, psfile, /maxs
  clr = getcolor(/load)

  for qq=0,1 do begin
     ;; Plot MgII
     yrng=[-0.05, 2.8]
     xrng=[2786., 2812]
     plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
           ytitle='Normalized Flux', ymarg=[5,5], xmarg=[7,1], $
           xtitle='Wavelength (Ang)', yrange=yrng, thick=thk, xthic=thk, ythic=thk, $
           xrange=xrng, ystyle=1, xstyle=9, psym=1, /nodata
     
     oplot, wv, fx, color=clr.black, psym=10, thick=thk
     
     oplot, replicate(2796.352,2), yrng, color=clr.blue, linesty=2
     oplot, replicate(2803.531,2), yrng, color=clr.blue, linesty=2
     xlbl = 0.05
     ylbl = 0.90
     xyouts, xrng[0]+xlbl*(xrng[1]-xrng[0]), yrng[1]*ylbl, $
             'MgII', color=clr.black, charsiz=lsz
     oplot, xrng, [1., 1.], color=clr.green, linestyle=1, thick=5

     if qq GT 0 then oplot, wv, noscatt_fx, color=clr.darkgray, thick=thk 
     
     xrng2 = (xrng/2796.352 - 1)*3e5
     axis, xaxis=1, charsiz=csz, xsty=1, xrang=xrng2, xtitl='Velocity (km/s) Relative to MgII 2796', xthic=thk

  endfor

  if keyword_set( PSFILE ) then x_psclose
  !p.multi = [0,1,1]

  return

end
