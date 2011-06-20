pro fig_lls_models

  if not keyword_set(CSZ) then csz = 2.5
  if not keyword_set(lSZ) then lsz = 2.5

  ;; Lines
  nlin = 20
  NHI = 16. + findgen(nlin)*0.2

  ;; Model
  lls_model_fil = getenv('XIDL_DIR')+'/SDSS/LLS/nhi16_19b30.fits'
  nh1 = xmrdfits(lls_model_fil,0)
  npnh = n_elements(nh1)
  wv_mod = 10^(2.7d + dindgen(npnh)*1e-4)

  ;; Emission Plot
  x_psopen, 'fig_lls_models.ps', /maxs
  !p.multi = [0,1,1]
  clr = getcolor(/load)
  lclr = clr.white

  xmrg = [7,1]
  ymrg = [7,1]
  yrng=[0.0, 1.15]
  xrng=[700., 1230]

  ytit = 'Normalized Flux'
  xtit = 'Wavelength (Ang)'

  for qq=0,nlin-1 do begin

     ;; Get the spectrum
     fx = xmrdfits(lls_model_fil, qq)

     ;; Plot
     plot, [0], [0], color=lclr, background=clr.black, charsize=csz,$
           xmargin=xmrg, ymargin=ymrg, ytitle=ytit, $
           xtitle=xtit, yrange=yrng, thick=9, ytickinter=1.0, $
           xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata

     oplot, wv_mod, fx, color=lclr, thick=9, psym=10

     xyouts, 1050., 1.04, 'N!dHI!N = 10!u'+string(NHI[qq], format='(f4.1)')+'!N cm!u-2!N', $
             color=clr.yellow, charsi=lsz, align=0.

;     if ncloud[ii] LE 10 then cstr = strtrim(round(ncloud[ii]),2) $
;     else cstr = '10!u'+strtrim(round(alog10(ncloud[ii])),2)+'!N'
;     xyouts, xrng[0]+0.3, 0.2, 'n!dc!N = '+cstr, $
;             color=clr.yellow, charsiz=lsz
  endfor

  ;; End
  x_psclose
  !p.multi = [0,1,1]


  return

end
