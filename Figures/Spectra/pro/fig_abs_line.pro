pro fig_abs_line

  if not keyword_set(CSZ) then csz = 2.5
  if not keyword_set(lSZ) then lsz = 2.5

  npix = 200
  wave = 1210. + 0.05 * findgen(npix)
  ncloud = [0, 1, 2, 10, 1e6]

  ;; Emission Plot
  x_psopen, 'fig_abs_line.ps', /maxs
  !p.multi = [0,1,1]
  clr = getcolor(/load)


  xmrg = [7,1]
  ymrg = [15,1]
  yrng=[0.0, 1.15]
  xrng=1215.67 + [-1,1]*2.7

  lya = x_setline(1215.6701d)
  N1 = 13.5
  lya.b = 30.
  ytit = 'Intensity'
  xtit = 'Wavelength (Ang)'

  for qq=0,4 do begin

     ;; Create the spectrum
     fx = replicate(1., npix)
     if ncloud[qq] GT 0 then begin
        lya.N = N1 + alog10(ncloud[qq])
        fx = x_voigt(wave, lya, FWHM=3)
     endif

     ;; Plot
     plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
           xmargin=xmrg, ymargin=ymrg, ytitle=ytit, $
           xtitle=xtit, yrange=yrng, thick=5, ytickinter=1.0, $
           xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata, $
           xtickint=1.

     oplot, wave, fx, color=clr.black, thick=9, psym=10

     if qq GT 0 then begin
        xyouts, 1213.3, 0.30, 'Hydrogen', color=clr.black, charsi=lsz
        xyouts, 1213.3, 0.15, 'Lyman-!9a!X', color=clr.black, charsi=lsz
     endif

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
