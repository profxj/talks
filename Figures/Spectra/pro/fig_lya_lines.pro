pro fig_lya_lines

  if not keyword_set(CSZ) then csz = 3.0
  if not keyword_set(lSZ) then lsz = 3.0

  npix = 800
  wave = 1205. + 0.05 * findgen(npix)
  ncloud = [0.3, 3e3, 1e7]
  lbls = ['Ly!9a!X Forest', 'LLS', 'DLA']

  ;; Emission Plot
  x_psopen, 'fig_lya_lines.ps', /maxs
  !p.multi = [0,1,1]
  clr = getcolor(/load)
  lclr = clr.white


  xmrg = [7,1]
  ymrg = [4,1]
  yrng=[0.0, 1.15]
  xrng=1215.67 + [-1,1]*2.7

  lya = x_setline(1215.6701d)
  N1 = 13.5
  lya.b = 30.
  ytit = 'Intensity'
  xtit = 'Wavelength (Ang)'

  for qq=0,2 do begin

     ;; Create the spectrum
     fx = replicate(1., npix)
     if ncloud[qq] GT 0 then begin
        lya.N = N1 + alog10(ncloud[qq])
        fx = x_voigt(wave, lya, FWHM=3)
     endif
     if qq LT 2 then begin
        xlbl = 1213.1 
        ylbl = 0.26
     endif else begin
        xlbl = 1213.1
        ylbl = 0.8
        xrng=1215.67 + [-1,1]*9.7
     endelse

     ;; Plot
     plot, [0], [0], color=lclr, background=clr.black, charsize=csz,$
           xmargin=xmrg, ymargin=ymrg, ytitle=ytit, $
           xtitle=xtit, yrange=yrng, thick=9, ytickinter=1.0, $
           xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata;, $
;           xtickint=1.

     oplot, wave, fx, color=lclr, thick=9, psym=10

     xyouts, xlbl, ylbl, lbls[qq], color=lclr, charsi=lsz
     xyouts, xlbl, ylbl-0.1, 'HI Ly!9a!X', color=clr.tomato, charsi=lsz
     xyouts, xlbl, ylbl-0.2, 'N!dHI!N = 10!u'+string(lya.N,format='(f4.1)')+'!N cm!u-2!N', $
             color=clr.orange, charsi=lsz

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
