;; 
pro fig_pcygni, psfile

  lsz = 1.5
  csz = 2.5

  xrng = [2793.5, 2797.5]
  yrng = [-0.05, 1.25]

  wave = x_mkwvarray(2., xrng[0], xrng[1], NPIX=npix)
  fx = fltarr(npix)
  mgii_lin = x_setline(2796.352d)

  ;; 
  if not keyword_set(PSFILE) then psfile = 'fig_pcygni.ps'

  x_psopen, psfile, /maxs
  !p.multi=[0,1,1]
  clr = getcolor(/load)

  for qq=0,3 do begin 

     case qq of 
        0: fx[*] = 1.
        1: begin ;; v = 0 km/s
           mgii_lin.N = 14.
           mgii_lin.b = 20.
           fx = x_voigt(wave, mgii_lin, FWHM=4)
           ew = total(1.-fx)*(wave[1]-wave[0])
        end
        2: begin ;; v = -200 km/s
           mgii_lin.zabs = -200./3e5 
           fx = x_voigt(wave, mgii_lin, FWHM=4)
        end
        3: begin ;; v = -200 km/s
           gdp = where(fx LT 0.99, ngd)
           mn = min(fx, imn)
           mn2 = min(abs(wave-2796.352), imn2)
           fx[imn2-ngd/2+lindgen(ngd)] = fx[imn2-ngd/2+lindgen(ngd)] + (1-fx[gdp])
           yrng[1] = 2.1
        end
        else: stop
     endcase
     
     ;;;;;;;;;;;;;;;;;;;;;;;
     ;; Plot
     plot, [0],  [0],  color=clr.black, $
           background=clr.black, charsize=csz,$
           xmargin=[7,2.5], ymargin=[3.5,4.5], xtitle='Wavelength (Ang)', $
           ytitle='Normalized Flux', /nodata, xthick=7, ythick=7, xstyle=9, ystyle=1, $
           yr=yrng, xr=xrng
     
     if qq GT 0 then oplot, replicate(2796.352,2), $
                            yrng, color=clr.cyan, thick=5, linesty=2
     oplot, xrng, [0.,0], color=clr.yellow, thick=5, linesty=1
     oplot, wave, fx, color=clr.black, thick=7, psym=10

     xrng2 = (xrng/2796.352 - 1)*3e5
     axis, xaxis=1, charsiz=csz, xsty=1, xrang=xrng2, xtitl='!9d!Xv (km/s)'
     
     
  endfor
  x_psclose
  !p.multi=[0,1,1]

  return
end
