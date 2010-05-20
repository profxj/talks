pro fig_clouds

  if not keyword_set(CSZ) then csz = 3.1
  if not keyword_set(lSZ) then lsz = 2.5

  npix = 200
  wave = 1210. + 0.05 * findgen(npix)
  ncloud = [0, 1, 2, 3, 5, $
            10, 100, 1000, 1e5, 1e6]

  ;; Emission Plot
  x_psopen, 'fig_clouds.ps', /maxs
  !p.multi = [0,1,6]
  clr = getcolor(/load)


  xmrg = [8,1]
  ymrg = [0,0]
  yrng=[-0.05, 1.15]
  xrng=[1214., 1217]

  lya = x_setline(1215.6701d)
  N1 = 13.5
  lya.b = 30.

  for qq=0,4 do begin

     xspaces = replicate(' ',30) 
     xtit = ''
     !p.multi = [0,1,6]
      ;; Dummy
;      plot, [0], [0], color=clr.white, background=clr.white, xsty=4, ysty=4

     for ii=0,qq do begin

        if ii EQ 2 then ytit = 'Normalized Flux' else ytit= ''
        if ii EQ qq then begin
           xtit = 'Wavelength (Ang)'
           delvarx, xspaces
        endif
        ;; Create the spectrum
        fx = replicate(1., npix)
        if ncloud[ii] GT 0 then begin
           lya.N = N1 + alog10(ncloud[ii])
           fx = x_voigt(wave, lya, FWHM=3)
        endif

        ;; Plot
        plot, [0], [0], color=clr.lightgray, background=clr.black, charsize=csz,$
              xmargin=xmrg, ymargin=ymrg, ytitle=ytit, $
              xtitle=xtit, yrange=yrng, thick=5, ytickinter=0.5, $
              xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata, $
              xtickn=xspaces

        oplot, wave, fx, color=clr.yellow, thick=9
        if ncloud[ii] LE 10 then cstr = strtrim(round(ncloud[ii]),2) $
        else cstr = '10!u'+strtrim(round(alog10(ncloud[ii])),2)+'!N'
        xyouts, xrng[0]+0.3, 0.2, 'n!dc!N = '+cstr, $
                color=clr.yellow, charsiz=lsz
     endfor
  endfor

  xrng=[1212., 1219]
  for qq=0,4 do begin

     xspaces = replicate(' ',30) 
     xtit = ''
     !p.multi = [0,1,6]
      ;; Dummy
;      plot, [0], [0], color=clr.white, background=clr.white, xsty=4, ysty=4

     for ii=5,5+qq do begin

        if ii EQ 7 then ytit = 'Normalized Flux' else ytit= ''
        if ii EQ (qq+5) then begin
           xtit = 'Wavelength (Ang)'
           delvarx, xspaces
        endif
        ;; Create the spectrum
        fx = replicate(1., npix)
        if ncloud[ii] GT 0 then begin
           lya.N = N1 + alog10(ncloud[ii])
           fx = x_voigt(wave, lya, FWHM=3)
        endif

        ;; Plot
        plot, [0], [0], color=clr.lightgray, background=clr.black, charsize=csz,$
              xmargin=xmrg, ymargin=ymrg, ytitle=ytit, $
              xtitle=xtit, yrange=yrng, thick=5, ytickinter=0.5, $
              xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata, $
              xtickn=xspaces

        oplot, wave, fx, color=clr.yellow, thick=9
        if ncloud[ii] LE 10 then cstr = strtrim(round(ncloud[ii]),2) $
        else cstr = '10!u'+strtrim(round(alog10(ncloud[ii])),2)+'!N'
        xyouts, xrng[0]+0.3, 0.1, 'n!dc!N = '+cstr, $
                color=clr.yellow, charsiz=lsz
     endfor
  endfor

  ;; End
  x_psclose
  !p.multi = [0,1,1]


  return

end
