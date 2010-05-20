pro fig_ewdefine

  if not keyword_set(CSZ) then csz = 2.7
  if not keyword_set(lSZ) then lsz = 2.5

  npix = 80000L
  wave = 1204. + 0.001 * findgen(npix)
  dwv = wave - shift(wave,1)
  dwv[0] = dwv[1]
  vel = (wave-1215.6701)/1215.6701 * 3e5  ; km/s
  c = x_constants()

  ;; Emission Plot
  x_psopen, 'fig_ewdefine.ps', /maxs
  !p.multi = [0,1,1]
  clr = getcolor(/load)
  aclr = [clr.yellow]


  xmrg = [8,1]
  ymrg = [4,1]
  yrng=[0., 1.10]

  lya = x_setline(1215.6701d)
  NHI_vary = [13.1]
  b_vary = [10.]

  ytit = 'Normalized Flux' 
  xtit = 'Relative Velocity (km/s)'
  xrng=[-100, 100.]

  for qq=0,1 do begin
     ;; Plot
     plot, [0], [0], color=clr.lightgray, background=clr.black, charsize=csz,$
           xmargin=xmrg, ymargin=ymrg, ytitle=ytit, $
           xtitle=xtit, yrange=yrng, thick=5, ytickinter=0.5, $
           xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata
     
     ;; Create the spectrum
     fx = replicate(1., npix)
     lya.N = NHI_vary[0]
     lya.b = b_vary[0]
     fx = x_voigt(wave, lya, FWHM=3)
     
     oplot, vel, fx, color=aclr[0], thick=9, psym=10
        
     ;; Fill
     if qq EQ 1 then begin
        x_curvefill, vel, fx, replicate(1., npix), $
                     color=clr.cyan
        EW = total( (1.-fx) * dwv )
        ;; 
        xyouts, -90, 0.2, 'W!d!9l!X!N = '+$
                string(EW, format='(f4.2)')+' Ang', color=clr.cyan, $
                charsiz=lsz
     endif
  endfor

  ;; End
  x_psclose
  !p.multi = [0,1,1]


  return

end
