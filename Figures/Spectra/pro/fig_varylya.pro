pro fig_varylya

  if not keyword_set(CSZ) then csz = 2.7
  if not keyword_set(lSZ) then lsz = 2.5

  npix = 8000L
  wave = 1204. + 0.001 * findgen(npix)
  vel = (wave-1215.6701)/1215.6701 * 3e5  ; km/s
  c = x_constants()

  ;; Emission Plot
  x_psopen, 'fig_varylya.ps', /maxs
  !p.multi = [0,1,1]
  clr = getcolor(/load)
  aclr = [clr.yellow, clr.cyan, clr.orange, clr.green, clr.pink]


  xmrg = [8,1]
  ymrg = [4,1]
  yrng=[0., 1.10]

  lya = x_setline(1215.6701d)
  NHI_vary = [13.1, 12.1, 14.1, 16.1, 20.1] 
  nvary = n_elements(NHI_vary)
  b_vary = [10, 1., 30, 100, 300.]
  srt_N = sort(NHI_vary)
  srt_b = sort(b_vary)

  ytit = 'Normalized Flux' 
  xtit = 'Relative Velocity (km/s)'
  ;; Vary N_HI first
  for tt=0,nvary-1 do begin

     if tt EQ (nvary-1) then xrng = [-2700., 2700] else xrng=[-100, 100.]

     ;; Plot
     plot, [0], [0], color=clr.lightgray, background=clr.black, charsize=csz,$
           xmargin=xmrg, ymargin=ymrg, ytitle=ytit, $
           xtitle=xtit, yrange=yrng, thick=5, ytickinter=0.5, $
           xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata

     for qq=0L,tt do begin


        ;; Create the spectrum
        fx = replicate(1., npix)
        lya.N = NHI_vary[qq]
        lya.b = b_vary[0]
        fx = x_voigt(wave, lya, FWHM=3)
        
        oplot, vel, fx, color=aclr[qq], thick=9, psym=10
        
        ;; tau and N_HI
        tau0 = sqrt(!pi) * c.e^2 / c.me / c.c * (10.d^NHI_vary[qq] * lya.f * $
                                                 1215.6701 * 1e-8) / (lya.b*1e5)
        if tau0 LT 1 then tstr = string(tau0,format='(f3.1)')
        if tau0 GE 1 and tau0 LT 100 then tstr = strtrim(round(tau0),2)
        if tau0 Gt 100 then tstr = '10!u'+strtrim(round(alog10(tau0)),2)+'!N'
        yval = 0.38 - 0.08*srt_N[qq]
        xyouts, xrng[0]+0.05*(xrng[1]-xrng[0]), yval, 'log N!dHI!N = '+ $
                string(NHI_vary[qq], format='(f4.1)'), color=aclr[qq], charsiz=lsz, $
                align=0.
        xyouts, xrng[1]-0.2*(xrng[1]-xrng[0]), yval, '!9t!X!d0!N = '+ tstr, $
                color=aclr[qq], charsiz=lsz, align=0.
     endfor
  endfor

  ;; Vary b 
  for tt=0,nvary-1 do begin

     if tt EQ (nvary-1) then xrng = [-500., 500] else xrng=[-100, 100.]
     ;; Plot
     plot, [0], [0], color=clr.lightgray, background=clr.black, charsize=csz,$
           xmargin=xmrg, ymargin=ymrg, ytitle=ytit, $
           xtitle=xtit, yrange=yrng, thick=5, ytickinter=0.5, $
           xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata

     for qq=0L,tt do begin

        ;; Create the spectrum
        fx = replicate(1., npix)
        lya.N = NHI_vary[0]
        lya.b = b_vary[qq]
        fx = x_voigt(wave, lya, FWHM=3)
        
        oplot, vel, fx, color=aclr[qq], thick=9, psym=10
        
        ;; tau and N_HI
        tau0 = sqrt(!pi) * c.e^2 / c.me / c.c * (10.d^lya.N * lya.f * $
                                                 1215.6701 * 1e-8) / (lya.b*1e5)
        if tau0 LT 1 then tstr = string(tau0,format='(f3.1)')
        if tau0 GE 1 and tau0 LT 100 then tstr = strtrim(round(tau0),2)
        if tau0 Gt 100 then tstr = '10!u'+strtrim(round(alog10(tau0)),2)+'!N'
        yval = 0.38 - 0.08*srt_b[qq]
        xyouts, xrng[0]+0.05*(xrng[1]-xrng[0]), yval, 'b = '+ $
                strtrim(round(b_vary[qq]),2)+' km/s', color=aclr[qq], charsiz=lsz, $
                align=0.
        xyouts, xrng[1]-0.2*(xrng[1]-xrng[0]), yval, '!9t!X!d0!N = '+ tstr, $
                color=aclr[qq], charsiz=lsz, align=0.
     endfor
  endfor

  ;; End
  x_psclose
  !p.multi = [0,1,1]


  return

end
