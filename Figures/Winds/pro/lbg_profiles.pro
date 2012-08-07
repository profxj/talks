pro lbg_profiles

  if not keyword_set(CSZ) then csz = 1.7
  if not keyword_set(lSZ) then lsz = 1.5
  if not keyword_set(lSZ2) then lsz2 = 1.9

  npix = 200L
  wave = 1200. + 0.3 * findgen(npix)
  vel = (wave-1215.6701)/1215.6701 * 3e5  ; km/s
  c = x_constants()

  ;; File
  x_psopen, 'lbg_profiles.ps', /portrait
  !p.multi = [0,1,2]
  clr = getcolor(/load)
  aclr = [clr.red, clr.black]

  xmrg = [8,1]
  ymrg = [4,1]
  yrng=[0., 1.10]

  ;; HI
  lya = x_setline(1215.6701d)
  NHI_vary = [16.0, 18.5]
  lbls = ['High N!dHI!N', 'Low N!dHI!N']
  nvary = n_elements(NHI_vary)
  b_vary = [55, 30]
  Cf_vary = [1., 0.9]
  srt_N = sort(NHI_vary)
  srt_b = sort(b_vary)

  ytit = 'Normalized Flux' 
  xtit = 'Relative Velocity (km/s)'

  ;; Plot
  xrng=[-500, 500.]
  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, ytitle=ytit, $
        xtitle=xtit, yrange=yrng, thick=5, ytickinter=0.5, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata

  ;; Loop
  for tt=0,nvary-1 do begin

     ;; Create the spectrum
     fx = replicate(1., npix)
     lya.N = NHI_vary[tt]
     lya.b = b_vary[tt]
     fx = x_voigt(wave, lya, FWHM=2, COVERING=Cf_vary[tt])
        
     oplot, vel, fx, color=aclr[tt], thick=9, psym=10
        
     yval = 0.38 - 0.08*srt_N[tt]
     xyouts, xrng[0]+0.05*(xrng[1]-xrng[0]), yval, lbls[tt], $
             color=aclr[tt], charsiz=lsz,  align=0.
  endfor
  
  xyouts, xrng[1]-0.20*(xrng[1]-xrng[0]), 0.15, 'Ly!9a!X', $
          color=clr.black, charsiz=lsz2,  align=0.

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; CII
  npix = 200L
  wrest = 1334.5323d
  wave = wrest - 15 + 0.3 * findgen(npix)
  vel = (wave-wrest)/wrest * 3e5  ; km/s

  line = x_setline(wrest)
  metals = [0., -2.]
  NHI_vary = [16.0, 18.5]  + metals - 3.6 + 2.
  nvary = n_elements(NHI_vary)
  b_vary = [55, 30]
  Cf_vary = [1., 0.9]
  srt_N = sort(NHI_vary)
  srt_b = sort(b_vary)

  ytit = 'Normalized Flux' 
  xtit = 'Relative Velocity (km/s)'

  ;; Plot
  xrng=[-500, 500.]
  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, ytitle=ytit, $
        xtitle=xtit, yrange=yrng, thick=5, ytickinter=0.5, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata

  ;; Loop
  for tt=0,nvary-1 do begin

     ;; Create the spectrum
     fx = replicate(1., npix)
     line.N = NHI_vary[tt]
     line.b = b_vary[tt]
     fx = x_voigt(wave, line, FWHM=2, COVERING=Cf_vary[tt])
        
     oplot, vel, fx, color=aclr[tt], thick=9, psym=10
        
     yval = 0.38 - 0.08*srt_N[tt]
     xyouts, xrng[0]+0.05*(xrng[1]-xrng[0]), yval, lbls[tt], $
             color=aclr[tt], charsiz=lsz,  align=0.
  endfor

  xyouts, xrng[1]-0.20*(xrng[1]-xrng[0]), 0.15, 'CII 1334', $
          color=clr.black, charsiz=lsz2,  align=0.

  ;; End
  x_psclose
  !p.multi = [0,1,1]


  return

end
