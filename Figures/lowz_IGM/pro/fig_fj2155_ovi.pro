pro fig_fj2155_ovi

  if not keyword_set(CSZ) then csz = 2.5
  if not keyword_set(lSZ) then lsz = 4.5

  scale = 1d14

  ;; Emission Plot
  x_psopen, 'fig_fj2155_ovi.ps', /maxs
  !p.multi = [0,1,1]
  clr = getcolor(/load)

  tmp = { $
        wvmnx: fltarr(2), $
        ymnx: fltarr(2), $
        zabs: 0.d $
        }
  nsys = 4
  all_sys = replicate(tmp, nsys)
  all_sys[0].zabs = 0.0515
  all_sys[0].wvmnx = [1090.80, 1091.37]
  all_sys[0].ymnx = [-0.3, 8.2]

  all_sys[1].zabs = 0.0808
  all_sys[1].wvmnx = [1115.10, 1115.67]
  all_sys[1].ymnx = [-0.3, 9.6]

  all_sys[2].zabs = 0.1324
  all_sys[2].wvmnx = [1168.13, 1168.95]
  all_sys[2].ymnx = [-0.3, 9]
        
  all_sys[3].zabs = 0.1549
  all_sys[3].wvmnx = [1191.48, 1192.09]
  all_sys[3].ymnx = [-0.3, 8.05]
        
  ;; Read in data
  fx1 = x_readspec('/Users/xavier/MLSS/data/PHL1811/STIS/E140M/PHL1811_ST_E140M_F.fits',wav=wave1)
  fx2 = x_readspec('/Users/xavier/MLSS/data/PHL1811/FUSE/PHL1811lif2a.fits', $
                   wav=wave2, infl=3)
  fx1 = fx1 * scale
  fx2 = fx2 * 1d14


  ;; Plot

  xmrg = [7,1]
  ymrg = [4,1]

  ytit = 'Relative Flux'
  xtit = 'Wavelength (Ang)'

;  for qq=0,0 do begin
  for qq=0,nsys-1 do begin

     if qq LT 2 then begin
        fx = fx2
        wave = wave2
     endif else begin
        fx = fx1
        wave = wave1
     endelse
     ;;
     yrng = all_sys[qq].ymnx
     xrng = all_sys[qq].wvmnx

     ;; Plot
     plot, [0], [0], color=clr.white, background=clr.black, charsize=csz,$
           xmargin=xmrg, ymargin=ymrg, ytitle=ytit, $
           xtitle=xtit, yrange=yrng, thick=8, $
           xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata

     oplot, xrng, [0., 0.], color=clr.tomato, linesty=1
     oplot, wave, fx, color=clr.white, thick=9, psym=10

     xxy = all_sys[qq].wvmnx[0] + (all_sys[qq].wvmnx[1]-all_sys[qq].wvmnx[0])*0.1
     yxy = all_sys[qq].ymnx[0] + (all_sys[qq].ymnx[1]-all_sys[qq].ymnx[0])*0.1
     yxy2 = all_sys[qq].ymnx[0] + (all_sys[qq].ymnx[1]-all_sys[qq].ymnx[0])*0.20
     xyouts, xxy, yxy2, 'OVI 1031', color=clr.yellow, charsi=lsz
     xyouts, xxy, yxy, 'z='+string(all_sys[qq].zabs,format='(f6.4)'), color=clr.yellow, charsi=lsz

  endfor

  ;; End
  x_psclose
  !p.multi = [0,1,1]


  return

end
