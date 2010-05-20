;; 
pro fig_civdoublet, psfile

  lsz = 1.5
  csz = 2.5

 ;; 
  if not keyword_set(PSFILE) then psfile = 'fig_civdoublet.ps'


  x_psopen, psfile, /maxs
  !p.multi=[0,1,1]
  clr = getcolor(/load)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; GRB 050401
  datfil = '/u/xavier/Keck/HIRES/RedData/FJ0812+32/FJ0812+32_f.fits'
  fx = x_readspec(datfil, wav=wave)
  xrng = [5580., 5640]
  yrng = [-0.05, 1.25]

  ;; Plot
  plot, [0],  [0],  color=clr.lightgray, $
    background=clr.black, charsize=csz,$
    xmargin=[7,2.5], ymargin=[3.5,0.5], xtitle='Wavelength (Ang)', $
    ytitle='Normalized Flux', /nodata, xthick=5, ythick=5, xstyle=1, ystyle=1, $
    yr=yrng, xr=xrng

  oplot, wave, fx, color=clr.lightgray, thick=5, psym=10

  ;; Label
  oplot, [5592., 5592., 5609, 5609], [1.1, 1.2, 1.2, 1.1], $
         color=clr.cyan, thick=7
  oplot, [5612.5, 5612.5, 5627.6, 5627.6], [1.1, 1.2, 1.2, 1.1], $
         color=clr.cyan, thick=7

  x_psclose
  !p.multi=[0,1,1]

  return
end
