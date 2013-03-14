;; Shows Best estimates of D/H and thereby Omega b
pro fig_mucol_fj0812

  if not keyword_set(CSZ) then csz = 1.6
  if not keyword_set(lSZ) then lsz = 2.0
  if not keyword_set(lSZ2) then lsz2 = 1.5

  ;; Read
  fj0812 = x_readspec('~/Keck/HIRES/RedData/FJ0812+32/FJ0812+32_f.fits',  /struct)
  mucol = x_readspec('~/Keck/HIRES/RedData/FJ0812+32/FJ0812+32_f.fits',  /struct)

  ;; Emission Plot
  x_psopen, 'fig_mucol_fj0812.ps', /maxs
  !p.multi = [0,1,2]
  clr = getcolor(/load)

  ;; Eta formulation from Steidgman 2007
  xmrg = [9,1]
  ymrg = [3,0.1]

  ytit = 'Normalized Flux'
  xtit = 'Wavelength (Ang)'
  yrng=[0., 1.1]
  xrng=[1520., 2030]

  ;; Plot
  for ss=0,1 do begin
     case ss of
        0: begin
           data = fj0812
           zabs = 2.6263
        end
        1: begin
           data = mucol
           zabs = 0.
        end
        else: stop
     endcase
     
     plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
           xmargin=xmrg, ymargin=ymrg, ytitle=ytit, $
           xtitle=xtit, yrange=yrng, thick=5, $
           xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata
     
     ;; Data
     oplot, data.wv/(1+zabs), data.fx, psym=10, color=clr.black, thick=3
     
  endfor
  x_psclose

  !p.multi = [0,1,1]

  return

end
