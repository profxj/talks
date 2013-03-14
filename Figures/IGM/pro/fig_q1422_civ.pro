;; Shows Best estimates of D/H and thereby Omega b
pro fig_q1422_civ

  if not keyword_set(CSZ) then csz = 1.6
  if not keyword_set(lSZ) then lsz = 1.6
  if not keyword_set(lSZ2) then lsz2 = 1.5

  ;; Read
  q1422 = x_readspec('~/Keck/HIRES/RedData/Q1422+2309/Q1422+2309.fits', $
                     infl=2, /struct)

  ;; Fix sky
  sky = where( q1422.wv GT 6301.4 and q1422.wv LT 6302.9)
  q1422.fx[sky] = 1662.
  
  zciv = [2.6828, 2.72011, 2.74896, 2.77163, 2.78618, 2.79674, 2.81837, $
          2.89550, 2.94754, 2.96197, 2.97135, 2.97614, 3.00484, 2.99921, $
          3.03504, 3.06433, 3.08666, 3.09003, 3.13440, 3.13708, 3.15600, $
          3.19144, 3.23329, 3.24048, 3.25704, 3.26570, 3.27595, 3.31800, $
          3.33410, 3.38167, 3.39740, 3.41146, 3.51471]
  nciv = n_elements(zciv)

  ;; Emission Plot
  x_psopen, 'fig_q1422_civ.ps', /maxs
  !p.multi = [0,1,1]
  clr = getcolor(/load)

  ;; Eta formulation from Steidgman 2007
  xmrg = [9,7]
  ymrg = [7,17]
  yrng=[0., 2787]
  xrng=[5684, 6858.]

  ytit = 'Relative Flux'
  xtit = 'Wavelength (Ang)'

  asz = 1.9

  ;; Plot
  for ss=0,1 do begin
     plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
           xmargin=xmrg, ymargin=ymrg, ytitle=ytit, $
           xtitle=xtit, yrange=yrng, thick=5, $
           xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata
     
     ;; Data
     oplot, q1422.wv, q1422.fx, psym=10, color=clr.black, thick=3
     
     ;; Mark CIV
     plotsym, 1, 1.8, thick=6
     if ss GT 0 then begin
        for kk=0L,nciv-1 do begin
           wvp = 1549.*(1+zciv[kk]) 
           pix = where( abs( q1422.wv-wvp ) LT 2)
           oplot, [1549.*(1+zciv[kk])], [median(q1422.fx[pix])*1.15], $
                  psym=8, color=clr.blue
        endfor
        xyouts, 6090., 2300., 'CIV', color=clr.blue, charsi=lsz, align=0.5
     endif
     
  endfor
  x_psclose
  !p.multi = [0,1,1]

  return

end
