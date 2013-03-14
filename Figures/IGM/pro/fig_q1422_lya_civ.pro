;; Shows Best estimates of D/H and thereby Omega b
pro fig_q1422_lya_civ

  if not keyword_set(CSZ) then csz = 1.6
  if not keyword_set(lSZ) then lsz = 2.0
  if not keyword_set(lSZ2) then lsz2 = 1.5

  ;; Read
  q1422 = x_readspec('~/Keck/HIRES/RedData/Q1422+2309/Q1422+2309N.fits', $
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

  zplt = [2.7, 3.17]

  ;; Emission Plot
  x_psopen, 'fig_q1422_lya_civ.ps', /maxs
  !p.multi = [0,1,2]
  clr = getcolor(/load)

  ;; Eta formulation from Steidgman 2007
  xmrg = [9,1]
  ymrg = [3,0.1]



  ytit = 'Normalized Flux'
  xtit = 'Wavelength (Ang)'
  yrng=[0., 1.1]

  asz = 1.9

  ;; Plot
  for ss=0,1 do begin
     case ss of
        0: wrest = 1548.195
        1: wrest = 1215.6701d
        else: stop
     endcase
     
     xrng=(1+zplt)*wrest
     plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
           xmargin=xmrg, ymargin=ymrg, ytitle=ytit, $
           xtitle=xtit, yrange=yrng, thick=5, $
           xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata
     
     ;; Data
     oplot, q1422.wv, q1422.fx, psym=10, color=clr.black, thick=3
     
     ;; Mark CIV
     plotsym, 1, 1.8, thick=6
     if ss EQ 1 then pclr = clr.red else pclr = clr.blue
     for kk=0L,nciv-1 do begin
        wvp = wrest*(1+zciv[kk]) 
        pix = where( abs( q1422.wv-wvp ) LT 2)
        oplot, [wrest*(1+zciv[kk])], [1.08], psym=8, color=pclr
     endfor

     ;; Label
     yplt = 0.2
     if ss EQ 1 then begin
        xplt = 6090.*1215.67/1548.19
        x_curvefill, xplt+[-40, 40], [0.14, 0.14], [0.32, 0.32], color=clr.white
        xyouts, 6090.*1215.67/1548.19, yplt, 'HI Ly!9a!X', color=clr.red, charsi=lsz, align=0.5 
     endif else $
        xyouts, 6090., yplt, 'CIV 1548', color=clr.blue, charsi=lsz, align=0.5
     
  endfor
  x_psclose
  !p.multi = [0,1,1]

  return

end
