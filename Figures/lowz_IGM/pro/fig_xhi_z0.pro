; Shows the neutral fraction vs. n_H for z=0 EUVB

pro fig_xhi_z0

  ;; Plotting
  if not keyword_set( PSFILE ) then psfile = 'fig_xhi_z0.ps'
  if not keyword_set( CSZ ) then csz = 2.1
  if not keyword_set( CSZ2 ) then csz2 = 1.9
  if not keyword_set( LSZ ) then lsz = 2.1

  c= x_constants()

  ;; Cloudy output
  cldy = xmrdfits('~/Cloudy/Grid/Output/z0_igm_qg.fits',1)
  srt = sort(cldy.nH)
  cldy =  cldy[srt]
  nH = 10.d^cldy.nH
  xHI = 10.d^cldy.X[1,1]

  ;; PLOT
  x_psopen, psfile, /maxs
  !p.multi=[0,1,1]
  clr = getcolor(/load)
  lclr = clr.black
  
  xtit = 'n!dH!N (cm!u-3!N)'
  ytit = 'x!dHI!N'
          
  xrng = [1e-7, 1.]
  yrng = [1e-6, 1.]
      ;;
  plot, [0.], [0.], color=lclr, background=clr.white, charsize=csz,$
        xmargin=[8.0,8], ymargin=[4,4.], xtitle=xtit, $
        /nodata, xrange=xrng, ystyle=9, ytitl=ytit, $
        yrange=yrng, xstyle=9, /xlog, /ylog

  xyouts, 8e-7, 2e-1, 'z=0; EUVB', color=lclr, charsi=lsz
  oplot, nH, xHI, color=clr.yellow

  ;; Label top axis
  xyouts, 10.^(-3.7), 10.^(0.6), '!9dr/r!X', charsiz=csz, color=lclr

  ;; Side axis
  oplot, replicate(xrng[1]*0.9999, 2), yrng, color=lclr, thick=5
  dval = [1,10., 100., 1000., 10000., 1e5, 1e6,1e7]
  NHI = 1d15
  nd = n_elements(dval)
  for ii=0,nd-1 do begin
     xHIval = interpol(xHI, NHI/xHI/nH/c.pc, dval[ii]) 
     xyouts, xrng[1]*1.5, xHIval*0.8, '10!u'+string(alog10(dval[ii]),format='(i1)'), $
             color=lclr, charsiz=csz, align=0.0
     oplot, [xrng[1],xrng[1]*0.70], replicate(xHIval,2), color=lclr, thick=5
  endfor
  xyouts, 0.97, 0.5, 'size (pc; N!dHI!N = 10!u15!N cm!u-2!N)', color=lclr, charsi=csz, $
          align=0.5, orient=90, /norma

  ;; Top axis
;  rhob = c.rhoc * (0.72)^2 * 0.04  ;; Mass density in baryons g/cm^3
;  nH_b = rhob / (c.mp*1.3)  ;; H per cm^3
  nH_b = c.nhb

  xrho = xrng / nH_b
  axis, xaxis=1, color=lclr, /xlog, charsi=csz, xrang=xrho, xsty=1



  x_psclose
  !p.multi=[0,1,1]
  print, 'fig_beta:  All done!'

  return
end
      
      
