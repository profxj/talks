; Shows the neutral fraction vs. n_H for z=0 EUVB

pro fig_xhi_vs_nh_vs_nhi, CLDY=cldy

  ;; Plotting
  if not keyword_set( PSFILE ) then psfile = 'fig_xhi_vs_nh_vs_nhi.ps'
  if not keyword_set( CSZ ) then csz = 2.1
  if not keyword_set( CSZ2 ) then csz2 = 1.9
  if not keyword_set( LSZ ) then lsz = 2.1

  c= x_constants()

  ;; Cloudy output
  if not keyword_set(CLDY) then $
     cldy = xmrdfits('~/Cloudy/Grid/Output/cgm_z0.2_qg.fits',1)
     ;cldy = xmrdfits('~/Cloudy/Grid/Output/z0_igm_qg.fits',1)
     ;cldy = xmrdfits('~/Dropbox/COS-Halos/lowions/finercldygrid_logn_coshalos.fits',1)

  ;; PLOT
  x_psopen, psfile, /maxs
  !p.multi=[0,1,1]
  clr = getcolor(/load)
  bclr = clr.black
  fclr = clr.lightgray
  
  xtit = 'n!dH!N (cm!u-3!N)'
  ytit = 'x!dHI!N'
          
  xrng = [1e-7, 1.]
  yrng = [1e-6, 1.]
      ;;
  plot, [0.], [0.], color=fclr, background=bclr, charsize=csz,$
        xmargin=[8.0,8], ymargin=[4,4.], xtitle=xtit, $
        /nodata, xrange=xrng, ystyle=9, ytitl=ytit, $
        yrange=yrng, xstyle=9, /xlog, /ylog

  xyouts, 5e-7, 2e-1, 'z=0.2; EUVB only', color=fclr, charsi=lsz

  NHI_val = [15., 17., 19]
  ;NHI_val = [15.]
  nNHI = n_elements(NHI_val)
  NHI_clr = [clr.tomato, clr.cyan, clr.yellow]

  for ss=0L,nNHI-1 do begin
     idx = where(abs(cldy.NHI - NHI_val[ss]) LT 1e-3)
     scldy = cldy[idx]
     srt = sort(scldy.nH)
     scldy =  scldy[srt]
     nH = 10.d^scldy.nH
     xHI = 10.d^scldy.X[1,1]
     oplot, nH, xHI, color=NHI_clr[ss]
     ;; Label
     xyouts, 7e-3, 5e-6*(1.5^ss), 'N!dHI!N = 10!u'+strtrim(round(NHI_val[ss]),2)+'!N cm!u-2!N', $
             color=NHI_clr[ss], charsi=lsz
  endfor

  ;; Label top axis
  xyouts, 10.^(-3.7), 10.^(0.6), '!9dr/r!X', charsiz=csz, color=fclr

  ;; Side axis
  oplot, replicate(xrng[1]*0.9999, 2), yrng, color=fclr, thick=5
  dval = [1,10., 100., 1000., 10000., 1e5, 1e6,1e7]
  NHI = 1d15
  nd = n_elements(dval)
  ;for ii=0,nd-1 do begin
  ;   xHIval = interpol(xHI, NHI/xHI/nH/c.pc, dval[ii]) 
  ;   xyouts, xrng[1]*1.5, xHIval*0.8, '10!u'+string(alog10(dval[ii]),format='(i1)'), $
  ;           color=lclr, charsiz=csz, align=0.0
  ;   oplot, [xrng[1],xrng[1]*0.70], replicate(xHIval,2), color=lclr, thick=5
  ;endfor
  xyouts, 0.97, 0.5, 'size (pc; N!dHI!N = 10!u15!N cm!u-2!N)', color=fclr, charsi=csz, $
          align=0.5, orient=90, /norma

  ;; Top axis
;  rhob = c.rhoc * (0.72)^2 * 0.04  ;; Mass density in baryons g/cm^3
;  nH_b = rhob / (c.mp*1.3)  ;; H per cm^3
  nH_b = c.nhb

  xrho = xrng / nH_b
  axis, xaxis=1, color=fclr, /xlog, charsi=csz, xrang=xrho, xsty=1



  x_psclose
  !p.multi=[0,1,1]
  print, 'fig_beta:  All done!'

  return
end
      
      
