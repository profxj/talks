; Shows the neutral fraction vs. n_H for z=0 EUVB

pro fsi3_vs_u, CLDY=cldy

  ;; Plotting
  if not keyword_set( PSFILE ) then psfile = 'fsi3_vs_u.ps'
  if not keyword_set( CSZ ) then csz = 2.1
  if not keyword_set( CSZ2 ) then csz2 = 1.9
  if not keyword_set( LSZ ) then lsz = 2.1

  c= x_constants()

  ;; Cloudy output
  if not keyword_set(CLDY) then $
     cldy = xmrdfits('~/Cloudy/Grid/Output/cgm_z0.2_qg.fits',1)

  ;; PLOT
  x_psopen, psfile, /maxs
  !p.multi=[0,1,1]
  clr = getcolor(/load)
  ;bclr = clr.black
  ;fclr = clr.lightgray
  bclr = clr.white
  fclr = clr.black
  
  xtit = 'U'
  ytit = 'f(Si!u++!N)'
          
  xrng = [1e-6, 1.]
  yrng = [1e-3, 1.]
      ;;

  plot, [0.], [0.], color=fclr, background=bclr, charsize=csz,$
        xmargin=[8.0,1.5], ymargin=[4,1], xtitle=xtit, $
        /nodata, xrange=xrng, ystyle=1, ytitl=ytit, $
        yrange=yrng, xstyle=1, /xlog, /ylog, ytickformat='x_logticks'
  
  xyouts, 5e-6, 1e2, 'EUVB', color=fclr, charsi=lsz
     
  NHI_val = [15., 17., 18., 19]
                                ;NHI_val = [15.]
  nNHI = n_elements(NHI_val)
  NHI_clr = [clr.red, clr.blue, clr.darkgreen, clr.orange]
  
  for ss=0L,nNHI-1 do begin
     idx = where(abs(cldy.NHI - NHI_val[ss]) LT 1e-3)
     scldy = cldy[idx]
     srt = sort(scldy.U)
     scldy =  scldy[srt]
     
     yplt = 10.d^(scldy.X[14,3])
     oplot, 10.^scldy.U, yplt, color=NHI_clr[ss]
     ;; Label
     xyouts, 1e-3, 2e-3*(1.5^ss), 'N!dHI!N = 10!u'+strtrim(round(NHI_val[ss]),2)+'!N cm!u-2!N', $
             color=NHI_clr[ss], charsi=lsz, align=0.5
  endfor

  x_psclose
  !p.multi=[0,1,1]
  print, 'fsi3_vs_u.:  All done!'

  return
end
      
      
