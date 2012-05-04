; Shows the neutral fraction vs. n_H for z=0 EUVB

pro euvb_si3_si2_vs_u, CLDY=cldy

  ;; Plotting
  if not keyword_set( PSFILE ) then psfile = 'euvb_si3_si2_vs_u.ps'
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
  
  xtit = 'U'
  ytit = 'Si!u++!N / Si!u+!N'
          
  xrng = [1e-6, 1.]
  yrng = [1e-3, 1e3]
      ;;

  for qq=0,2 do begin
     plot, [0.], [0.], color=fclr, background=bclr, charsize=csz,$
           xmargin=[8.0,1.5], ymargin=[4,1], xtitle=xtit, $
           /nodata, xrange=xrng, ystyle=1, ytitl=ytit, $
           yrange=yrng, xstyle=1, /xlog, /ylog, ytickformat='x_logticks'
     
     xyouts, 5e-6, 1e2, 'EUVB', color=fclr, charsi=lsz
     
     NHI_val = [15., 17., 18., 19]
     ;NHI_val = [15.]
     nNHI = n_elements(NHI_val)
     NHI_clr = [clr.tomato, clr.green, clr.cyan, clr.yellow]
     
     for ss=0L,nNHI-1 do begin
        idx = where(abs(cldy.NHI - NHI_val[ss]) LT 1e-3)
        scldy = cldy[idx]
        srt = sort(scldy.U)
        scldy =  scldy[srt]
        
        yplt = 10.d^(scldy.X[14,3] - scldy.X[14,2])
        oplot, 10.^scldy.U, yplt, color=NHI_clr[ss]
        ;; Label
        xyouts, 2e-2, 2e-3*(2.5^ss), 'N!dHI!N = 10!u'+strtrim(round(NHI_val[ss]),2)+'!N cm!u-2!N', $
                color=NHI_clr[ss], charsi=lsz
     endfor

     ;; Lines
     if qq GT 0 then begin
        oplot, xrng, [0.5, 0.5], color=clr.gray, linesty=2, thick=4
        oplot, xrng, replicate(30., 2), color=clr.gray, linesty=2, thick=4
     endif
     if qq GT 1 then begin
        oplot, replicate(3e-4,2), yrng, color=clr.gray, linesty=2, thick=4
        oplot, replicate(1e-2,2), yrng, color=clr.gray, linesty=2, thick=4
     endif
  endfor

  x_psclose
  !p.multi=[0,1,1]
  print, 'euvb_si3.:  All done!'

  return
end
      
      
