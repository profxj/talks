; Shows the neutral fraction vs. n_H for z=0 EUVB

pro nh_vs_nhi_si3, CLDY=cldy

  ;; Plotting
  if not keyword_set( PSFILE ) then psfile = 'nh_vs_nhi_si3.ps'
  if not keyword_set( CSZ ) then csz = 2.1
  if not keyword_set( CSZ2 ) then csz2 = 1.9
  if not keyword_set( LSZ ) then lsz = 2.1

  if not keyword_set( Si3val) then si3val = 13.5

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
  
  xtit = 'N!dHI!N'
  ytit = 'N!dH!N'
          
  xrng = [15., 20.]
  yrng = [15., 25]
      ;;

  plot, [0.], [0.], color=fclr, background=bclr, charsize=csz,$
        xmargin=[8.0,1.5], ymargin=[4,1], xtitle=xtit, $
        /nodata, xrange=xrng, ystyle=1, ytitl=ytit, $
        yrange=yrng, xstyle=1, /xlog, /ylog, ytickformat='x_logticks'
  
  xyouts, 5e-6, 1e2, 'EUVB', color=fclr, charsi=lsz
     
  uni_nhi = cldy[uniq(cldy.NHI, sort(cldy.NHI))].NHI
  npts = n_elements(uni_NHI)
  
  for ss=0L,npts-1 do begin
     ;; Grab all at N_HI
     idx = where(abs(cldy.NHI - uni_nhi[ss]) LT 1e-3)
     scldy = cldy[idx]
     srt = sort(scldy.U)
     scldy =  scldy[srt]
     
     NH = uni_nhi[ss] - scldy.X[1,1]
     logf = scldy.X[14,3]
     NSi3 = NH - 12 + 7.5 + logf  ;; Assumes Solar
     mn = min(abs(NSi3-Si3val), imn)
     ;NHval = interpol(NH, NSi3, 13.5)
     if mn LT 0.1 then psym = 1 else psym=2

     oplot, [uni_nhi[ss]], [NH[imn]], psym=psym, color=clr.blue
  endfor

  x_psclose
  !p.multi=[0,1,1]
  print, 'fsi3_vs_u.:  All done!'

  return
end
      
      
