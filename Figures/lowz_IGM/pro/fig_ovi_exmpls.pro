;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro fig_ovi_exmpls

  compile_opt strictarr

  ;; Get structure if necessary
  if not keyword_set( PSFILE ) then psfile = 'fig_ovi_exmpls.ps'
  if not keyword_set( CSZ ) then csz = 1.9
  if not keyword_set( LSZ ) then lsz = 1.9
  if not keyword_set( BLSZ ) then blsz = 2.0

  if not keyword_set( MAXCHI ) then MAXCHI = 1.3
  if not keyword_set( WVMIN ) then wvmin = 3900. 
  if not keyword_set(OMEGA_M) then omega_m = 0.3
  if not keyword_set(H0) then H0 = 72.
  if not keyword_set(ZQSO) then zqso = 3.8
  if not keyword_set(FWHM) then fwhm = 2.


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; PLOT
  x_psopen, psfile, /maxs
  clr = getcolor(/load)
  lclr = clr.white
  !p.multi=[0,1,2]
  
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; PKS0312-770 first

  xmrg = [8,2]
  ymrg = [3.2,0]
  xrng=[1231., 1249.3]
  yrng=[-0.1, 1.9]

  thk = 5

  ;; Plot
  plot, [0], [0], color=lclr, background=clr.white, $
        charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, xtitle='Wavelength (Angstroms)', $
        ytitle='Flux (10!u-14!N cgs)', yrange=yrng, thick=4, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata ;, /ylog

  oplot, xrng, [0,0], color=clr.gray, linesty=2

  ;; Read in Data
  zovi = 0.202660
  fx = x_readspec('~/MLSS/data/B0312-770/STIS/E140M/B0312-770__E140M_F.fits',$
                      wav=wav, /auto, sig=sig)
  fx = fx * 1d14
     
  oplot, wav, fx, color=lclr, thick=4, psym=10
  xyouts, 1232., 1.7, 'PKS0312-770', color=lclr, charsi=lsz, align=0.
  xyouts, 1232., 1.55, 'z!dOVI!N = 0.203', color=lclr, charsi=lsz, align=0.

  xyouts, 1233.2, 1.3, 'HI Ly!9b!X', color=clr.orange, charsi=lsz, align=0.5
  xyouts, 1241.1, 1.69, 'OVI', color=clr.cyan, charsi=lsz, align=0.5
  xyouts, 1241.1, 1.5, '1031', color=clr.cyan, charsi=lsz, align=0.5
  xyouts, 1248.0, 1.69, 'OVI', color=clr.cyan, charsi=lsz, align=0.5
  xyouts, 1248.0, 1.5, '1037', color=clr.cyan, charsi=lsz, align=0.5
  xyouts, 1246.4, 1.69, 'CII', color=clr.yellow, charsi=lsz, align=0.5
  xyouts, 1246.4, 1.5, '1036', color=clr.yellow, charsi=lsz, align=0.5

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; H1821
  xrng=[1298.1, 1315.2]
  yrng=[-0.1, 5.5]

  ;; Plot
  plot, [0], [0], color=lclr, background=clr.white, $
        charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, xtitle='Wavelength (Angstroms)', $
        ytitle='Flux (10!u-14!N cgs)', yrange=yrng, thick=4, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata ;, /ylog

  oplot, xrng, [0,0], color=clr.gray, linesty=2

  ;; Blend
  x_curvefill, [1301.1, 1304.9], replicate(yrng[0],2), replicate(yrng[1],2), $
               color=clr.gray, /line_fil, orient=45
  ;; Read in Data
  zovi = 0.26656
  fx = x_readspec('~/MLSS/data/H1821+643/STIS/E140M/H1821+643__E140M_F.fits',$
                      wav=wav, /auto, sig=sig)
  fx = fx * 1d14
     
  oplot, wav, fx, color=lclr, thick=4, psym=10
  xyouts, 1311., 1.0, 'H1821+640', color=lclr, charsi=lsz, align=0.
  xyouts, 1311., 0.5, 'z!dOVI!N = 0.267', color=lclr, charsi=lsz, align=0.

  xyouts, 1299.1, 2.6, 'HI', color=clr.orange, charsi=lsz, align=0.5
  xyouts, 1299.1, 2.0, 'Ly!9b!X', color=clr.orange, charsi=lsz, align=0.5
  xyouts, 1307.0, 1.9, 'OVI', color=clr.cyan, charsi=lsz, align=0.5
  xyouts, 1307.0, 1.3, '1031', color=clr.cyan, charsi=lsz, align=0.5
  xyouts, 1314.2, 2.90, 'OVI', color=clr.cyan, charsi=lsz, align=0.5
  xyouts, 1314.2, 2.3, '1037', color=clr.cyan, charsi=lsz, align=0.5
;  xyouts, 1246.4, 1.69, 'CII', color=clr.yellow, charsi=lsz, align=0.5
;  xyouts, 1246.4, 1.5, '1036', color=clr.yellow, charsi=lsz, align=0.5

  ;; Close Ps
  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]

  print, 'fig_sdss:  All done!'
       
  return
end
      
      
