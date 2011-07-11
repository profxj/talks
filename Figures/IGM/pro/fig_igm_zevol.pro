;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro fig_igm_zevol

  compile_opt strictarr

  ;; Get structure if necessary
  if not keyword_set( PSFILE ) then psfile = 'fig_igm_zevol.ps'
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
  
  xrng=[900.0, 1215.]
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Low z first

  xmrg = [8,2]
  ymrg = [3.2,0]

  thk = 5

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; H1821
  yrng=[-0.05, 1.0]
  plot, [0], [0], color=lclr, background=clr.white, $
        charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, xtitle='Rest Wavelength (Angstroms)', $
        ytitle='Flux (10!u-13!N cgs)', yrange=yrng, thick=4, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata ;, /ylog

  oplot, xrng, [0,0], color=clr.gray, linesty=2

  ;; Blend
  x_curvefill, [1301.1, 1304.9], replicate(yrng[0],2), replicate(yrng[1],2), $
               color=clr.gray, /line_fil, orient=45
  ;; Read in Data
  zqso = 0.297
  fx = x_readspec('~/MLSS/data/H1821+643/STIS/E140M/H1821+643__E140M_F.fits',$
                      wav=wav, /auto, sig=sig)
  fx = fx * 1d13
     
  oplot, wav/(1+zqso), fx, color=lclr, thick=3, psym=10
  xlbl = 910.
  xyouts, xlbl, 0.9, 'H1821+640', color=lclr, charsi=lsz, align=0.
  xyouts, xlbl, 0.82, 'z!dQSO!N = 0.297', color=lclr, charsi=lsz, align=0.

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; PSS0209+0517
  yrng=[-0.1, 5.6]
  plot, [0], [0], color=lclr, background=clr.white, $
        charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, xtitle='Rest Wavelength (Angstroms)', $
        ytitle='Flux (cgs)', yrange=yrng, thick=4, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata ;, /ylog

  oplot, xrng, [0,0], color=clr.gray, linesty=2

  ;; Blend
  x_curvefill, [1301.1, 1304.9], replicate(yrng[0],2), replicate(yrng[1],2), $
               color=clr.gray, /line_fil, orient=45
  ;; Read in Data
  zqso = 4.17 
  fx = x_readspec('~/Keck/ESI/RedData/PSS0209+0517/PSS0209+0517_xF.fits', $
                  wav=wav, /auto, sig=sig)
     
  oplot, wav/(1+zqso), fx, color=lclr, thick=3, psym=10
  xyouts, xlbl, 5.1, 'PSS0209+0517', color=lclr, charsi=lsz, align=0.
  xyouts, xlbl, 4.6, 'z!dQSO!N = 4.17', color=lclr, charsi=lsz, align=0.

  ;; Close Ps
  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]

  print, 'fig_sdss:  All done!'
       
  return
end
      
      
