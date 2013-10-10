;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro cgm_comp_ovi_hiloz

  compile_opt strictarr

  ;; Get structure if necessary
  if not keyword_set( PSFILE ) then psfile = 'cgm_comp_ovi_hiloz.ps'
  if not keyword_set( CSZ ) then csz = 1.5
  if not keyword_set( LSZ ) then lsz = 1.5
  if not keyword_set( BLSZ ) then blsz = 2.0

  if not keyword_set( MAXCHI ) then MAXCHI = 1.3
  if not keyword_set( WVMIN ) then wvmin = 3900. 
  if not keyword_set(OMEGA_M) then omega_m = 0.3
  if not keyword_set(H0) then H0 = 72.
  if not keyword_set(ZQSO) then zqso = 3.8
  if not keyword_set(FWHM) then fwhm = 2.


  ;; Read in Data
  fx_hires = x_readspec('~/Dropbox/koa/Analysis/Q1700+6416/Q1700+6416a_f.fits', $
                        wav=wv_hires)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; PLOT
  x_psopen, psfile, /portrait
  clr = getcolor(/load)
  !p.multi=[0,1,2]
  fclr=clr.white
  
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; HS1700 first
  xmrg = [8,1]
  ymrg = [4.0,0.5]
  xrng=[3540, 3571.]

  thk = 5
  yrng=[-0.05, 1.13]

  ;; Start plot
  plot, [0], [0], color=fclr, background=clr.white, $
        charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, xtitle='Wavelength (Angstroms)', $
        ytitle='Relative Intensity', yrange=yrng, thick=4, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata ;, /ylog
     
     
  oplot, wv_hires, fx_hires, color=fclr, thick=4, psym=10
  xyouts, 3541., 0.1, 'HS1700+6416', color=clr.yellow, charsi=lsz, align=0.
  xyouts, 3541., 0.03, 'Simcoe+06', color=clr.yellow, charsi=lsz, align=0.
  oplot, xrng, [0., 0.], color=clr.gray, linesty=1, thick=2

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; J1009

  fx_cos = x_readspec( $
           '~/Dropbox/cos-halos/Targets/J1009+0713/Data/J1009+0713_nbin3_coadd.fits', $
           wav=wv_cos, inflg=3)
  xrng=[1265.2, 1275.16]

  thk = 5
  yrng=[-0.0005, 0.0064]

  ;; Start plot
  plot, [0], [0], color=fclr, background=clr.white, $
        charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, xtitle='Wavelength (Angstroms)', $
        ytitle='Relative Intensity', yrange=yrng, thick=4, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata ;, /ylog
     
     
  oplot, wv_cos, fx_cos, color=fclr, thick=4, psym=10

  ;; Label
  xyouts, 1266., 0.006, 'J1009+0713', color=clr.yellow, charsi=lsz, align=0.
  xyouts, 1266., 0.0055, 'Tumlinson+10', color=clr.yellow, charsi=lsz, align=0.
  oplot, xrng, [0., 0.], color=clr.gray, linesty=1, thick=2

  ;; Close Ps
  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]

  print, 'fig_sdss:  All done!'
       
  return
end
      
      
