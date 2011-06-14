;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro fig_cos_j1016

  compile_opt strictarr

  ;; Get structure if necessary
  if not keyword_set( PSFILE ) then psfile = 'fig_cos_j1016.ps'
  if not keyword_set( CSZ ) then csz = 1.6
  if not keyword_set( LSZ ) then lsz = 1.9
  if not keyword_set( BLSZ ) then blsz = 2.0

  if not keyword_set( MAXCHI ) then MAXCHI = 1.3
  if not keyword_set( WVMIN ) then wvmin = 3900. 
  if not keyword_set(OMEGA_M) then omega_m = 0.3
  if not keyword_set(H0) then H0 = 72.
  if not keyword_set(ZQSO) then zqso = 3.8
  if not keyword_set(FWHM) then fwhm = 2.


  ;; Read in Data
  flux = x_readspec('~/Dropbox/COS-Halos/Targets/J1016+4706/Data/J1016+4706_nbin6_coadd.fits', $
                    wav=wave, infl=3)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; PLOT
  x_psopen, psfile, /maxs
  clr = getcolor(/load)
  !p.multi=[0,1,1]
  
  xmrg = [8,2]
  ymrg = [15.0,4]

  thk = 5
;  yrng=[-0.02, 3.3]
  yrng=[0.0, 9.]
  xrng=[1150., 1800]

  ;; Start plot
  plot, [0], [0], color=clr.black, background=clr.white, $
        charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, xtitle='Wavelength (Angstroms)', $
        ytitle='Counts (10!u-3!N s!u-1!N)', yrange=yrng, thick=4, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata ;, /ylog
     
     
  oplot, wave, flux*1e3, color=clr.black, thick=3, psym=10

  xyouts, 4450., 1.1, 'J1016+4706 (HST/COS)', color=clr.red, charsi=lsz, align=0.

  ;; Close Ps
  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]

  print, 'fig_sdss:  All done!'
       
  return
end
      
      
