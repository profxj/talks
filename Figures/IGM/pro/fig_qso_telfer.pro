;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro fig_qso_telfer

  compile_opt strictarr

  ;; Get structure if necessary
  if not keyword_set( PSFILE ) then psfile = 'fig_qso_telfer.ps'
  if not keyword_set( CSZ ) then csz = 2.3
  if not keyword_set( LSZ ) then lsz = 1.6
  if not keyword_set( BLSZ ) then blsz = 2.0

  if not keyword_set( MAXCHI ) then MAXCHI = 1.3
  if not keyword_set( WVMIN ) then wvmin = 3900. 
  if not keyword_set(OMEGA_M) then omega_m = 0.3
  if not keyword_set(H0) then H0 = 72.
  if not keyword_set(ZQSO) then zqso = 3.8
  if not keyword_set(FWHM) then fwhm = 2.

  c = x_constants()
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; PLOT
  x_psopen, psfile, /maxs
  clr = getcolor(/load)
  !p.multi=[0,1,1]
  
  ;; FIRST PLOT ::  Data and best model
  xmrg = [8,2]
  ymrg = [4.0,4]
  xrng1=[800., 1500.]

  thk = 5
  yrng=[-0.02, 4.0]

  ;; Start plot
  plot, [0], [0], color=clr.black, background=clr.white, $
        charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, xtitle='Rest Wavelength (Angstroms)', $
        ytitle='Relative Flux', yrange=yrng, thick=4, $
        xrange=xrng1, ystyle=1, xstyle=1, psym=1, /nodata ;, /ylog

  ;; Read in Telfer
  readcol, getenv('LLSPAP')+'/taueff/Analysis/Telfer/hst_comp01_rq.asc', $ 
           wave, flux, format='F,F'
  gd = where(wave GT xrng1[0] and wave LT xrng1[1])

  oplot, wave[gd], flux[gd], color=clr.black, thick=5, psym=10

  ;; Close Ps
  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]


  print, 'fig_sdss:  All done!'
       
  return
end
      
      
