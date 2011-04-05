;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro fig_wfc3_stack

  compile_opt strictarr

  ;; Get structure if necessary
  if not keyword_set( PSFILE ) then psfile = 'fig_wfc3_stack.ps'
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
  ymrg = [7.0,7]
  xrng1=[600., 1400.]

  thk = 5
  yrng=[0.0, 1.7]

  ;; Start plot
  plot, [0], [0], color=clr.black, background=clr.white, $
        charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, xtitle='Wavelength (Angstroms)', $
        ytitle='Relative Flux', yrange=yrng, thick=4, $
        xrange=xrng1, ystyle=1, xstyle=1, psym=1, /nodata ;, /ylog

  ;; Read in stack
  fx = x_readspec(getenv('LLSPAP')+'/HST/WFC3/Analysis/wfc3_stack_qsos.fits', $ 
           wav=wv, infl=2)

  oplot, wv, fx, color=clr.black, thick=5, psym=10

  ;; Close Ps
  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]


  print, 'fig_sdss:  All done!'
       
  return
end
      
      
