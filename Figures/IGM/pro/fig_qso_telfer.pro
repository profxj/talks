;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro fig_qso_telfer

  compile_opt strictarr

  ;; Get structure if necessary
  if not keyword_set( PSFILE ) then psfile = 'fig_qso_telfer.ps'
  if not keyword_set( CSZ ) then csz = 1.8
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
  yrng=[-0.02, 3.0]

  ;; Start plot
  plot, [0], [0], color=clr.black, background=clr.white, $
        charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, xtitle='!9l!X!drest!N (Angstroms)', $
        ytitle='f!d!9l!X!N', yrange=yrng, thick=4, $
        xrange=xrng1, ystyle=1, xstyle=9, psym=1, /nodata ;, /ylog
  oplot, [-1e9,1e9], [0.0, 0], color=clr.darkgray, lines=1
      

  ;; Read in Telfer

  ;; Close Ps
  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]


  print, 'fig_sdss:  All done!'
       
  return
end
      
      
