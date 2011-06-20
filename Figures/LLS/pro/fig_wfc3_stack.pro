;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro fig_wfc3_stack, PSFILE=psfile, TILT=tilt

  compile_opt strictarr

  ;; Get structure if necessary
  if not keyword_set( PSFILE ) then psfile = 'fig_wfc3_stack.ps'
  if not keyword_set( CSZ ) then csz = 2.4
  if not keyword_set( LSZ ) then lsz = 1.9

  if not keyword_set( MAXCHI ) then MAXCHI = 1.3
  if not keyword_set( WVMIN ) then wvmin = 2000. 
  if not keyword_set( MAX_DZ ) then max_dz = 0.6
  if not keyword_set(OMEGA_M) then omega_m = 0.3
  if not keyword_set(H0) then H0 = 72.

  ;; Stack
  spec_fil = getenv('LLSPAP')+'/HST/WFC3/Analysis/wfc3_stack_qsos.fits'
  flux = xmrdfits(spec_fil,0,head,/sile)
  sig = xmrdfits(spec_fil,1,/sile)
  wave_rest = xmrdfits(spec_fil,2,/sile)
  npix = n_elements(flux)

  ;; Boot
  boot_fil = getenv('LLSPAP')+'/HST/WFC3/Analysis/wfc3_boot_stack.fits'
  boot_spec = xmrdfits(boot_fil,1)
  sigma_boot = fltarr(npix)
  avg_boot = fltarr(npix)
  ;; Stats
  for ii=0L,npix-1 do begin
     djs_iterstat, boot_spec[ii,*], sigma=sigm, mean=avg
     avg_boot[ii] = avg
     sigma_boot[ii] = sigm
  endfor

  ;; Read in Telfer
  readcol, getenv('LLSPAP')+'/taueff/Analysis/Telfer/hst_comp01_rq.asc', t_wv, t_fx

  ;; normalize
  t_1460 = where(abs(t_wv-1460.) LT 10)
  t_nrm = median(t_fx[t_1460])
  t_fx = t_fx / t_nrm

  ;; Use first QSO as a template 
  qso = xmrdfits(getenv('LLSPAP')+'/HST/WFC3/Analysis/wfc3_specstrct.fits',1)
  wv_wfc3 = qso[0].wave
  zem = qso[0].zem
  t_wv_obs = t_wv*(1+zem)
  ;; Rebin onto WFC3 grid
  x_specrebin, t_wv_obs, t_fx, wv_wfc3, t_fx_wfc3, /flambda
;  x_splot, t_wv_obs, t_fx, xtwo=wv_wfc3, ytwo=t_fx_wfc3, /bloc
  ;; Smooth by FWHM=5 pixels
  nsmooth = (5./2.3548)
  kernel = gauss_kernel(nsmooth)
  t_smooth_fx = convol(t_fx_wfc3, kernel)
;  x_splot, t_wv_obs, t_fx, xtwo=wv_wfc3, ytwo=t_smooth_fx, /bloc

  ;; Shift back
  t_wfc3_rest = wv_wfc3/(1+zem)
  ;; Bin onto stack grid
  x_specrebin, t_wfc3_rest, t_smooth_fx, wave_rest, t_fx_stack, /flam
;  x_splot, t_wfc3_rest, t_smooth_fx, xtwo=wave_rest, ytwo=t_fx_stack, /block, psym2=10

  ;; Normalize again
  t_1460 = where(abs(wave_rest-1460.) LT 10)
  t_nrm2 = median(t_fx_stack[t_1460])
  t_fx_stack = t_fx_stack / t_nrm2
  
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; PLOT
  x_psopen, psfile, /maxs
  clr = getcolor(/load)
  !p.multi=[0,1,1]
  lclr = clr.white

  ;; FIRST PLOT ::  Data and best model
  xmrg = [8,2]
  ymrg = [4.0,0.5]
;  xrng=[3900., 5020.]
  xrng=[600., 1500.]
  yrng=[0.0, 2.0]
  plot, [0], [0], color=lclr, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, xtitle='!9l!X!drest!N (Angstroms)', $
        ytitle='f!d!9l!X!N/f!d1450!N', yrange=yrng, thick=9, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata
;  oplot, [-1e9,1e9], [0.0, 0], color=clr.gray, lines=1
  
  ;; Plot
  x_curvefill, wave_rest, avg_boot+sigma_boot, avg_boot-sigma_boot, $
               color=clr.darkgray
;  oplot, wave_rest, avg_boot, color=clr.gray

  ;;
  oplot, wave_rest, flux, color=lclr, psym=10, thick=9

  ;;;;;
  ;; Plot Telfer
  gd = where(t_wv GT xrng[0] and t_wv LT xrng[1]) 
;  oplot, t_wv[gd], smooth(t_fx[gd],5), color=clr.red, thick=3
  oplot, wave_rest, t_fx_stack, color=clr.tomato;, linesty=1
  ;; Tilt
  if keyword_set(TILT) then begin
     piv_wv = 1450.         ; Ang
     alpha = -0.7
     oplot, wave_rest, t_fx_stack*(wave_rest/piv_wv)^alpha, color=clr.tomato, linesty=2
     alpha = 0.1
     oplot, wave_rest, t_fx_stack*(wave_rest/piv_wv)^alpha, color=clr.tomato, linesty=2
  endif

  ;; Close Ps
  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]


  print, 'fig_boot_stack:  All done!'
       
  return
end
      
      
