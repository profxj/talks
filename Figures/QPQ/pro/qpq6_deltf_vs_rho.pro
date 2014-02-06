;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  Plots delta F vs R 
pro qpq6_deltf_vs_rho, STRCT=strct,  NOSTOP=nostop, XLOG=xlog

  compile_opt strictarr

  ;cd, '../Analysis/AbsLines/pro/', curr=curr
  ;resolve_routine, 'qpq6_mkdata', /compile_full_file, /is_function
  ;resolve_routine, 'qpq6_grab_idlines', /compile_full_file, /is_function
  ;cd, curr

  ;; Get structure if necessary
  if not keyword_set( PSFILE ) then psfile = 'qpq6_deltf_vs_rho.ps'
  if keyword_set(XLOG) then psfile = 'fig_deltf_vs_rho_xlog.ps'
  if not keyword_set( CSZ ) then csz = 1.9
  if not keyword_set( LSZ ) then lsz = 1.9

  if not keyword_set(SIG_FX) then sig_fx = 0.05 ;; Error in tau values
  if not keyword_set( MAXCHI ) then MAXCHI = 1.3
  if not keyword_set( NRNG ) then nrng = 50L

  qsodir = getenv('QSO_DIR')

  ;; Read in pair structure
  qpq_fil = '~/Dropbox/QSOPairs/qpq6_final_cut.fits'
  qpq_strct = xmrdfits(qpq_fil, 1)

  ;; Redshift cut
  gdz = where(qpq_strct.z_fg GT 2.0 and qpq_strct.z_fg LT 3.); and $
;              qpq_strct.R_phys LT 500.)
  gd_strct = qpq_strct[gdz]


  ;; S2N cut
  S2N = 8.
  if keyword_set(S2N) then begin
     gds2 = where(gd_strct.S2N_Lya GT S2N)
     gd_strct=  gd_strct[gds2]
  endif
  npair = n_elements(gd_strct)
  print, 'Npair = ', npair

  ;; Delta <F>
  VINT = 1000.  ; km/s
  ewi = round( median( where(abs(gd_strct.vEW - VINT) LT 1.)  $
                       MOD n_elements(gd_strct[0].vEW) ) )
  ;EWLya = gd_strct.vEWLya[ewi]
  deltF = gd_strct.deltFLya[ewi]
  avgF = gd_strct.avgFLya[ewi]

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; PLOT
  x_psopen, psfile, /maxs
  clr = getcolor(/load)
  !p.multi=[0,1,1]
  fclr = clr.white
  pclr = clr.cyan
  
  xmrg = [8,15]
  ymrg = [4.0,0.5]

  xrng = [0., 1000]
  if keyword_set(XLOG) then xrng = [30., 1000]
  dLam = 1215.6701 * (2*VINT) / 3e5 ; Ang
  yrng = [-1.0,1.0]

  thisletter = byte(94)
  perpletter = '!9' + string(thisletter) + '!X'
  xtit='R!d'+perpletter+'!N (kpc)'
  ;ytit='!9d!X!d<F>'
  ytit='Excess HI Ly!9a!X Absorption'
  plot, [0], [0], color=fclr, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, xtitle=xtit, $
        ytitle=ytit, yrange=yrng, thick=5, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata, xlog=xlog
  
  oplot, xrng, [0., 0.], color=clr.gray, linesty=2

  plotsym, 0, 0.5, /fill
  oplot, gd_strct.R_phys, deltF, psym=8, color=fclr

  ;; Linear fit
  fit = linfit(gd_strct.R_phys, deltF, sigma=sig) 
  xplt = findgen(1001)
  yplt = fit[0] + xplt*fit[1]
  oplot, xplt, yplt, color=clr.cyan, linesty=2, thick=4
  printcol, fit, sig

  ;; Label 
  if keyword_set(XLOG) and xrng[0] GT 10 then $
     xyouts, xrng[0], yrng[0]-0.1, strtrim(round(xrng[0]),2), $
             color=clr.black, charsiz=csz, align=0.5

  ;; Correlate
  print, 'Spearman: ', r_correlate(gd_strct.R_phys, deltF)
  print, 'Kendall: ', r_correlate(gd_strct.R_phys, deltF, /kendall)
  
  ;; Close Ps
  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]

  print, 'fig_deltf_vs_rho:  All done!'
       
  return
end
      
      
