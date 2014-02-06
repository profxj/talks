;+ 
; NAME:
; fuse_velplt
;  V1.1
;
;   fig_ew_vs_rho, 1238.821
;   fig_ew_vs_rho, 1334.5323
;   fig_ew_vs_rho, 1548.195
;   fig_ew_vs_rho, 2796.352
;; /NO_STAT -- Don't perform ASURV stats
;; /CLIP_HIGH -- Show <EW> restricted to data with <0.5Ang
;------------------------------------------------------------------------------
pro qpq7_ew_vs_rho, wrest, FOREST=forest, SIGEWlim=SIGEWlim, $
                   NOPS=nops, CSIZE=csize, LBL=lbl, NO_STATS=no_stats, $
                   CLIP_HIGH=clip_high

  compile_opt strictarr

  cd, getenv('QSO_DIR')+'tex/QPQ7/Analysis/AbsLines/pro/', curr=curr
  resolve_routine, 'qpq7_mkdata', /compile_full_file, /is_function
  resolve_routine, 'qpq7_sigewlim', /compile_full_file, /is_function
  cd, curr
  cd, getenv('QSO_DIR')+'tex/QPQ7/Analysis/pro/', curr=curr
  resolve_routine, 'run_asurv', /compile_full_file
  cd, curr

  if not keyword_set( LSZ ) then lsz = 1.8
  if not keyword_set( LTHICK ) then lthick = 4.
  if not keyword_set( PTHICK ) then pthick = 2.
  if not keyword_set( BLSZ ) then blsz = 1.7
  if not keyword_set( CSIZE ) then csize = 2.5
  if not keyword_set( SSIZE ) then ssize = 0.9

  if not keyword_set( CONTI ) then conti = 0.8

  ;; Grab transition info
  if not keyword_set(WREST) then wrest = 1548.195d ;
  getion, wrest, ioni, elm, Z=ionz, NM=strioni

  ;; sigEWlim
  if not keyword_set( SIGEWlim ) then $
     sigEWlim = qpq7_sigewlim(wrest)

  ;; QPQ7
  qpq_fil = '~/Dropbox/QSOPairs/qpq7_pairs.fits'
  qpq_strct = xmrdfits(qpq_fil, 1)

  if not keyword_set(FOREST) then begin ;; Cut out Lya forest
     good = where( (qpq_strct.z_fg+1)*wrest GT ((qpq_strct.z_bg+1)*1215.6701 + 20.), ngd)
     qpq_strct = qpq_strct[good]
  endif

  fwr = fix(wrest)
  cfwr = strtrim(fix(wrest),2)
  if not keyword_set(PSFILE) then psfile = 'fig_ew_'+cfwr+'_vs_rho.ps'

  ;; Read in LBG
  if keyword_set(LBG) then begin
     lls_struct, lbg, '/u/xavier/LLS/Lists/lbg_lls.lst', ROOT='/u/xavier/LLS/', /ew
     ion = where(abs(lbg.systems[0].ion[ionZ].state[ioni,*].lambda-wrest) LT 1e-3, nlbg)
     lbg_idx = ion / 30
     nobj = nlbg  + 2 ;; Stacked spectra
  endif
  
  bins = [ [30., 100], [100., 200], [200, 500], [500., 1000]] ; kpc
  sz_bins = size(bins, /dim)
  avgR = fltarr(sz_bins[1])
  avgEW = fltarr(sz_bins[1])
  avgsigEW = fltarr(sz_bins[1])
  low_avgEW = fltarr(sz_bins[1])
  low_avgsigEW = fltarr(sz_bins[1])

  ;; Plot
  if keyword_set( PSFILE ) and not keyword_set(NOPS) then begin
     x_psopen, psfile, /maxs
     !p.multi=[0,1,1]
  endif 
  if keyword_set(NOPS) then begin
     xmrg = [7, 2]
     ymrg = [17, 1]
  endif else begin
     xmrg = [8, 14]
     ymrg = [4, 1]
  endelse
  clr = getcolor(/load)
  lbgc = clr.red
  fclr = clr.lightgray
  pclr = clr.white
  aclr = clr.cyan


  ;; ;;;;;;;;;;;;;;;;;;
  ;; Histograms
  yrng = [0.01, 5.] ; Ang
  xrng = [30., 1000] ; Number
  thisletter = byte(94)
  perpletter = '!9' + string(thisletter) + '!X'
  xtitle='R!d'+perpletter+'!N (kpc)'
  plot, [0], [0], xrange=xrng, $
        yrange=yrng, /NODATA, xmarg=xmrg, ymarg=ymrg, $
        charsize=csize, psym=10, background=clr.white, color=fclr, $
        xstyle=1, ystyle=1, thick=pthick, ytickinterval=5.,$
        ytitle='W!d'+cfwr+'!N (Ang)', xtit=xtitle, /ylog, /xlog

  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; QPQ7
  gd = where( abs(qpq_strct.metal_wrest-wrest) LT 1e-3 AND $
              qpq_strct.metal_sigEW LT sigEWlim AND $
              (qpq_strct.flg_metal_eye MOD 256) GE 128 AND $
              qpq_strct.flg_metal_EW EQ 1, ngd)
  print, 'Number of measurements = ', ngd
  sz = size(qpq_strct.metal_wrest,/dimen)
  gdi = gd / sz[0]
  QPQ7_EW = (qpq_strct.metal_EW)[gd] > 0.013
  if ngd gt 0 then begin
     plotsym, 8, ssize, /fill
     oplot, qpq_strct[gdi].R_phys, QPQ7_EW, color=pclr, psym=8
  endif
  sv_gd_EW = QPQ7_EW
  sv_gd_R = qpq_strct[gdi].R_phys

  low = where(QPQ7_EW LT 0.05,nlow)
  ;if nlow GT 0 then stop

  ;; ;;;;;;;
  ;; Limits
  ulim = where( abs(qpq_strct.metal_wrest-wrest) LT 1e-3 AND $
                qpq_strct.metal_sigEW LT sigEWlim AND $
                (qpq_strct.flg_metal_eye MOD 256) GE 128 AND $
                qpq_strct.flg_metal_EW EQ 3, nulim)
  print, 'Number of upper limits = ', nulim
  if nulim GT 0 then begin
     ui = ulim / sz[0]
     QPQ7_EW = (qpq_strct.metal_EW)[ulim] 
     QPQ7_sigEW = (qpq_strct.metal_sigEW)[ulim] 
     QPQ7_lim_EW = QPQ7_EW > (2.* QPQ7_sigEW)
                                ;stop ;; check that many of the EWs are negative (not as many as I expected..)
     plotsym, 8, ssize
     oplot, qpq_strct[ui].R_phys, QPQ7_lim_EW > 0.013, color=pclr, psym=8
     plotsym, 1, ssize
     oplot, qpq_strct[ui].R_phys, QPQ7_lim_EW > 0.013, color=pclr, psym=8
  endif

  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Asurv
  if not keyword_set(NO_STATS) then begin
     all_EW = [sv_gd_EW, QPQ7_lim_EW]
     all_R = [sv_gd_R, qpq_strct[ui].R_phys]
     flg_sens = [replicate(0,n_elements(sv_gd_EW)), replicate(-1, nulim)]
     run_asurv, all_EW, flg_sens, var2=all_R, /KT, prob=prob 
     print, 'Asurv probability = ', prob
  endif
  

  ;; 2 sigma??
  ;tsig = where( QPQ7_EW GT 2.*QPQ7_sigEW, ntsig) 
  ;if ntsig GT 0 then oplot, [qpq_strct[ui[tsig]].R_phys], [QPQ7_EW[tsig]], color=clr.gray, psym=2

  ;; Binned evaluations
  all = where( abs(qpq_strct.metal_wrest-wrest) LT 1e-3 AND $
               qpq_strct.metal_sigEW LT sigEWlim AND $
              (qpq_strct.flg_metal_eye MOD 256) GE 128 AND $
              qpq_strct.flg_metal_EW GT 0, nall)
  alli = all / sz[0]
  Rperp = qpq_strct[alli].R_phys
  EW = (qpq_strct.metal_EW)[all] 
  sigEW = (qpq_strct.metal_sigEW)[all] 
  for ii=0L,sz_bins[1]-1 do begin
     gdR = where( Rperp GT bins[0,ii] and Rperp LT bins[1,ii], ngdR)
     ;; Stats
     avgR[ii] = mean(Rperp[gdR])
     
     stats = moment(EW[gdR])
     avgEW[ii] = stats[0]
     avgsigEW[ii] = sqrt(stats[1]/float(ngdR))

     ;; Clip the large EWs
     lowerEW = where( EW[gdR] LT 0.5, nlower )
     stats2 = moment(EW[gdR[lowerEW]])
     low_avgEW[ii] = stats2[0]
     low_avgsigEW[ii] = sqrt(stats2[1]/float(nlower))
  endfor
  plotsym, 0, 1.5, /fill
  oploterror, avgR, avgEW, avgsigEW, color=aclr, errcol=aclr, psym=8
  if keyword_set(CLIP_HIGH) then $
     oploterror, avgR, low_avgEW, low_avgsigEW, color=clr.red, errcol=clr.red, psym=8

  ;; LBL ;;;;;;;;;;;;;;;;;;
  if keyword_set(LBL) then xyouts, 35., 0.015, LBL, color=fclr, charsi=1.5

  if keyword_set( PSFILE ) and not keyword_set(NOPS) then begin
     x_psclose
     !p.multi=[0,1,1]
  endif
  
  return
end
