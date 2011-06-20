; .com ../Analysis/pro/dla_proxdx
; .com ../Analysis/pro/dla_proxstat
; fig_lox 

function fnhi_fig_lls_powerfn, logNHI
;; calculate d^Num / dN dz

common cmmn_powerfn, fn_strct, dxdz, teff_cross, beta_plls

  fNHI = dblarr(n_elements(logNHI))

  ;; Lya
  lya = where(logNHI LE 14.5, nlya)
  if nlya GT 0 then $
    fNHI[Lya] = 10.^(fn_strct.k_Lya + fn_strct.a_Lya*logNHI[lya]) 

  ;; pLLS
  plls_1 = where(logNHI GT 14.5 and logNHI LT fn_strct.N_pLLS, nplls1)
  if nplls1 GT 0 then $
    fNHI[plls_1] = 10.^(fn_strct.k_plls + fn_strct.a_pLLS*logNHI[plls_1])

  ;; LLS (and beyond)
  lls = where(logNHI GE fn_strct.N_pLLS and logNHI LT 19, nnlls)
  if nnlls GT 0 then $
    fNHI[lls] = 10.^(fn_strct.k_lls + logNHI[lls]*fn_strct.a_LLS)
  
  ;; SLLS
  slls = where(logNHI GE 19. and logNHI LT 20.3, nnslls) 
  if nnslls GT 0 then $
    fNHI[slls] = 10.^fn_strct.k_slls * 10.^(logNHI[slls]*fn_strct.a_SLLS)

  ;; DLA
  lodla = where(logNHI GE 20.3 and logNHI LT fn_strct.dla_str.Nd, nndla)
  if nndla GT 0 then fNHI[lodla] = fn_strct.dla_str.k3 * $
                10.^((logNHI[lodla] - fn_strct.dla_str.Nd)*fn_strct.dla_str.a3)
  hidla = where(logNHI GE fn_strct.dla_str.Nd, nndla)
  if nndla GT 0 then fNHI[hidla] = fn_strct.dla_str.k3 * $
    10.^((logNHI[hidla] - fn_strct.dla_str.Nd)*fn_strct.dla_str.a4)
  
  return, fNHI * dxdz
end

;;;;;;;  
function int_fig_lls_powerfn, NHI

  common cmmn_powerfn
  logNHI = alog10(NHI)
  
  hi = where(logNHI GT 19., nnhi, complement=low, ncomplem=nnlow)
  expt = dblarr(n_elements(NHI))
  if nnhi GT 0 then expt[hi] = 1.
  if nnlow GT 0 then expt[low] = 1.d - exp(-NHI[low]*teff_cross)
  
  return, fnhi_fig_lls_powerfn(logNHI) * expt 

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro fig_lls_powerfn, PSFILE=psfile, GZFIL=gzfil, DLALST=dlalst, $
                     RSDSS=rsdss, RTOT=rtot, TIMBIN=timbin, STRCT=strct, $
                     BINS=bins, OUTFIL=outfil, VPROX=vprox, $
                     VMAX=vmax, PLLS=plls, PROX=prox, QSOFIL=qsofil, $
                     LLSFIL=llsfil, DR5_FNSTR=dr5_fnstr, DLA_STR=dla_str, $
                     SLLS_STRCT=slls_strct, $
                     NOSLLS=noslls, PLT_EXP=plt_EXP, BOOST=boost, NOPS=nops, $
                     RESULTS=results, THK=thk, NOMFP=nomfp

  common cmmn_powerfn
  compile_opt strictarr

  ;; Get structure if necessary
  if not keyword_set( PSFILE ) then psfile = 'fig_lls_powerfn.ps'
  if keyword_set( NOMFP ) then psfile = 'fig_lls_powerfn_nomfp.ps'
  cd, getenv('LLSPAP')+'/SDSS/DR7/Analysis/pro/', curr=curr
  RESOLVE_ROUTINE, 'max_lozpower'
  lls_dr7_initparm, init
  cd, curr

  cd, getenv('LLSPAP')+'/SLLS/PaperI/Figures/pro/', curr=curr
  RESOLVE_ROUTINE, 'fig_fnfit'
  cd, '../../Analysis/pro'
  RESOLVE_ROUTINE, 'slls_stat', /is_function
  RESOLVE_ROUTINE, 'slls_dx', /is_function
  cd, curr

  if not keyword_set( LLSFIL ) then llsfil = init.llsfil
  if not keyword_set( QSOFIL ) then qsofil = init.qsofil 
  if not keyword_set( VPROX ) then vprox = init.vprox
  if not keyword_set( XZFIL ) then xzfil = init.xzfil
  if not keyword_set( BINS ) then BINS = init.bins
  if not keyword_set( MAXDZ ) then maxdz = init.maxoff 
  if not keyword_set( ZEM_MIN ) then zem_min = init.zem_min
  if not keyword_set( CSZ ) then csz = 2.1
  if not keyword_set( LSZ ) then lsz = 2.3
  if not keyword_set( STP ) then stp = 0.2
  if not keyword_set( SPL_SYM ) then spl_sym = 2
  if not keyword_set(NINT) then nint = 200L 

  if not keyword_set(OBS_MFP) then obs_mfp = 45.
  if not keyword_set(SIG_MFP) then sig_mfp = 7.
  if not keyword_set(H0) then H0 = 72.
  if not keyword_set(OMEGA_M) then omega_m = 0.3

  if not keyword_set( PLT_EXP ) then PLT_EXP = 0L
  if not keyword_set(BOOST) then boost = 1.

  if not keyword_set(THK) then thk = 5

  if BOOST GT 1. then stop

  ;; Output out file
  close, /all

  ;; DLA
  if not keyword_set( GZFIL ) then $
    gzfil = '~/SDSS/DR5_QSO/dr5_dlagz_s2n4.fits'
  bins = [ [2.4,2.6], [3.4,4.0]]
  if not keyword_set(DLA_STR) then $
  sdss_dlalox, GZFIL=gzfil, STRCT=dla_str, BINS=bins
  if not keyword_set(DR5_FNSTR) then $
    sdss_fndla, GZFIL=GZFIL, STRCT=dr5_fnstr,  STP=STP, XZFIL=xzfil, $
                /SALL, ALLDR=alldr, CL=0.95, NPLT=nplt, ZBIN=[3.4, 4.0]
  lX_DLA = dla_str.dndx[1]

  ;; MFP
  c = x_constants()
  teff_engy = c.Ryd/c.eV
  teff_cross = x_photocross(1,1,teff_engy)

  zqso = 3.7
  drdz = 3e5/H0/sqrt(omega_m) / (1+zqso)^(2.5) ;; Mpc

  ;; f(N) structure
  dXdz = abs(cosm_xz(3.75)-cosm_xz(3.65)) / 0.1
  fn_strct = { $
             a_Lya: -1.5, $
             k_Lya: 9.12 - alog10(dXdz), $
             a_pLLS: 0., $
             k_pLLS: 0., $
             N_pLLS: 0., $
             a_LLS: 0., $
             k_LLS: 0., $
             a_SLLS: 0., $
             k_SLLS: 0., $
             dla_str: dr5_fnstr $
          }
  

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; PLOT
  if keyword_set( PSFILE ) and not keyword_set(NOPS) then $
    x_psopen, psfile, /maxs
  clr = getcolor(/load)
  lclr = clr.white
  !p.multi=[0,1,1]
  
  xmrg = [8,2]
  ymrg = [4.0,1]
  case PLT_EXP of 
      0: begin
          yrng=[-26, -9.]
          ytit='log [f(N!dHI!N,X)]'
      end
      1: begin
          yrng=[-5, 3]
          ytit='log [N!dHI!N f(N!dHI!N,X)]'
      end
      2: begin
          yrng=[-5, 3]
          ytit='log [N!dHI!N!u2!N f(N!dHI!N,X)]'
      end
      else: stop
  endcase
  if not keyword_set(XRNG) then xrng=[12.0, 22.5]
  plot, [0], [0], color=lclr, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, xtitle='log N!dHI!N', $
        ytitle=ytit, yrange=yrng, thick=9, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata, $
        xthi=thk, ythi=thk

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Lya forest 
  xpt = [12., 14.5]
  npt = n_elements(xpt)


  log_fLyaz = 9.12-1.5*xpt ;; f(N,z)
  log_fLyaX = alog10(10.d^log_fLyaz / dXdz)

  ypt = [log_fLyaX]
  oplot, xpt, log_fLyaX+xpt*(PLT_EXP), color=clr.lightgray, thick=thk

  xyouts, 13.2, -10.9, 'Ly!9a!X', color=clr.lightgray, charsiz=lsz

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;
  ;; DLA
  xval = dr5_fnstr.xval
  yplt = dr5_fnstr.yplt
  bins = dr5_fnstr.bins
  yerr1 = dr5_fnstr.yerr1
  yerr2 = dr5_fnstr.yerr2

  val = where(yplt GT -99., complement=lim, ncomplement=nlim)
  oploterror, [xval[val]], [yplt[val]], [xval[val]-bins[val]],  $
              [yplt[val]-yerr2[val]], psym=3, thick=6, $
              color=clr.tomato, /lobar, errcolor=clr.tomato, errthi=6
  oploterror, [xval[val]], [yplt[val]], [bins[val]+stp-xval[val]],  $
              [yerr1[val]-yplt[val]], psym=3, thick=6, $
              color=clr.tomato, /hibar, errcolor=clr.tomato, errthi=6
  if nlim NE 0 then begin
      plotsym,1, 2., thick=3
      oploterror, [xval[lim]], [yerr1[lim]], [xval[lim]-bins[lim]], $
                  replicate(0.,nlim), $
                  color=clr.tomato, psym=8, /lobar, errcolor=clr.tomato, thick=6, errthi=6
      oploterror, [xval[lim]], [yerr1[lim]], [bins[lim]+stp-xval[lim]], $
                  replicate(0., nlim), $
                  color=clr.tomato, psym=8, /hibar, errcolor=clr.tomato, thick=6, errthi=6
  endif


  xyouts, 20.85, -21.8, 'DLA', color=clr.tomato, charsiz=lsz

  fplt = 20.3 + (22-20.3)*findgen(1000L)/1000.
  thy = fltarr(n_elements(fplt))
  ns10 = 10.^dr5_fnstr.Nd
  lw = where(fplt LT dr5_fnstr.Nd, complement=hi, ncomplement=numhi)
  thy[lw] = dr5_fnstr.k3 * (10^fplt[lw]/ns10)^dr5_fnstr.a3
  if numhi NE 0 then thy[hi] = dr5_fnstr.k3 * $
    (10^fplt[hi]/ns10)^dr5_fnstr.a4
  oplot, fplt, alog10(thy) + fplt*PLT_EXP, color=lclr, thick=thk

  lX_DLA = dla_str.dndx[1]
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; SLLS
  if not keyword_set(lX_SLLS) then lX_SLLS = init.lX_SLLS 
  a_SLLS = -1.2    ;; Full sample
  N1 = 1d19
  N2 = 2d20
  k_SLLS = alog10(lX_SLLS * (1+a_SLLS) / (N2^(a_SLLS+1) - N1^(1+a_SLLS)))
  fplt = 19.0 + 1.3*findgen(1000)/1000.
  thy0 = 10.^k_SLLS * (10.^fplt)^a_SLLS

;  lX_SLLS = 10.^k_SLLS * (N2^(a_SLLS+1) - N1^(1+a_SLLS)) / $
;            (1+a_SLLS)
;  print, 'SLLS lX: ', lX_SLLS
  a1_SLLS = -1.4
  k1_SLLS = alog10(lX_SLLS * (1+a1_SLLS) / (N2^(a1_SLLS+1) - N1^(1+a1_SLLS)))
  thy1 = 10.^k1_SLLS * (10.^fplt)^a1_SLLS
;  oplot, fplt, alog10(thy1), color=clr.green, linesty=1
  a2_SLLS = -1.01
  k2_SLLS = alog10(lX_SLLS * (1+a2_SLLS) / (N2^(a2_SLLS+1) - N1^(1+a2_SLLS)))
  thy2 = 10.^k2_SLLS * (10.^fplt)^a2_SLLS
;  oplot, fplt, alog10(thy2), color=clr.green, linesty=2
  x_curvefill, [fplt[0],fplt[999]], $
               alog10([thy1[0],thy2[999]]) + [fplt[0],fplt[999]]*PLT_EXP, $
               alog10([thy2[0],thy1[999]]) + [fplt[0],fplt[999]]*PLT_EXP, $
               color=clr.yellow
  oplot, fplt, alog10(thy0) + fplt*PLT_EXP, color=lclr, thick=thk

  xyouts, 19.5, -20.3, 'SLLS', color=clr.yellow, charsiz=lsz

  logfn_19 = alog10(10.^k_SLLS * (10.^19.)^a_SLLS)

  fn_strct.a_SLLS = a_SLLS
  fn_strct.k_SLLS = k_SLLS


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; LLS
  lX_LLS = 0.5 * BOOST   ;; Good for z=3.7; Error is 0.1
;  avec = -0.3 - 0.65*findgen(1000)/1000.
  Nt2 = 10.d^(17.5)
  N19 = 1d19
  avec = -0.10005 - 1.5*dindgen(100)/100.

  alls = fltarr(3,3)
  alls = fltarr(3,3)
  klls = fltarr(3,3)
  klls = fltarr(3,3)

  ;; MFP stuff
  nmfp = 15
  NHI_MFP = 15. + (2.5)*findgen(nmfp)/(nmfp-1)
  Nmin = 1d12
  Nmax = 1d22
  dlnN = alog(Nmax/Nmin) / (nint-1)
  Nval = Nmin * exp( dlnN * dindgen(nint) )

  sv_LXS = fltarr(3,3)
  sv_LXt = fltarr(3,3)
  sv_f19 = fltarr(3,3)
  ;;;;;;;
  ;; LLS constraint only
  for ii=1,2 do begin
      ;; Set SLLS and kvec
      if ii EQ 1 then f_SLLS = logfn_19-0.2 else f_SLLS = logfn_19+0.2
      kvec = f_SLLS - avec*19
      ;; Vary LX_LLS
      for jj=1,2 do begin
          if jj EQ 1 then LX_tLLS = 0.15 else LX_tLLS = 0.35 ;; 17.5-19 only
          ;; BOOST
          LX_tLLS = LX_tLLS + (lX_LLS - 0.5)
          ;;
          lvec = 10.^kvec / (1+avec) * ( N19^(1+avec) - Nt2^(1+avec) )
          mn = min(abs(lvec-LX_tLLS), imn)

          ;; Save 
          alls[ii,jj] = avec[imn]
          klls[ii,jj] = kvec[imn]
          sv_lXt[ii,jj] = lX_tLLS
          sv_f19[ii,jj] = f_SLLS
      endfor
  endfor

  ;; Central
  lX_LLS = 0.5 * BOOST   ;; Good for z=3.7; Error is 0.1
  kvec = k_SLLS + 19*(a_SLLS-avec)
  lvec = 10.^kvec / (1+avec) * ( N19^(1+avec) - Nt2^(1+avec) )
  lX_tLLS = LX_LLS - lX_SLLS - lX_DLA
  mn = min(abs(lvec-LX_tLLS), imn)
  ;; Best
  a_LLS = avec[imn]
  k_LLS = kvec[imn]

  ;; Save
  alls[0,0] = a_LLS
  klls[0,0] = k_LLS
  sv_lXt[0,0] = lX_tLLS
  sv_f19[0,0] = alog10(10.^k_SLLS * (10.^19.)^a_SLLS)

  ;; RESULTS
  results = { $
            a_LLS: aLLS, $
            k_LLS: kLLS, $
            lXt: sv_lXt, $
            f19: sv_f19, $
            bplls: fltarr(3,3,nmfp), $
            Nplls: fltarr(3,3,nmfp), $
            Cplls: fltarr(3,3,nmfp), $
            msk_plls: bytarr(3,3,nmfp) $
            }


  ;;;;;;;;;
  ;; Plot
  npt = 101L
  fplt = 17.5 + 1.5*findgen(npt)/(npt-1)
  up_lls = replicate(-9999, npt)
  low_lls = replicate(1., npt)
  for ii=1,2 do begin
      for jj=1,2 do begin
          if aLLS[ii,jj] LT (-0.7) then begin
              up_lls = up_lls > (kLLS[ii,jj] + fplt*aLLS[ii,jj]) 
              low_lls = low_lls < (kLLS[ii,jj] + fplt*aLLS[ii,jj]) 
          endif
      endfor
  endfor
  x_curvefill, fplt, up_lls+fplt*PLT_EXP, $
               low_lls+fplt*PLT_EXP, color=clr.cyan

  oplot, fplt, kLLS[0,0] + fplt*aLLS[0,0] + fplt*PLT_EXP, $ 
         color=lclr, thick=thk

  xyouts, 18.1, -18.5, 'LLS', color=clr.cyan, charsiz=lsz

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; MFP
  fplt = 14.5 + 3.0*findgen(npt)/(npt-1)
  all_gdbeta = [0]
  for ii=0,2 do begin
      for jj=0,2 do begin
          if jj GT 0 and ii EQ 0 then continue

          fn_strct.a_LLS = aLLS[ii,jj]
          fn_strct.k_LLS = kLLS[ii,jj]

          sv_mfp = fltarr(nmfp)
          sv_beta = fltarr(nmfp)
          sv_kpLLS = fltarr(nmfp)
          for tt=0L,nmfp-1 do begin
              ;; Set break and alpha
              fn_strct.N_pLLS = NHI_MFP[tt]
              f_break = kLLS[ii,jj] + NHI_MFP[tt]*aLLS[ii,jj]
              beta = (f_break-log_fLyaX[1]) / (NHI_MFP[tt]-14.5)
              fn_strct.a_pLLS = beta
              fn_strct.k_pLLS = log_fLyaX[1] - beta*14.5

              ;; Get the MFP
              intgrnd = int_fig_lls_powerfn(Nval) * Nval 
              dteffdz = total(intgrnd*dlnN)
              sv_mfp[tt] = drdz / dteffdz  ;; Mpc
              sv_beta[tt] = beta
              sv_kpLLS[tt] = log_fLyaX[1] - beta*14.5
          endfor
          gd_mfp = where(abs(sv_mfp-OBS_MFP) LT SIG_MFP, ngd_mfp)
          if ngd_mfp GT 0 then begin
              all_gdbeta = [all_gdbeta, sv_beta[gd_mfp]]
              mn = min(abs(sv_mfp-OBS_MFP),imn)
              print, 'Best: ', aLLS[ii,jj], kLLS[ii,jj], sv_beta[imn]
              ;; Save
              results.msk_plls[ii,jj,gd_mfp] = 1B
              results.bplls[ii,jj,gd_mfp] = sv_beta[gd_mfp]
              results.Nplls[ii,jj,gd_mfp] = NHI_MFP[gd_mfp]
              results.CpLLS[ii,jj,gd_mfp] = log_fLyaX[1] - sv_beta[gd_mfp]*14.5
          endif
;          printcol, NHI_MFP, sv_beta, sv_mfp
;          stop
          ;; Plot
          if ii GT 0 and jj GT 0 then begin
              for kk=0L,ngd_mfp-1 do begin
                  tt = gd_mfp[kk]
                  fn_strct.N_pLLS = NHI_MFP[tt]
                  fn_strct.a_pLLS = sv_beta[tt]
                  fn_strct.k_pLLS = log_fLyaX[1] - sv_beta[tt]*14.5
                  fval = alog10(fnhi_fig_lls_powerfn(fplt) / dxdz)
                  ;;
                  if not keyword_set(NOMFP) then $
                     oplot, fplt, fval + fplt*PLT_EXP, color=clr.orange, $
                            thick=2, linesty=(ii-1)
              endfor
          endif 
          if ii EQ 0 and jj EQ 0 then begin
              tt = median(gd_mfp)
              fn_strct.N_pLLS = NHI_MFP[tt]
              fn_strct.a_pLLS = sv_beta[tt]
              fn_strct.k_pLLS = log_fLyaX[1] - sv_beta[tt]*14.5
              fval = alog10(fnhi_fig_lls_powerfn(fplt) / dxdz)
              if not keyword_set(NOMFP) then $
                 oplot, fplt, fval + fplt*PLT_EXP, color=lclr, thick=thk
          endif
      endfor
  endfor
;  all_gdbeta = all_gdbeta[1:*]
  if not keyword_set(NOMFP) then $
     print, 'Good beta values: ', all_gdbeta


  if not keyword_set(NOMFP) then $
     xyouts, 16.0, -14.9, '!9l!X!dmfp!N', color=clr.orange, charsiz=lsz

;  printcol, alls, klls

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Labels

  ;;; beta
  xyouts, 20.8, -24.5,  '!9b!X=!9-!X1.8', color=lclr, charsiz=lsz, align=0.5
  xyouts, 19.5, -22.,  '!9b!X=!9-!X1.2', color=lclr, charsiz=lsz, align=0.5
  xyouts, 18.0, -21.,  '!9b!X=!9-!X0.8', color=lclr, charsiz=lsz, align=0.5
  if not keyword_set(NOMFP) then $
     xyouts, 15.0, -18.,  '!9b!X=(!9-1.9,!9-!X5)', color=lclr, $
             charsiz=lsz, align=0.5
  xyouts, 13.5, -13.3,  '!9b!X=!9-!X1.5', color=lclr, charsiz=lsz, align=0.5

  ;; redshift
  xyouts, 20.5, -11.,  'z~3.7', color=lclr, charsiz=lsz*1.5, align=0.0

  ;; Close Ps
  if keyword_set( PSFILE ) then x_psclose
  print, 'fig_lls_powerfn:  All done!'
  close, /all

  strct = fn_strct

  return
end
      
      
