;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  Plots OThick absorbers
pro qpq6_corr_llsonly, COVFACT = COVFACT, psfile = psfile, outfile = outfile, NOPS=nops, $
                  CSZ=csz, LSZ=lsz, XLBL=xlbl, XMRG=xmrg, RVIR=rvir, LBLX=lblx, $
                  LBG=lbg, NO_RAN=no_ran, XLOG=xlog, $
                  XRNG=xrng, STRCT=strct, LIT_H=LIT_H

  compile_opt strictarr

  cd, getenv('QSO_DIR')+'tex/QPQ6/Analysis/Clustering/pro/', curr=curr
  resolve_routine, 'xcorr_line', /compile_full_file, /is_function
  cd, curr

  if not keyword_set(LIT_H) then LIT_H = 0.70  ;; (Hubble)

  ;; Get structure if necessary
  if not keyword_set( PSFILE ) then psfile = 'qpq6_corr_llsonly.ps'
  if not keyword_set( CSZ ) then csz = 3.2
  if not keyword_set( LSZ ) then lsz = 1.3

  if not keyword_set(SIG_FX) then sig_fx = 0.05 ;; Error in tau values
  if not keyword_set( MAXCHI ) then MAXCHI = 1.3
  ;  if not keyword_set( WVMIN ) then wvmin = 700. 
  if not keyword_set(OMEGA_M) then omega_m = 0.3

  if not keyword_set( NRNG ) then nrng = 50L
  if n_elements(xlog) EQ 0 then xlog = 1

  c = x_constants()

  ;; Read in pair stuff
  qpq_fil = '~/Dropbox/QSOPairs/qpq6_final_cut.fits'
  qpq6 = xmrdfits(qpq_fil, 1)
  gds = where(round(qpq6.s2n_lya) GE 10 and qpq6.z_fg LT 3.)

  qpq_strct = qpq6[gds]
  npair = n_elements(qpq_strct)

  ;; Correlation function setup
  meanz = mean(qpq_strct.z_fg)
  anow = 1.0D/(1.0D + meanz)
  Hz = cosm_hubble(meanz, /init, /w05map)
  vmax = 1500. ;; km/s

  ;; IGM
  fn_strct = xmrdfits(getenv('XIDL_DIR')+'/IGM/fN_empirical/fn_z2.4_5param.fits', 1,/silen)

  ;; R steps
  if not keyword_set(XRNG) then xrng = [50.,4000]/1e3  ; Mpc

  nstep = 5
  Rmin = 80.
  Rmax = 4000.
  Rstep = (Rmax/Rmin)^(1./nstep)


  ;; Write
  ;writecol, 'qpq5_smpl.dat', qpq_strct.z_fg, qpq_strct.r_phys, qpq_strct.flg_OThick

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; PLOT
  thisletter = byte(94)
  perpletter = '!9' + string(thisletter) + '!X'
  x_psopen, psfile, /portrait
  !p.multi=[0,1,4]
  xmrg = [7,14]
  ymrg = [0.0,0.0]
  clr = getcolor(/load)
  ;pclr = clr.white
  pclr = clr.gray
  lclr = pclr
  eclr = lclr

  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Correlation functions
  correl= xmrdfits(getenv('QSO_DIR')+'tex/QPQ6/Analysis/Clustering/qpq6_clustering.fits',1)
  gamma = 1.6

  for ss=2,2 do begin
     
     case ss of
        0: begin
           lbl = 'DLA'
           logN_max = 30.
           yrng = [1e-1, 2e2]
           ;; Get answers
           mt = where(abs(correl.gamma - gamma) LT 1e-3 and $
                      strmatch(strtrim(correl.lbl,2),lbl,/fold_case), nmt)
           if nmt NE 1 then stop
           r0 = correl[mt].r0
           sig_r0 = correl[mt].sig_r0
           logN_min = correl[mt].Nmin
        end
        1: begin
           lbl = 'SLLS'
           logN_max = 20.3
           yrng = [1e-1, 2e2]
           ;; Get answers
           mt = where(abs(correl.gamma - gamma) LT 1e-3 and $
                      strmatch(strtrim(correl.lbl,2),lbl,/fold_case), nmt)
           if nmt NE 1 then stop
           r0 = correl[mt].r0
           sig_r0 = correl[mt].sig_r0
           logN_min = correl[mt].Nmin
        end
        2: begin
           lbl = 'LLS'
           logN_max = 30.
           yrng = [1e-1, 2e2]
           ;; Get answers
           mt = where(abs(correl.gamma - gamma) LT 1e-3 and $
                      strmatch(strtrim(correl.lbl,2),lbl,/fold_case), nmt)
           if nmt NE 1 then stop
           r0 = correl[mt].r0
           sig_r0 = correl[mt].sig_r
           logN_min = correl[mt].Nmin
        end
     endcase

     if ss LT 2 then xspaces = replicate(' ',30) else delvarx, xspaces
     ;; y label
     if ss EQ 2 then ytit = '!9c!X!d'+perpletter+'!N (R,!9D!Xv)' else ytit = ''
     if ss EQ 2 then xtit='R!S!d'+perpletter+'!N!R!ucom!N (h!u-1!N Mpc)' else xtit = ''

     plot, [0], [0], color=lclr, background=clr.white, charsize=csz,$
           xmargin=xmrg, ymargin=ymrg, xtitle=xtit, $
           ytitle=ytit, yrange=yrng, thick=4, xtickn=xspaces, $
           xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata, $
           /XLOG, /ylog, ytickformat='x_logticks'

     ;; Generate the correlation function
     vdimless=vmax/Hz/anow/R0   ;

     ;; Correlate
     nplt = 30L
     Rval = (xrng[0] + (xrng[1]-xrng[0])*findgen(nplt) / (nplt-1))*1e3 ;; kpc

     ;; Error
     sig_model_chiP0 = fltarr(nplt)
     sig_model_chiP1 = fltarr(nplt)
     for jj=0L,nplt-1 do begin
        if (r0-sig_r0[0]) LT 0.1 then sig_model_chiP0[jj] = 1e-5 else $
           sig_model_chiP0[jj] = xcorr_line((Rval[jj]/1e3)/(r0-sig_r0[0]), gamma, vdimless)
        sig_model_chiP1[jj] = xcorr_line((Rval[jj]/1e3)/(r0+sig_r0[1]), gamma, vdimless)
     endfor
     ;x_curvefill, Rval/1e3, sig_model_chiP0, sig_model_chiP1, color=clr.lightgray

     model_chiP = fltarr(nplt)
     for jj=0L,nplt-1 do $
        model_chiP[jj] = xcorr_line((Rval[jj]/1e3)/r0, gamma, vdimless)
     oplot, Rval/1e3, model_chiP, color=clr.cyan

        
     ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     ;; Loop on Radii (Data)
     for kk=0L,nstep-1 do begin
        x1 = Rmin * ((Rstep)^kk)
        x2= Rmin * ((Rstep)^(kk+1))
        
        ;; Data -- Comoving
        gdi = where(qpq_strct.R_phys*(1+qpq_strct.z_fg)*LIT_H GE x1 $
                    and qpq_strct.R_phys*(1+qpq_strct.z_fg)*LIT_H LT x2, npts)
        print, npts
        case ss of
           2: gdsys = where(qpq_strct[gdi].flg_othick EQ 1, QA)
           else: gdsys = where(qpq_strct[gdi].NHI GE logN_min and qpq_strct[gdi].flg_NHI EQ 1 $
                               and qpq_strct[gdi].NHI LT logN_max, QA)
        endcase
        
        ;; Random
        dNdz = igm_calc_lox(fn_strct, qpq_strct[gdi].z_fg, logN_min, LOGN_MAX) * $
               cosm_dxdz(qpq_strct[gdi].z_fg, /w05map)
        dz = 2.0D*(1.0D + qpq_strct[gdi].z_fg)*(vmax/(c.c/1e5))
        QR = total( dNdz * dz )
        
        ;; Stats
        npoiss = x_poisscl(QA,sigma=1)
        xbin = mean([x1,x2])
        
        ;; Evaluate
        chiP = (QA/QR - 1.) > yrng[0]
        chiSig = npoiss/QR - 1.
        
        ;stop
        ;; Plot
        plotsym, 8, 1.2, /fill
        oploterror, [xbin]/1e3, [chiP], [abs(xbin-x1)]/1e3, [(chiP - chiSig[1]) > 1e-5], $
                    psym=8, errcolor=eclr, errthick=5, /lobar, color=eclr
        oploterror, [xbin]/1e3, [chiP], [abs(xbin-x2)]/1e3, [abs(chiP - chiSig[0])], $
                    psym = 8, errcolor = eclr, errthick = 5, /hibar, color = eclr
        print, xbin, chiP
     
        ;; Save
                                ;strct[kk].R1 = x1
                                ;strct[kk].R2 = x2
                                ;strct[kk].npair = npts
                                ;strct[kk].fC = cov_frac
                                ;strct[kk].sig_fC = abs([wilson[2],wilson[0]]-cov_frac)
                                ;strct[kk].fC_IGM = NLLS/float(npts)
        
     endfor

     ;; Label
     xyouts, 0.06, 10.^( alog10(yrng[0]) + 0.15*(alog10(yrng[1])-alog10(yrng[0]))), $
             lbl, color=pclr, charsi=lsz, align=0.
     xyouts, 0.70, 10.^( alog10(yrng[1]) - 0.15*(alog10(yrng[1])-alog10(yrng[0]))), $
             'r!S!d0!R!u'+lbl+'!N = '+string(R0,format='(f4.1)')+' h!u-1!N Mpc', $
             color=pclr, charsi=lsz, align=0.


  endfor
  
  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]

  print, 'fig_trans_corr:  All done!'
  return
end
      
      
