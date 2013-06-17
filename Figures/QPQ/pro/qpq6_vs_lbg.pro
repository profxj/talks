;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  Plots Average/Median EW_Lya vs Physical/Comoving Mpc
pro qpq6_vs_lbg, NO_LINE=no_line, FSTRCT=fstrct

  compile_opt strictarr

  ;; Get structure if necessary
  if not keyword_set( PSFILE ) then psfile = 'qpq6_vs_lbg.ps'
  if not keyword_set( CSZ ) then csz = 2.2
  if not keyword_set( LSZ ) then lsz = 1.9
  if not keyword_set( ASZ ) then asz = 2.0
  if not keyword_set( SYMS ) then syms = 1.2

  if not keyword_set(SIG_FX) then sig_fx = 0.05 ;; Error in tau values
  if not keyword_set( MAXCHI ) then MAXCHI = 1.3
  ;  if not keyword_set( WVMIN ) then wvmin = 700. 
  if not keyword_set(OMEGA_M) then omega_m = 0.3

  if not keyword_set( NRNG ) then nrng = 50L

  compile_opt strictarr

  cd, getenv('QSO_DIR')+'tex/QPQ6/Figures/pro/', curr=curr
  resolve_routine, 'fig_ewstack_vs_rho', /compile_full_file 
  cd, curr


  ;; Read in pair stuff
  qpq_fil = '~/Dropbox/QSOPairs/qpq5_pairs.fits'
  qpq5_strct = xmrdfits(qpq_fil, 1)
  npair = n_elements(qpq5_strct)

  wrest = [1215.6701]
  nwr = n_elements(wrest)


  COVFACT = 1
  NO_LINE = 1

  ;; Read in LBG
  lls_struct, lbg, '/u/xavier/LLS/Lists/lbg_lls.lst', ROOT='/u/xavier/LLS/', /ew

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; PLOT
  x_psopen, psfile, /maxs
  clr = getcolor(/load)
  !p.multi=[0,1,1]

  qclr = clr.white
  lclr = clr.cyan
  
  xmrg = [7,15]
  ymrg = [5,1]

  xrng = [10., 1000.]  ; kpc


  for ss=0,nwr-1 do begin
     
     ;; Label
     fwr = fix(wrest[ss])
     if fwr EQ 1215 then ytit='W!dLy!9a!X!N (Ang)'  $
     else ytit='W!d'+strtrim(fwr,2)+'!N (Ang)'

     if ss LT (nwr-1) then begin
        xspaces = replicate(' ',30) 
        xtit = ''
     endif else begin
        delvarx, xspaces
        thisletter = byte(94)
        perpletter = '!9' + string(thisletter) + '!X'
        xtit='R!d'+perpletter+'!N (kpc)'
     endelse

     case ss of
        0: begin
           yrng = [0.0, 4.]     ; Ang
           lbl = '(a) Ly!9a!X'
        end
        1: begin
           yrng = [0.0, 1.3]    ; Ang
           lbl = '(b) CII 1334'
        end
        2: begin
           yrng = [0.0, 0.95]   ; Ang
           lbl = '(c) CIV 1548'
        end
        else: stop
     endcase
     
     
     plot, [0], [0], color=clr.white, background=clr.white, charsize=csz,$
           xmargin=xmrg, ymargin=ymrg, xtitle=xtit, xtickn=xspaces, $
           ytitle=ytit, yrange=yrng, thick=4, $
           xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata, /xlog

     xlbl = 200.
     ;xyouts, xlbl, yrng[1]*0.9, lbl, color=clr.black,charsi=lsz, align=0.

     ;; Label the Samples
     if ss EQ 0 then begin
        xsmpl = 300.
        ysmpl = 3.5
        ;scl = 1.5
        scl = 0.34
        ;xyouts, xsmpl, ysmpl/(scl^1), 'z~0; Dwarf', color=clr.red, charsi=lsz
        xyouts, xsmpl, ysmpl-(0*scl), 'z~2; QSO', color=qclr, charsi=lsz, align=0.
        xyouts, xsmpl, ysmpl-(1*scl), 'z~2; LBG', color=lclr, charsi=lsz, align=0.
        ;xyouts, xsmpl, ysmpl-(2*scl), 'z~0; L*', color=clr.blue, charsi=lsz, align=0.
        ;xyouts, xsmpl, ysmpl/(scl^0), 'z~2; QSO', color=clr.black, charsi=lsz, align=0.
        ;xyouts, xsmpl, ysmpl/(scl^1), 'z~2; LBG', color=clr.red, charsi=lsz, align=0.
        ;xyouts, xsmpl, ysmpl/(scl^2), 'z~0; L*', color=clr.blue, charsi=lsz, align=0.
     endif

     ;; Label the Stats
     ;if (ss MOD 2) EQ 0 then slbl = 'Average' else slbl='Median' 
     ;xyouts, 12., yrng[1]/2., slbl, color=clr.black, charsiz=lsz

     ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     ;; z~0 L*
     if keyword_set(z0) then begin
        pclr = clr.blue
        
        ;; Initialize
        cd, getenv('PCOSHALOS')+'Analysis/pro/', curr=curr
        RESOLVE_ROUTINE, 'coshalos_lowmetals_initparm', /compile_full, /either
        RESOLVE_ROUTINE, 'coshalos_cut_mega', /compile_full, /either
        RESOLVE_ROUTINE, 'mega_parse_trans', /compile_full, /either
        coshalos_lowmetals_initparm, coshalos_init
        cd, curr
        
        ;; ;;;;;;;;;;;
        ;; COS-Halos
        ldir = getenv('DROPBOX_DIR')+'/COS-Halos/lowions/'
        
        all_mega_fil = ldir+coshalos_init.all_mega_fil 
        restore, all_mega_fil 
        
        ;;  Cut on luminosity and R
        megastruct = coshalos_cut_mega(megastruct)
        
        ;; Bins
        cos_bins = [ [0., 75], [75, 150]] ; kpc
        sz_cos = size(cos_bins,/dimen)
        
        ;; Parse
        lambda = wrest[ss]      ;1215.6701d
        getion, lambda, ioni, elm, Z=ionz, NM=strioni
        ewstrct = mega_parse_trans(lambda, mega=megastruct)
        gdew = where(ewstrct.flgEW GT 0)
        
        cos_Rperp = megastruct[gdew].rhofinal 
        cos_EW = ewstrct[gdew].EW / 1e3  ;; Ang
        
        ;; Co-Moving?
                                ;if ss GT 1 then cos_Rperp = cos_Rperp * (1+megastruct[gdew].zfinal) 
        
        Rperp = cos_Rperp
        EW = cos_EW
        ;; Loop on bins
        R_lowz = fltarr(sz_cos[1])
        EW_lowz = fltarr(sz_cos[1])
        sigEW_lowz = fltarr(sz_cos[1])
        for ii=0L,sz_cos[1]-1 do begin
           gdR = where( Rperp GT cos_bins[0,ii] and Rperp LT cos_bins[1,ii], ngdR)
           ;; Stats
           R_lowz[ii] = mean(Rperp[gdR])
           
           stats = moment(EW[gdR])
           EW_lowz[ii] = stats[0]
           sigEW_lowz[ii] = sqrt(stats[1]/float(ngdR))
        endfor
        
        ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        ;; LCO Survey
        if fix(wrest[ss]) EQ 1215 then begin
           lco_bins = [ [150., 300]] ;, [300, 500], [500., 1000]] ; kpc
                                ;sz_lco = size(lco_bins, /dimen)
           sz_lco = [2,1]
           
           summ_fil = getenv('POVI')+'/Galaxies/Analysis/lstar_galabs_1Mpc_strct.fits'
           struct = xmrdfits(summ_fil,1)
           ngal = n_elements(struct)
           
           vidx = where(struct.mag[2] GT 0.)
           lco_Rperp = struct[vidx].dra
           lco_EW = struct[vidx].mag[4]/1e3 ; A
           
           Rperp = [Rperp, lco_Rperp]
           EW = [EW, lco_EW]
           
           ;; Co-moving
                                ;if ss GT 1 then lco_Rperp = lco_Rperp * (1+struct[vidx].z)
           
           for ii=0L,sz_lco[1]-1 do begin
              gdR = where( Rperp GT lco_bins[0,ii] and Rperp LT lco_bins[1,ii], ngdR)
              ;; Stats
              R_lowz = [R_lowz, mean(Rperp[gdR])]
              stats = moment(EW[gdR])
              EW_lowz = [EW_lowz, stats[0]]
              sigEW_lowz = [sigEW_lowz, sqrt(stats[1]/float(ngdR))]
           endfor
        endif
        
        ;; Overwrite CIV from Chen et al. 2001
        if fix(wrest[ss]) EQ 1548 then begin
           R_lowz = findgen(300.)                                ; kpc
           EW_lowz = 2* 0.006 * sqrt( (96./0.72)^2 - R_lowz^2) / 2. ; Full doublet
        endif
        
        ;; Plot
        plotsym, 8, syms, /fill
        if fix(wrest[ss]) EQ 1548 then begin
           oplot, R_lowz, EW_lowz, color=pclr, linesty=1
        endif else begin
           oplot, R_lowz, EW_lowz, color=pclr, psym=8
           if not keyword_set(NO_LINE) then oplot, R_lowz, EW_lowz, color=pclr
           oploterror, R_lowz, EW_lowz, sigEW_lowz, color=pclr, errcol=pclr, psym=8
        endelse
     endif
     
     ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     ;; LBG
     pclr = lclr
     s10_Rperp = [31., 63, 103] ; kpc
     case fix(wrest[ss]) of 
        1215: begin
           s10_EW = [2.01, 1.23, 0.92]       ; A
           s10_sigEW = [0.15, 0.20, 0.12]       ; A
        end
        1334: begin
           s10_EW = [0.90, 0.67, 0.12]       ; A
           s10_sigEW = [0.08, 0.12, 99.99]       ; A
        end
        1548: begin
           s10_EW = [2.13, 1.18, 0.13]/2. ; A
           s10_sigEW = [0.15, 0.15, 0.05]/2. ; A
        end
        else: stop
     endcase

     ;; Need to add in Olivera stuff here!!

     ;if ss GT 1 then s10_Rperp = s10_Rperp * (1 + 2.3)

     R_lbg = s10_Rperp
     EW_lbg = s10_EW
     sigEW_lbg = s10_sigEW

     ;; Plot
     plotsym, 3, 0.8
     oplot, R_lbg, EW_lbg, color=pclr, psym=8
     if not keyword_set(NO_LINE) then oplot, R_lbg, EW_lbg, color=pclr, linesty=1, thick=2
     gdE = where(sigEW_lbg LT 99, complement=limit, ncomplement=nlimit)
     oploterror, [R_lbg[gdE]], [EW_lbg[gdE]], [sigEW_lbg[gdE]], color=pclr, errcol=pclr, errsty=1, thick=2, psym=8
     if nlimit GT 0 then begin
        plotsym, 1, asz, thick=4
        oplot, [R_LBG[limit]], [EW_LBG[limit]], psym=8, color=pclr
     endif

     ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     ;; LBG from high dispersion
     ;; EW from Rakic+12
     R_LBG = [0.0905958,     0.152221,     0.215018,     0.30372,  0.429, $ 
              0.606002,     0.856001,      1.20913] * 1e3                           ; kpc
     EW_LBG = [1.0730036,      0.41661176,      0.65806370,      0.44781657, $    ; Ang
               0.42758677,      0.39706513,      0.45419831,      0.40758605]
     sigEW_high = [1.2479758,      0.64263709,      0.98686527,      0.60663553, $
                   0.53847043,      0.48711898,      0.51980822,      0.45615305]
     sigEW_low = [0.90953773,      0.17523460,      0.31262733,      0.26702346, $
                  0.30904305,      0.28374102,      0.37978157,      0.35363075]
     ;; Plot high-dispersion
     plotsym, 3, syms, /fill
     oplot, R_LBG, EW_LBG, color=pclr, psym=8
     if not keyword_set(NO_LINE) then oplot, R_LBG, EW_LBG, color=pclr

     if fwr EQ 1215 then begin
        oploterror, R_LBG, EW_LBG, [sigEW_high-EW_LBG], color=pclr, errcol=pclr, /hibar, psym=8
        oploterror, R_LBG, EW_LBG, [EW_LBG-sigEW_low], color=pclr, errcol=pclr, /lobar, psym=8
     endif else oploterror, R_LBG, EW_LBG, sigEW_LBG, color=pclr, errcol=pclr, psym=8

     ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     ;; QSOs
     pclr = qclr
     if not keyword_set(FSTRCT) then begin
        fig_ewstack_vs_rho, FSTRCT=qpq_str, /NOPS
        fstrct = qpq_str
     endif else qpq_str = fstrct
        ;;
     R_QSO = qpq_str.Rval 
     EW_QSO = qpq_str.EWval
     sigEW_QSO = qpq_str.sigEW

     ;; Plot
     plotsym, 0, syms, /fill
     ;oplot, R_QSO, EW_QSO, color=pclr, psym=8
     ;if not keyword_set(NO_LINE) then oplot, R_QSO, EW_QSO, color=pclr
     oploterror, R_QSO, EW_QSO, sigEW_QSO, color=pclr, errcol=pclr, psym=8
  
  endfor

  ;; Close Ps
  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]

  print, 'fig_compare_cgm:  All done!'
  return
end
      
      
