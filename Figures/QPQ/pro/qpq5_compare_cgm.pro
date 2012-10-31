;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  Plots Average/Median EW_Lya vs Physical/Comoving Mpc
pro qpq5_compare_cgm, NO_LINE=no_line

  compile_opt strictarr

  ;; Get structure if necessary
  if not keyword_set( PSFILE ) then psfile = 'qpq5_compare_cgm.ps'
  if not keyword_set( CSZ ) then csz = 2.5
  if not keyword_set( LSZ ) then lsz = 1.3
  if not keyword_set( ASZ ) then asz = 2.0
  if not keyword_set( SYMS ) then syms = 1.2

  if not keyword_set(SIG_FX) then sig_fx = 0.05 ;; Error in tau values
  if not keyword_set( MAXCHI ) then MAXCHI = 1.3
  ;  if not keyword_set( WVMIN ) then wvmin = 700. 
  if not keyword_set(OMEGA_M) then omega_m = 0.3

  if not keyword_set( NRNG ) then nrng = 50L

  cd, getenv('QSO_DIR')+'tex/QPQ5/Figures/pro/', curr=curr
  resolve_routine, 'fig_boot_stack', /compile_full_file
  cd, curr

  ;; Read in pair stuff
  qpq_fil = '~/Dropbox/QSOPairs/qpq5_pairs.fits'
  qpq5_strct = xmrdfits(qpq_fil, 1)
  npair = n_elements(qpq5_strct)

  wrest = [1215.6701, 1334.5323, 1548.195d]
  wrest = [1215.6701, 1334.5323];, 1548.195d]
  nwr = n_elements(wrest)


  NO_LINE = 1
  COVFACT = 1

  ;; Read in LBG
  lls_struct, lbg, '/u/xavier/LLS/Lists/lbg_lls.lst', ROOT='/u/xavier/LLS/', /ew

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; PLOT
  x_psopen, psfile, /portrait
  clr = getcolor(/load)
  !p.multi=[0,1,4]
  fclr = clr.white
  qclr = clr.yellow
  
  xmrg = [7,30]
  ymrg = [0,0]

  xrng = [0., 310]  ; kpc


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
        xtit='R!dphys!N (kpc)'
     endelse

     case ss of
        0: begin
           yrng = [0.0, 3.]     ; Ang
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
     
     
     plot, [0], [0], color=fclr, background=clr.white, charsize=csz,$
           xmargin=xmrg, ymargin=ymrg, xtitle=xtit, xtickn=xspaces, $
           ytitle=ytit, yrange=yrng, thick=4, $
                                ;ytickformat='x_logticks', $
           xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata ;, /ylog

     xlbl = 200.
     xyouts, xlbl, yrng[1]*0.9, lbl, color=fclr,charsi=lsz, align=0.

     ;; Label the Samples
     if ss EQ 1 then begin
        xsmpl = 210.
        ysmpl = 1.0
        ;scl = 1.5
        scl = 0.14
        ;xyouts, xsmpl, ysmpl/(scl^1), 'z~0; Dwarf', color=clr.red, charsi=lsz
        xyouts, xsmpl, ysmpl-(0*scl), 'z~2; QSO', color=qclr, charsi=lsz, align=0.
        xyouts, xsmpl, ysmpl-(1*scl), 'z~2; LBG', color=clr.red, charsi=lsz, align=0.
        xyouts, xsmpl, ysmpl-(2*scl), 'z~0; L*', color=clr.cyan, charsi=lsz, align=0.
        ;xyouts, xsmpl, ysmpl/(scl^0), 'z~2; QSO', color=clr.black, charsi=lsz, align=0.
        ;xyouts, xsmpl, ysmpl/(scl^1), 'z~2; LBG', color=clr.red, charsi=lsz, align=0.
        ;xyouts, xsmpl, ysmpl/(scl^2), 'z~0; L*', color=clr.cyan, charsi=lsz, align=0.
     endif

     ;; Label the Stats
     ;if (ss MOD 2) EQ 0 then slbl = 'Average' else slbl='Median' 
     ;xyouts, 12., yrng[1]/2., slbl, color=clr.black, charsiz=lsz

     ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     ;; z~0 L*
     pclr = clr.cyan

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
     cos_bins = [ [0., 75], [75, 150]]  ; kpc
     sz_cos = size(cos_bins,/dimen)

     ;; Parse
     lambda = wrest[ss] ;1215.6701d
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
        R_lowz = findgen(300.) ; kpc
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
     
     ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     ;; z~0 <0.1L*
     pclr = clr.red

     ;; ;;;;;;;;;;;
     ;; COS-Dwarfs
     ldir = getenv('DROPBOX_DIR')+'/COS-Dwarfs/lowions/'

     restore, ldir+'/cosmetals_megastructure.sav' 

     ;; Bins
     cos_bins = [ [0., 75], [75, 150]]  ; kpc
     sz_cos = size(cos_bins,/dimen)

     ;; Parse
     lambda = wrest[ss] ; 1215.6701d
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
     for ii=0L,sz_cos[1]-1 do begin
        gdR = where( Rperp GT cos_bins[0,ii] and Rperp LT cos_bins[1,ii])
        ;; Stats
        R_lowz[ii] = mean(Rperp[gdR])
        EW_lowz[ii] = mean(EW[gdR]) 
     endfor

     ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     ;; LCO Survey
     if fix(wrest[ss]) EQ 1215 then begin
        lco_bins = [ [150., 300], [300, 500], [500., 1000]] ; kpc
        sz_lco = size(lco_bins, /dimen)
        
        summ_fil = getenv('POVI')+'/Galaxies/Analysis/dwarf_galabs_1Mpc_strct.fits'
        struct = xmrdfits(summ_fil,1)
        ngal = n_elements(struct)
;
        vidx = where(struct.mag[2] GT 0.)
        lco_Rperp = struct[vidx].dra
        lco_EW = struct[vidx].mag[4]/1e3 ; A
        
        Rperp = [Rperp, lco_Rperp]
        EW = [EW, lco_EW]
        
        ;; Co-moving
        ;if ss GT 1 then lco_Rperp = lco_Rperp * (1+struct[vidx].z)

        for ii=0L,sz_lco[1]-1 do begin
           gdR = where( Rperp GT lco_bins[0,ii] and Rperp LT lco_bins[1,ii], ngdR)
           if ngdR EQ 0 then begin
              flg = 0
              print, 'No galaxies!!', lco_bins[*,ii]
           endif else flg = 1
           if flg EQ 0 then continue
           ;; Stats
           R_lowz = [R_lowz, mean(Rperp[gdR])]
           EW_lowz = [EW_lowz, mean(EW[gdR])] 
        endfor
     endif

     ;; Plot
     if keyword_set(DWARFS) then begin
        plotsym, 0, syms, thick=6
        oplot, R_lowz, EW_lowz, color=pclr, psym=8
        if not keyword_set(NO_LINE) then oplot, R_lowz, EW_lowz, color=pclr
     endif
     
     ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     ;; LBG
     pclr = clr.red
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
     if ss EQ 0 then begin
        ;; EW from Rakic+12
        R_LBG = [0.0905958,     0.152221,     0.215018,     0.30372] * 1e3 ; kpc
        EW_LBG = [1.0730036,      0.41661176,      0.65806370,      0.44781657] ; Ang
        sigEW_high = [1.2479758,      0.64263709,      0.98686527,      0.60663553]
        sigEW_low = [0.90953773,      0.17523460,      0.31262733,      0.26702346]
     endif else begin
        
        LBG_bins = [ [0., 100], [100, 300]]
        sz_LBG = size(LBG_bins,/dimen)

        ;; Parse R and EWs
        ion = where(abs(lbg.systems[0].ion[ionZ].state[ioni,*].lambda-wrest[ss]) LT 1e-3, nlbg)
        lbg_idx = ion / 30
        
        ;; EW
        LBG_Rperp = lbg[lbg_idx]. srvy_mag
        LBG_EW = (lbg.systems[0].ion[ionZ].state[ioni,*].clm)[ion]  
        LBG_sigEW = (lbg.systems[0].ion[ionZ].state[ioni,*].sigclm)[ion]  

        ;; Cut on R
        Rperp = LBG_Rperp
        EW = LBG_EW
        R_LBG = fltarr(sz_LBG[1])
        EW_LBG = fltarr(sz_LBG[1])
        for ii=0L,sz_LBG[1]-1 do begin
           gdR = where( Rperp GT LBG_bins[0,ii] and Rperp LT LBG_bins[1,ii], ngdR)
           ;; Stats
           R_LBG[ii] = mean(Rperp[gdR])

           stats = moment(EW[gdR])
           EW_LBG[ii] = stats[0]
           sigEW_LBG[ii] = sqrt(stats[1]/float(ngdR))
        endfor
     endelse
     ;; Plot high-dispersion
     plotsym, 3, syms, /fill
     oplot, R_LBG, EW_LBG, color=pclr, psym=8
     if not keyword_set(NO_LINE) then oplot, R_LBG, EW_LBG, color=pclr

     if fwr EQ 1215 then begin
        oploterror, R_LBG, EW_LBG, [sigEW_high-EW_LBG], color=pclr, errcol=pclr, /hibar, psym=8
        oploterror, R_LBG, EW_LBG, [EW_LBG-sigEW_low], color=pclr, errcol=pclr, /lobar, psym=8
     endif else oploterror, R_LBG, EW_LBG, sigEW_LBG, color=pclr, errcol=pclr, psym=8

     ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     ;; QSOs
     pclr = clr.yellow
     if fix(wrest[ss]) EQ 1215 then begin
        fig_boot_stack, /noplot, EW_vals=ew_vals
        ;QPQ5_Rperp = [71., 142, 251] ; kpc
        ;QPQ5_EW = [3.50, 1.90, 1.80] ; A
        ;;
        R_QSO = EW_vals[2,1:*]
        EW_QSO = EW_vals[0,1:*]
        sigEW_QSO = EW_vals[1,1:*]
     endif else begin ;;  Metals
        QSO_bins = [ [0., 100], [100, 150], [150,300]] ; kpc
        sz_QSO = size(QSO_bins,/dimen)

        ;; Require the lines be outside the Forest
        if not keyword_set(FOREST) then begin ;; Cut out Lya forest
           good = where( (qpq5_strct.z_fg+1)*wrest[ss] GT ((qpq5_strct.z_bg+1)*1215.6701 + 20.), ngd)
           qpq_strct = qpq5_strct[good]
        endif

        ;; Parse R and EWs
        gdew = where( abs(qpq_strct.metal_wrest-wrest[ss]) LT 1e-3 AND $
                      qpq_strct.flg_metal_EW GT 0)
        sz = size(qpq_strct.metal_wrest,/dimen)
        gdi = gdew / sz[0]
        QSO_Rperp = qpq_strct[gdi].R_phys
        QSO_EW = (qpq_strct.raw_metal_EW)[gdew] 

        ;; Cut on R
        Rperp = QSO_Rperp
        EW = QSO_EW
        R_QSO = fltarr(sz_QSO[1])
        EW_QSO = fltarr(sz_QSO[1])
        for ii=0L,sz_QSO[1]-1 do begin
           gdR = where( Rperp GT QSO_bins[0,ii] and Rperp LT QSO_bins[1,ii], ngdR)
           ;; Stats
           R_QSO[ii] = mean(Rperp[gdR])

           stats = moment(EW[gdR])
           EW_QSO[ii] = stats[0]
           sigEW_QSO[ii] = sqrt(stats[1]/float(ngdR))
        endfor

     endelse

     ;; Plot
     plotsym, 0, syms, /fill
     oplot, R_QSO, EW_QSO, color=pclr, psym=8
     if not keyword_set(NO_LINE) then oplot, R_QSO, EW_QSO, color=pclr

     case fwr of
        1215: oploterror, R_QSO, EW_QSO, sigEW_QSO, color=pclr, errcol=pclr, psym=8
        1334: begin
           nQSO = n_elements(R_QSO)
           ;;
           oploterror, R_QSO[0:nQSO-2], EW_QSO[0:nQSO-2], sigEW_QSO[0:nQSO-2], color=pclr, errcol=pclr, psym=8
           ;; Limit
           plotsym, 1, asz, thick=4
           oplot, [R_QSO[nQSO-1]], [EW_QSO[nQSO-1]], psym=8, color=pclr
        end
        1548: oploterror, R_QSO, EW_QSO, sigEW_QSO, color=pclr, errcol=pclr, psym=8
        else:
     endcase
  endfor
  

  ;; Close Ps
  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]

  print, 'fig_compare_cgm:  All done!'
  return
end
      
      
