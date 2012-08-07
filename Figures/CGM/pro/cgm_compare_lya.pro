;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  Plots Average/Median EW_Lya vs Physical/Comoving Mpc
pro cgm_compare_lya

  compile_opt strictarr

  ;; Get structure if necessary
  if not keyword_set( PSFILE ) then psfile = 'cgm_compare_lya.ps'
  if not keyword_set( CSZ ) then csz = 2.5
  if not keyword_set( LSZ ) then lsz = 2.3
  if not keyword_set( ASZ ) then asz = 2.0

  if not keyword_set(SIG_FX) then sig_fx = 0.05 ;; Error in tau values
  if not keyword_set( MAXCHI ) then MAXCHI = 1.3
  ;  if not keyword_set( WVMIN ) then wvmin = 700. 
  if not keyword_set(OMEGA_M) then omega_m = 0.3

  if not keyword_set( NRNG ) then nrng = 50L

  ;; Read in QPQ stuff
  qpq_fil = '~/Dropbox/QSOPairs/qpq5_pairs.fits'
  qpq5_strct = xmrdfits(qpq_fil, 1)
  npair = n_elements(qpq5_strct)

  cd, getenv('QSO_DIR')+'tex/QPQ5/Figures/pro/', curr=curr
  resolve_routine, 'fig_boot_stack', /compile_full_file
  cd, curr

  wrest = 1215.6701

  COVFACT = 1

  ;; Read in LBG
  lls_struct, lbg, '/u/xavier/LLS/Lists/lbg_lls.lst', ROOT='/u/xavier/LLS/', /ew

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; PLOT
  x_psopen, psfile, /maxs
  clr = getcolor(/load)
  lclr = clr.white
  !p.multi=[0,1,1]
  
  xmrg = [8,5]
  ymrg = [4,1]

  xrng = [0., 310]  ; kpc

  for qq=0L,4 do begin

     ;; Label
     fwr = fix(wrest)
     ytit='W!dLy!9a!X!N (Ang)'  
     xtit='R!dphys!N (kpc)'
     yrng = [0.0, 3.]           ; Ang

     plot, [0], [0], color=lclr, background=clr.white, charsize=csz,$
           xmargin=xmrg, ymargin=ymrg, xtitle=xtit, xtickn=xspaces, $
           ytitle=ytit, yrange=yrng, thick=4, $
           xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata ;, /ylog

     xsmpl = 230.
     ysmpl = 2.7
     scl = 0.25
     xyouts, xsmpl, ysmpl-(0*scl), 'z~2; LBG', color=clr.red, charsi=lsz, align=0.
     if qq GT 1 then xyouts, xsmpl, ysmpl-(1*scl), 'z~0; L*', color=clr.cyan, charsi=lsz, align=0.
     if qq GT 2 then xyouts, xsmpl, ysmpl-(2*scl), 'z~2; QSO', color=clr.yellow, charsi=lsz, align=0.
     if qq GT 3 then xyouts, xsmpl, ysmpl-(3*scl), 'z~0; Dwarfs', color=lclr, charsi=lsz, align=0.

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
     if qq EQ 0 then begin 
        restore, all_mega_fil 
        lstar_mega = megastruct
     endif
        
     ;;  Cut on luminosity and R
     megastruct = coshalos_cut_mega(lstar_mega)

     ;; Bins
     cos_bins = [ [0., 75], [75, 150]]  ; kpc
     sz_cos = size(cos_bins,/dimen)

     ;; Parse
     lambda = wrest ;1215.6701d
     getion, lambda, ioni, elm, Z=ionz, NM=strioni
     ewstrct = mega_parse_trans(lambda, mega=megastruct)
     gdew = where(ewstrct.flgEW GT 0)

     cos_Rperp = megastruct[gdew].rhofinal 
     cos_EW = ewstrct[gdew].EW / 1e3  ;; Ang

     ;; Co-Moving?
     
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
     if fwr EQ 1215 then begin
        lco_bins = [ [150., 300] , [300, 500]];, [500., 1000]] ; kpc
        sz_lco = size(lco_bins, /dimen)
        ;sz_lco = [2,1]
        
        summ_fil = getenv('POVI')+'/Galaxies/Analysis/lstar_galabs_1Mpc_strct.fits'
        struct = xmrdfits(summ_fil,1)
        ngal = n_elements(struct)
        
        vidx = where(struct.mag[2] GT 0.)
        lco_Rperp = struct[vidx].dra
        lco_EW = struct[vidx].mag[4]/1e3 ; A
        
        Rperp = [Rperp, lco_Rperp]
        EW = [EW, lco_EW]
        
        ;; Co-moving
        
        for ii=0L,sz_lco[1]-1 do begin
           gdR = where( Rperp GT lco_bins[0,ii] and Rperp LT lco_bins[1,ii])
           ;; Stats
           R_lowz = [R_lowz, mean(Rperp[gdR])]
           EW_lowz = [EW_lowz, mean(EW[gdR])] 
        endfor
     endif

     ;; Overwrite CIV from Chen et al. 2001
     if fwr EQ 1548 then begin
        R_lowz = findgen(300.) ; kpc
        EW_lowz = 2* 0.006 * sqrt( (96./0.72)^2 - R_lowz^2) / 2. ; Full doublet
     endif

     ;; Plot
     if qq GT 1 then begin
        plotsym, 8, 1.5, /fill
        if fwr EQ 1548 then begin
           oplot, R_lowz, EW_lowz, color=pclr, linesty=1
        endif else begin
           oplot, R_lowz, EW_lowz, color=pclr, psym=8
           oplot, R_lowz, EW_lowz, color=pclr
        endelse
     endif
     
     ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     ;; z~0 <0.1L*
     pclr = lclr

     ;; ;;;;;;;;;;;
     ;; COS-Dwarfs
     ldir = getenv('DROPBOX_DIR')+'/COS-Dwarfs/lowions/'
     if qq EQ 0 then begin
        restore, ldir+'/cosmetals_megastructure.sav' 
        dwarf_mega = megastruct
     endif
     megastruct = dwarf_mega

     ;; Bins
     cos_bins = [ [0., 75], [75, 150]]  ; kpc
     sz_cos = size(cos_bins,/dimen)

     ;; Parse
     lambda = wrest ; 1215.6701d
     getion, lambda, ioni, elm, Z=ionz, NM=strioni
     ewstrct = mega_parse_trans(lambda, mega=megastruct)
     gdew = where(ewstrct.flgEW GT 0)

     cos_Rperp = megastruct[gdew].rhofinal 
     cos_EW = ewstrct[gdew].EW / 1e3  ;; Ang

     ;; Co-Moving?
     
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
     if fwr EQ 1215 then begin
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

     if qq GT 3 then begin
        ;; Plot
        plotsym, 0, 1.5, thick=6
        oplot, R_lowz, EW_lowz, color=pclr, psym=8
        oplot, R_lowz, EW_lowz, color=pclr
     endif
     
     ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     ;; LBG
     pclr = clr.red
     s10_Rperp = [31., 63, 103] ; kpc
     case fwr of
        1215: s10_EW = [2.01, 1.23, 0.92]  ; A
        1334: s10_EW = [2.13, 1.18, 0.13]/2. ; A
        1548: s10_EW = [0.90, 0.67, 0.12] ; A
        else: stop
     endcase

     ;; Need to add in Olivera stuff here!!


     R_lbg = s10_Rperp
     EW_lbg = s10_EW

     ;; Plot
     plotsym, 3, 1.0, thick=4
     oplot, R_lbg, EW_lbg, color=pclr, psym=8
     oplot, R_lbg, EW_lbg, color=pclr, linesty=1, thick=4

     ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     ;; LBG from high dispersion
     if fwr EQ 1215 then begin
        ;; EW from Rakic+12
        R_LBG = [0.0905958,     0.152221,     0.215018,     0.30372] * 1e3 ; kpc
        EW_LBG = [1.0730036,      0.41661176,      0.65806370,      0.44781657] ; Ang
     endif else begin
        
        LBG_bins = [ [0., 100], [100, 300]]
        sz_LBG = size(LBG_bins,/dimen)

        ;; Parse R and EWs
        ion = where(abs(lbg.systems[0].ion[ionZ].state[ioni,*].lambda-wrest) LT 1e-3, nlbg)
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
           gdR = where( Rperp GT LBG_bins[0,ii] and Rperp LT LBG_bins[1,ii])
           ;; Stats
           R_LBG[ii] = mean(Rperp[gdR])
           EW_LBG[ii] = mean(EW[gdR]) 
        endfor
     endelse
     ;; Plot high-dispersion
     if qq GT 0 then begin
        plotsym, 3, 1.2, /fill
        oplot, R_LBG, EW_LBG, color=pclr, psym=8
        oplot, R_LBG, EW_LBG, color=pclr
     endif

     ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     ;; QSOs
     pclr = clr.yellow
     if fwr EQ 1215 then begin
        if qq EQ 0 then fig_boot_stack, /noplot, EW_vals=ew_vals
        ;;
        R_QSO = EW_vals[2,1:*]
        EW_QSO = EW_vals[0,1:*]
        sigEW_QSO = EW_vals[1,1:*]
     endif else begin ;;  Metals
        QSO_bins = [ [0., 100], [100, 150], [150,300]] ; kpc
        sz_QSO = size(QSO_bins,/dimen)

        ;; Require the lines be outside the Forest
        if not keyword_set(FOREST) then begin ;; Cut out Lya forest
           good = where( (qpq5_strct.z_fg+1)*wrest GT ((qpq5_strct.z_bg+1)*1215.6701 + 20.), ngd)
           qpq_strct = qpq5_strct[good]
        endif

        ;; Parse R and EWs
        gdew = where( abs(qpq_strct.metal_wrest-wrest) LT 1e-3 AND $
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
           gdR = where( Rperp GT QSO_bins[0,ii] and Rperp LT QSO_bins[1,ii])
           ;; Stats
           R_QSO[ii] = mean(Rperp[gdR])
           EW_QSO[ii] = mean(EW[gdR]) 
        endfor

     endelse

     if qq GT 2 then begin
        ;; Plot
        plotsym, 0, 1.5, /fill
        oplot, R_QSO, EW_QSO, color=pclr, psym=8
        oplot, R_QSO, EW_QSO, color=pclr
        if fwr EQ 1215 then oploterror, R_QSO, EW_QSO, sigEW_QSO, color=pclr, errcol=pclr
        
        ;; Limit on CII
        if fwr EQ 1334 then begin
           nQSO = n_elements(R_QSO)
           plotsym, 1, asz, thick=4
           oplot, [R_QSO[nQSO-1]], [EW_QSO[nQSO-1]], psym=8, color=pclr
        endif
     endif

  endfor
  

  ;; Close Ps
  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]

  print, 'fig_compare_cgm:  All done!'
  return
end
      
      
