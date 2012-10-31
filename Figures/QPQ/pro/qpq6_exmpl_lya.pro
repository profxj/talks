;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  Estimates and plots MFP with bootstrap
;;      Writes the fits to disk (if BOOT_STRCT_FIL not set)
;;
;; fig_mfp_obs, BOOT_STRCT_FIL='../Analysis/boot_mfp_wfc3.fits'
pro qpq6_exmpl_lya

  compile_opt strictarr

  cd, getenv('QSO_DIR')+'tex/QPQ6/Analysis/AbsLines/pro/', curr=curr
  resolve_routine, 'qpq6_mkdata', /compile_full_file, /is_function
  resolve_routine, 'qpq6_grab_idlines', /compile_full_file, /is_function
  cd, curr

  ;; Get structure if necessary
  if not keyword_set( PSFILE ) then psfile = 'qpq6_exmpl_lya.ps'
  if not keyword_set( CSZ ) then csz = 1.9
  if not keyword_set( LSZ ) then lsz = 1.6

  ;; Read in full pair info
;  savefile = '~/Dropbox/QSOPairs/QPQ_00.000-00.300_Thu-Aug-25-14:15:13-2011.sav' 
;  restore, savefile

  qsodir = getenv('QSO_DIR')

  ;; Read in pair structure
  qpq_fil = '~/Dropbox/QSOPairs/qpq6_final_cut.fits'
  qpq_strct = xmrdfits(qpq_fil, 1)

  ;gdp = where(qpq_strct.R_phys LT 300. AND $
  ;            qpq_strct.flg_EWLya GT 0, npair)  
  gd_strct = qpq_strct

  ;qso_plt = ['J002126.10-025222.0', $  ;; f/g QSO
  ;           'J022517.68+004821.9', $
  ;           'J023139.53+001758.4', $
  ;           'J025038.68-004739.1', $
  ;           'J074031.15+224616.1', $
  ;           'J081523.63+494011.9']
  qso_plt = ['BOSSJ0010-0051', $  ; 387kpc
             'BOSSJ0224-0539'] ; 315kpc
  nplt = n_elements(qso_plt)

  ;; ;;;;;;;;;;;;;;;;;
  ;; PLOT
  x_psopen, psfile, /portrait
  clr = getcolor(/load)
  fclr = clr.white
  !p.multi=[0,1,2]
  nwind = 6
  
  xtit='Wavelength (Ang)'
  ytit='Relative Flux'
  xmrg = [6.5,0.5]
  ymrg = [3.5,0.2]

  vmnx = [-2000., 2000]

  ;; Loop
  for qq=0L,nplt-1 do begin

     mt = where(strmatch(strtrim(gd_strct.qso,2), qso_plt[qq]), nmt)
     if nmt NE 1 then stop

;     print, qq, gd_strct[qq].qso, gd_strct[qq].lya_fil
     xrng = 1215.6701*(1+gd_strct[mt].z_fg) *(1+vmnx/3e5)
     xrng[0] = xrng[0] < (gd_strct[mt].wvmnx[0]-3.)
     xrng[1] = xrng[1] > (gd_strct[mt].wvmnx[1]+3.)

     ;; Read Spectrum
     if (gd_strct[mt].inflg_lyafil EQ 2) or (gd_strct[mt].inflg_lyafil EQ 11) then $
        auto = 0 else auto = 1
;     stop
     lyafil = strtrim(qsodir+gd_strct[mt].lya_fil,2)
     a = findfile(lyafil+'*', count=nfil)
     data = qpq6_mkdata(gd_strct[mt], flg=flg_data)
     flux_bg = data.flux_bg
     sig_bg = data.sig_bg
     wave_bg = data.wave

     ;if gd_strct[mt].qso EQ 'SDSSJ0225+0048' then stop

     IF djs_median(flux_bg LE 1d-10) THEN BEGIN
        scl = 1d17
        flux_bg = scl*flux_bg
        sig_bg = scl*sig_bg
     ENDIF ELSE scl = 1.

     ;;
     gdp = where(wave_bg GT xrng[0] and wave_bg LT xrng[1], ngdp)
     srt = sort(flux_bg[gdp])
     mx = flux_bg[gdp[srt[(round(0.9*ngdp)) < (ngdp-1)]]]

     ;; Voigt?  If so, use that continuum and
     pos2 = strpos(lyafil, '.fits')
     pos1 = strpos(lyafil, '/', /reverse_sear)
     voigt_fil = '../Analysis/N_HI/Lya_Fits/'+strmid(lyafil, pos1+1, pos2-pos1-1)+'_Lyafit.idl'
     fil = file_search(voigt_fil, count=nvoigt)
     if nvoigt EQ 1 then begin
        restore, fil
        cmx = max(conti[gdp])
     endif else cmx = 0.
     
     yrng = [-0.05, 1.4] * mx
     yrng[1] = yrng[1] > cmx*1.2

     plot, [0], [0], color=fclr, background=clr.white, charsize=csz,$
           xmargin=xmrg, ymargin=ymrg, xtitle=xtit, $
           ytitle=ytit, yrange=yrng, thick=4, $
           xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata

     ;; Shade the redshift uncertainty
     z_pm = x_relvel(gd_strct[mt].z_fg, [1,-1]*gd_strct[mt].z_fsig)
     x_curvefill, 1215.6701*(1+z_pm), replicate(yrng[0],2), replicate(yrng[1],2), $
                  color=clr.darkgray

     ;; Plot Continuum
     if nvoigt EQ 1 then $
        oplot, data.wave, CONTI, color=clr.cyan, linesty=2, thick=4 $
     else $
        oplot, data.wave, data.conti*SCL, color=clr.cyan, linesty=2, thick=4

     ;; Mark positions
;     oplot, replicate(1215.6701*(1+gd_strct[mt].z_fg),2), $
;            yrng, color=clr.blue, linesty=3
     ;oplot, replicate(1215.6701*(1+gd_strct[mt].z_Lya),2), $
     ;       yrng, color=clr.red, linesty=2
     ;oplot, replicate(gd_strct[qq].wvmnx[0],2), yrng, color=clr.orange, linesty=1
     ;oplot, replicate(gd_strct[qq].wvmnx[1],2), yrng, color=clr.orange, linesty=1

     ;; Velocity range
     wvmnx = gd_strct[mt].wvmnx
     ;oplot, replicate(wvmnx[0],2), yrng, color=clr.orange, linesty=1
     ;oplot, replicate(wvmnx[1],2), yrng, color=clr.orange, linesty=1
     EW = gd_strct[mt].EWlya

     ;; Voigt profile
     if nvoigt EQ 1 then begin
        fx = x_voigt(data.wave, lines, FWHM=4.)
        oplot, data.wave, fx*CONTI, color=clr.cyan, thick=2
     endif

     ;; Data
     oplot, wave_bg, flux_bg, color=fclr, psym=10, thick=5
     ;oplot, wave_bg, sig_bg, color=clr.green, psym=10, thick=2, linesty=1

     ;; Label
     x_radec, ras, decs, gd_strct[mt].rad, gd_strct[mt].decd, /flip
     rav = strsplit(ras, ':', /extract)
     decv = strsplit(decs, ':', /extract)
     qlbl = 'J'+rav[0]+rav[1]+rav[2]+decv[0]+decv[1]+decv[2]

     xyouts, xrng[0]+(xrng[1]-xrng[0])*0.5, yrng[1]-(yrng[1]-yrng[0])*0.08, $
             qlbl, color=fclr, charsiz=lsz, align=0.5

     if gd_strct[mt].flg_NHI EQ 3 then c_NHI = '<' else c_NHI = '='
     xyouts, xrng[0]+(xrng[1]-xrng[0])*0.5, yrng[1]-(yrng[1]-yrng[0])*0.16, $
             'R!dphys!N = '+string(round(gd_strct[mt].R_phys),format='(i3)')+$  
             'kpc;    z!dfg!N = '+string(gd_strct[mt].z_fg,format='(f5.3)'), $ 
             color=fclr, charsiz=lsz, align=0.5

  endfor

  ;; Close Ps
  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]

  print, 'fig_exmpl_lya:  All done!'
       
  return
end
      
      
