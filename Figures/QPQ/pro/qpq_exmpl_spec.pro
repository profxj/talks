;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  Estimates and plots MFP with bootstrap
;;      Writes the fits to disk (if BOOT_STRCT_FIL not set)
;;
;; fig_mfp_obs, BOOT_STRCT_FIL='../Analysis/boot_mfp_wfc3.fits'
pro qpq_exmpl_spec

  compile_opt strictarr

  cd, getenv('QSO_DIR')+'tex/QPQ6/Analysis/AbsLines/pro/', curr=curr
  resolve_routine, 'qpq6_mkdata', /compile_full_file, /is_function
  resolve_routine, 'qpq6_grab_idlines', /compile_full_file, /is_function
  cd, curr
  cd, getenv('QSO_DIR')+'tex/QPQ6/Analysis/pro/', curr=curr
  resolve_routine, 'qpq6_get_instrument', /compile_full_file, /is_function
  cd, curr

  ;; Get structure if necessary
  if not keyword_set( PSFILE ) then psfile = 'qpq_exmpl_spec.ps'
  if not keyword_set( CSZ ) then csz = 2.3
  if not keyword_set( LSZ ) then lsz = 1.3

  ;; Read in full pair info
;  savefile = '~/Dropbox/QSOPairs/QPQ_00.000-00.300_Thu-Aug-25-14:15:13-2011.sav' 
;  restore, savefile

  qsodir = getenv('QSO_DIR')

  ;; Read in pair structure
  qpq_fil = '~/Dropbox/QSOPairs/qpq6_pairs.fits'
  qpq_strct = xmrdfits(qpq_fil, 1)

  gdp = where(qpq_strct.R_phys LT 300. AND $
              qpq_strct.flg_EWLya GT 0, npair)  
  gd_strct = qpq_strct[gdp]

;  fg_qsos = ['J002126.10-025222.0', $  ;; BOSS
;             'J002801.18−104933.9', $ ;; ESI
;             'J003423.06−105002.0', $ ;; LRISb
;             'J081420.37+325016.1', $ ;; GMOS
;             'J085358.36−001108.0']  ;; MagE
  fg_qsos = ['BOSSJ0021-0252', $, ;; BOSS
             'SDSSJ0028-1049', $ ;; ESI
             'SDSSJ0034-1050', $ ;; LRISb
             'APOJ0814+3250', $ ;; GMOS
             'APOJ0853-0011' ]   ;;  MagE
  ymx = [8., 0.6, 12, 17, 5]
  lflg = [0, 0, 1, 1, 0]  ;; 0=left, 1=right
  apos = [5.5, 0.3, 4.5, 13,  2.0]
  nplt = n_elements(fg_qsos)

  npx = 1
  npy = nplt
  ;; ;;;;;;;;;;;;;;;;;
  ;; PLOT
  x_psopen, psfile, /portrait
  clr = getcolor(/load)
  ;lclr = clr.lightgray
  lclr = clr.white
  nwind = 6
  !p.multi=[npx*(npy+1),npx,npy+1,0,1]
  
  xmrg = [7.5,0.5]
  ymrg = [0.,0.]

  xrng = [3200., 5300] ;; Ang

  ;; Loop
  cnt = 0L
  for qq=0L,npair-1 do begin

     ;; QSO ?
     if total(strmatch(fg_qsos, strtrim(gd_strct[qq].qso,2))) EQ 0 then continue
     print, cnt, gd_strct[qq].qso

     if (cnt mod npy) EQ 0 then begin
        if cnt GT 0 then offy = (nplt MOD 2) else offy = 0
        !p.multi=[npx*(npy+1) - (npy+1)*(cnt/npy) - offy, $
                  npx,npy+1,0,1]
     endif

     ;; RA/DEC
     x_radec, ras, decs, gd_strct[qq].rad, gd_strct[qq].decd, /flip
     rav = strsplit(ras, ':', /extract)
     decv = strsplit(decs, ':', /extract)
     qlbl = 'J'+rav[0]+rav[1]+rav[2]+decv[0]+decv[1]+decv[2]

     ;; Read Spectrum
     if (gd_strct[qq].inflg_lyafil EQ 2) or (gd_strct[qq].inflg_lyafil EQ 11) then $
        auto = 0 else auto = 1
;     stop
     lyafil = strtrim(qsodir+gd_strct[qq].lya_fil,2)
     a = findfile(lyafil+'*', count=nfil)
     data = qpq6_mkdata(gd_strct[qq], flg=flg_data)
     flux_bg = data.flux_bg
     sig_bg = data.sig_bg
     wave_bg = data.wave

     ;if gd_strct[qq].qso EQ 'SDSSJ0225+0048' then stop

     IF djs_median(flux_bg LE 1d-10) THEN BEGIN
        scl = 1d17
        flux_bg = scl*flux_bg
        sig_bg = scl*sig_bg
     ENDIF ELSE scl = 1.

     gdp = where(wave_bg GT xrng[0] and wave_bg LT xrng[1], ngdp)
     srt = sort(flux_bg[gdp])
     mx = flux_bg[gdp[srt[(round(0.95*ngdp)) < (ngdp-1)]]]

     yrng = [-0.05, ymx[cnt]]

     if (cnt MOD npy) NE (npy-1) and (cnt NE nplt-1) then begin
        xspaces = replicate(' ',30) 
        xtit = ''
     endif else begin 
        if keyword_set(XSPACES) then delvarx, xspaces
        xtit='Wavelength (Ang)'
     endelse
     if (cnt GT npy-1) then yspaces = replicate(' ',30) 

     plot, [0], [0], color=lclr, background=clr.white, charsize=csz,$
           xmargin=xmrg, ymargin=ymrg, xtitle=xtit, $
           ytitle=ytit, yrange=yrng, thick=4, xtickn=xspaces, $
           xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata

     ;; Shade the redshift uncertainty
     ;z_pm = x_relvel(gd_strct[qq].z_fg, [1,-1]*gd_strct[qq].z_fsig)
     ;x_curvefill, 1215.6701*(1+z_pm), replicate(yrng[0],2), replicate(yrng[1],2), $
     ;             color=clr.lightgray

     ;; Data
     oplot, wave_bg, flux_bg, color=lclr, psym=10, thick=2
     oplot, wave_bg, sig_bg, color=clr.green, psym=10, thick=2, linesty=1

     ;; Plot Continuum
     oplot, data.wave, data.conti*SCL, color=clr.cyan, linesty=2, thick=4

     ;; Mark positions
     plotsym, 1, 2.0, thick=5
     ;wvv =  (1215.6701 + [-10, 10])*(1+gd_strct[qq].z_fg)
     ;gdp = where(wave_bg GT wvv[0] and wave_bg LT wvv[1], ngdp)
     ;srt = sort(flux_bg[gdp])
     ;yarr = flux_bg[gdp[srt[(round(0.95*ngdp)) < (ngdp-1)]]] + 0.1 * yrng[1]
     oplot, [1215.6701*(1+gd_strct[qq].z_fg)],  [apos[cnt]], color=clr.yellow, psym=8

     ;; Label
     x_radec, ras, decs, gd_strct[qq].rad, gd_strct[qq].decd, /flip
     rav = strsplit(ras, ':', /extract)
     decv = strsplit(decs, ':', /extract)
     qlbl = 'J'+rav[0]+rav[1]+rav[2]+decv[0]+decv[1]+decv[2]

     if lflg[cnt] EQ 0 then xlbl = 0.05 else xlbl = 0.6
     xyouts, xrng[0]+(xrng[1]-xrng[0])*xlbl, yrng[1]-(yrng[1]-yrng[0])*0.18, $
             qlbl, color=lclr, charsiz=lsz, align=0.0
     lyafil = strtrim(getenv('QSO_DIR')+gd_strct[qq].lya_fil, 2)
     spec = qpq6_get_instrument(lyafil)
     xyouts, xrng[0]+(xrng[1]-xrng[0])*xlbl, yrng[1]-(yrng[1]-yrng[0])*0.32, $
             spec, color=lclr, charsiz=lsz, align=0.0

     ;; Axes
;     for ii=0,1 do begin
;        if ii EQ 1 then spc = replicate(' ', 30) else spc = ''
;        axis, xaxis=ii, color=clr.black, charsize=csz, $
;              xmargin=xmrg, ymargin=ymrg, $
;              xthick=4, ythick=4, $
;              xrange=xrng, xstyle=1, $
;              xtickn=spc
;     endfor
  
     ;; Increment
     cnt++
  endfor

  xyouts, 0.05, 0.55, 'Normalized Flux', $
          alignment=0.5, ORIENTATION=90., /normal, charsize=lsz, color=lclr

  ;; Close Ps
  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]

  print, 'fig_exmpl_spec:  All done!'
       
  return
end
      
      
