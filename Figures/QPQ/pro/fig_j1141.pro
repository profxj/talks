;+ 
; NAME:
; fig_j1141
;  V1.1
;
; PURPOSE:
;    Velcoity plot of FUSE transitions for a given absorption 
;    system
; CALLING SEQUENCE:
;   
;   fuse_velplt, strct_fil, instr_list, vel_fil, NTOT=, CSIZE=,
;   LSIZE=, PSFILE=, XTINT=
;
; INPUTS:
;
; RETURNS:
;  strct_fil -- FITS file for the FUSE abs lin structure
;  instr_list -- List of instrument files
;  vel_fil -- Input file for velocity plot
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  NTOT -- Number of plots per page [default: 16]
;  LSIZE -- Label size [default: 1.8]
;  CSIZE -- Numbering character size [default: 1.8]
;  XTINT -- xtick interval
;  MULTI -- overplot multiple detections
;
; OPTIONAL OUTPUTS:
;  PSFILE -- Postscript filename
;  MSGFIL -- filename to receive messages from MULTI function
;
; COMMENTS:
;
; EXAMPLES:
;   fuse_calcewn, struct, fil_instr
;
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   12-Sep-2003 Written by JXP
;    1-Apr-2005 MULTI capability added by KLC
;   28-Apr-2005 Added zabs in upper right-hand corner of MULTI plots
;-
;------------------------------------------------------------------------------
;  grb_velplt, PSFILE='Figures/grb_velplt.ps'
pro fig_j1141, datfil, NTOT=NTOT, CSIZE=csize, PSFILE=psfile, XTINT=xtint, $
              FINE=fine, INFIL=infil, LSIZE=lsize, YPS=yps, LTHICK=lthick, $
              VSPLT=vsplt, TALK=talk

; lowzovi_prsdat -- Reads in DLA data to a structure

;if (N_params() LT 0) then begin 
;    print,'Syntax - ' + $
;      'fuse_velplt, strct, instr_list, vel_fil, NTOT=, CSIZE=, ' + $
;      'LSIZE=, PSFILE=, XTINT= (v1.1)' 
;    return
;endif 
  compile_opt strictarr
;  
  ;; Read data
  bg_qso = '~/Dropbox/QSOPairs/data/lya_spec_all/LRIS_GMOS_PRELIM/J1141+0724A.fits'
  fx= x_readspec(bg_qso, wav=wave, NPIX=npix, inflg=2) 
  bg_qso_conti = '~/Dropbox/QSOPairs/data/lya_spec_all/LRIS_GMOS_PRELIM/J1141+0724A_c.fits'
  conti= xmrdfits(bg_qso_conti)
  fx = fx / conti

  if not keyword_set(PSFILE) then psfile = 'fig_j1141.ps'
  if not keyword_set( NTOT ) then ntot = 16L
  if not keyword_set( LSIZE ) then lsize = 1.2
  if not keyword_set( LTHICK ) then lthick = 4.
  if not keyword_set( PTHICK ) then pthick = 2.
  if not keyword_set( INFIL ) then infil = 'Input/fig_j1141.inp'
  if not keyword_set( BLSZ ) then blsz = 1.7
  if not keyword_set( CSIZE ) then csize = 2.5

  ;FINE = 1

  c=x_constants()

  ;; Open vel_fil
  close, /all
  openr, 11, infil
  readf, 11, zabs, FORMAT='(f12.7)'
  readf, 11, vmin,vmax
  vmnx = [vmin,vmax]
  nlin = 0
  readf, 11, nlin

  ;; Number of pages
  new_nlin = nlin 
  npx = 2 
  npy = new_nlin/npx

  npsfile = psfile 
  if keyword_set( PSFILE ) then x_psopen, npsfile, /portrait

  ;;starting at p.multi[0], plot p.multi[1] rows and p.multi[2] columns
  ;;with p.multi[3] z dimensions and going from top to bottom (p.multi[4]
  !p.multi=[npx*(npy+1),npx,npy+1,0,1]
  clr = getcolor(/load)
  
  ;; Vel strct
  wrest = 0.d
  svinst = -1  ;;because SiC 1B is binary 1, which makes instr list # = 0
  nblnd = 0
  flg = 0
  
      
  ;; LOOP
  for ii=0L,new_nlin-1 do begin
     if (ii mod npy) EQ 0 then begin
        if ii GT 0 then offy = (new_nlin MOD 2) else offy = 0
        !p.multi=[npx*(npy+1) - (npy+1)*(ii/npy) - offy, $
                  npx,npy+1,0,1]
     endif
;      readf, 11, wrest, ymin, ymax, nblnd
          
     readf, 11, wrest, flg
     readf, 11, ymin, ymax
     ymnx = [ymin, ymax]
     
     readf, 11, nblnd
          
     ;; Set vel array
     x_pixminmax, wave, wrest, zabs, vmnx[0], vmnx[1], PIXMIN=pmn, $
                  PIXMAX=pmx, VELO=velo
     
          
     ;; Y-axis
     if ymnx[1]-ymnx[0] GT 0.8 then ytint = 0.5 else begin
        if ymnx[1]-ymnx[0] GT 0.4 then ytint = 0.2 else ytint = 0.1
     endelse
     
     ;; Plot
     ysty=1
     if (ii MOD npy) NE (npy-1) and (ii NE new_nlin-1) then begin
        xspaces = replicate(' ',30) 
        xtitle = ''
     endif else begin 
        if keyword_set(XSPACES) then delvarx, xspaces
        xtitle='Relative Velocity (km s!u-1!N)'
     endelse
     if (ii GT npy-1) then yspaces = replicate(' ',30) 
     
     plot, velo[pmn:pmx], fx[pmn:pmx], xrange=vmnx, $
           yrange=ymnx, xtickn=xspaces, xmargin=[6,0], $
           ymargin=[0,0], NODATA=nblnd, $ ;ytickn=yspaces, $
           charsize=csize, psym=10, background=clr.white, color=clr.black, $
           xstyle=1, ystyle=ysty, thick=pthick, ytickinterval=ytint,$
           xtitle=xtitle
     
     ;; Lines
     oplot, [0., 0.], ymnx, color=clr.blue, linestyle=2, thick=1
     oplot, [-10000., 10000.], [0.,0.], color=clr.darkgreen, linestyle=3,thick=1
     oplot, [-10000., 10000.], [1.,1.], color=clr.darkgreen, $
            linestyle=1,thick=2
     
     ;; BLENDS
     nlow = pmn
     if keyword_set(TALK) then bclr = clr.red else bclr = clr.orange
     for jj=0L,nblnd-1 do begin
        readf, 11, vmin, vmax, FORMAT='(f,f)'
        mn = min(abs(velo-vmin),pixmin)
        mn = min(abs(velo-vmax),pixmax)
        
        ;; Plot good
        nlow = nlow < pixmin
        oplot, velo[nlow:pixmin], fx[nlow:pixmin], $
               color=clr.black, psym=10, thick=pthick
        ;; Plot blend
        oplot, velo[pixmin:pixmax], fx[pixmin:pixmax], $
               color=bclr, $
               psym=10, linestyle=2, thick=1.
        ;; End
        if jj EQ nblnd-1 then begin
           if pixmax LT pmx then $
              oplot, velo[pixmax:pmx], fx[pixmax:pmx], color=clr.black, $
                     psym=10, thick=pthick
        endif
        ;; Set nlow
        nlow = pixmax
     endfor

     ;; Labels
     getfnam, wrest, fv, nam 
     nam = strtrim(nam,2)
     

     nam2 = ''
     if keyword_set(FINE) then begin
        abslin = x_setline(wrest, /override)
        if size(abslin,/type) NE 2 then begin
           ;; SiII
           if strmid(nam,0,4) EQ 'SiII' and abslin.j EQ 0 then begin
              abslin.j = 0.5
              nam2 = 'J=1/2'
           endif
           if strmid(nam,0,4) EQ 'SiII' and abslin.j EQ 1.5 then begin
              abslin.j = 0.5
              nam = 'SiII'+strmid(nam,5)
              nam2 = 'J=3/2'
           endif
           ;; CII
           if strmid(nam,0,3) EQ 'CII' and abslin.j EQ 0 then abslin.j = 0.5
           ;; OI
           if strmid(nam,0,2) EQ 'OI' and abslin.j EQ 0 $
              and abs(wrest-1306.0286) GT 0.1 then abslin.j = 2.0
           
           ;; FeII
           if strmid(nam,0,4) EQ 'FeII' then begin
              if strmid(nam,4,1) EQ '*' then nam = 'FeII'+strmid(nam,5)
              if abslin.j EQ 0 then abslin.j = 4.5
              nam2 = 'J='+strtrim(round(abslin.j*2),2)+'/2' 
           endif 
        endif else nam2 = ''
     endif
     
     case flg of 
        0: xyouts, 0.05*(vmnx[1]-vmnx[0])+vmnx[0], $
                   ymnx[0]+ (ymnx[1]-ymnx[0])*0.10, $
                   strtrim(nam,2), color=clr.black, charsize=LSIZE
        1: xyouts, 0.2*(vmnx[1]-vmnx[0])+vmnx[0], $
                   ymnx[0]+ (ymnx[1]-ymnx[0])*0.68, $
                   strtrim(nam,2), color=clr.black, charsize=LSIZE
        2: xyouts, 0.05*(vmnx[1]-vmnx[0])+vmnx[0], $
                   ymnx[0]+ (ymnx[1]-ymnx[0])*0.21, $
                   strtrim(nam,2), color=clr.black, charsize=LSIZE
        3: xyouts, 0.70*(vmnx[1]-vmnx[0])+vmnx[0], $
                   ymnx[0]+ (ymnx[1]-ymnx[0])*0.31, $
                   strtrim(nam,2), color=clr.black, charsize=LSIZE
        4: xyouts, 0.70*(vmnx[1]-vmnx[0])+vmnx[0], $
                   ymnx[0]+ (ymnx[1]-ymnx[0])*0.14, $
                   strtrim(nam,2), color=clr.black, charsize=LSIZE
        5: xyouts, 0.41*(vmnx[1]-vmnx[0])+vmnx[0], $
                   ymnx[0]+ (ymnx[1]-ymnx[0])*0.90, $
                   strtrim(nam,2), color=clr.black, charsize=LSIZE
        6: xyouts, 0.37*(vmnx[1]-vmnx[0])+vmnx[0], $
                   ymnx[0]+ (ymnx[1]-ymnx[0])*0.70, $
                   strtrim(nam,2), color=clr.black, charsize=LSIZE
        7: xyouts, 0.58*(vmnx[1]-vmnx[0])+vmnx[0], $
                   ymnx[0]+ (ymnx[1]-ymnx[0])*0.08, $
                   strtrim(nam,2), color=clr.black, charsize=LSIZE
        8: xyouts, 0.70*(vmnx[1]-vmnx[0])+vmnx[0], $
                   ymnx[0]+ (ymnx[1]-ymnx[0])*0.75, $
                   strtrim(nam,2), color=clr.black, charsize=LSIZE
        else: stop
     endcase
     if keyword_set(FINE) then begin 
;          if abslin.E NE 0 or strmid(nam,0,4) EQ 'FeII' then begin
        case flg of 
           0: xyouts, 0.2*(vmnx[1]-vmnx[0])+vmnx[0], $
                      ymnx[0]+ (ymnx[1]-ymnx[0])*0.52, $
                      strtrim(nam2,2), color=clr.black, charsize=LSIZE
           1: xyouts, 0.2*(vmnx[1]-vmnx[0])+vmnx[0], $
                      ymnx[0]+ (ymnx[1]-ymnx[0])*0.52, $
                      strtrim(nam2,2), color=clr.black, charsize=LSIZE
           2: xyouts, 0.05*(vmnx[1]-vmnx[0])+vmnx[0], $
                      ymnx[0]+ (ymnx[1]-ymnx[0])*0.10, $
                      strtrim(nam2,2), color=clr.black, charsize=LSIZE
           3: xyouts, 0.70*(vmnx[1]-vmnx[0])+vmnx[0], $
                      ymnx[0]+ (ymnx[1]-ymnx[0])*0.14, $
                      strtrim(nam2,2), color=clr.black, charsize=LSIZE
           else: stop
        endcase
     endif
  endfor 

  ;; Label
  xyouts, 0.03, 0.55, 'Normalized Flux', $
          alignment=0.5, ORIENTATION=90., /normal, charsize=blsz
  if not keyword_set(YPS) then begin
      if ntot EQ 18 then yps = 0.03 else yps = 0.09
  endif
;  xyouts, 0.57, yps, 'Relative Velocity  (km s!u-1!N)', $
;          alignment=0.5, /normal, charsize=blsz

  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]
  
  return
end
