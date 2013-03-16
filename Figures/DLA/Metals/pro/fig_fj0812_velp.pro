;+ 
; NAME:
; fuse_velplt
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
pro fig_fj0812_velp, datfil, NTOT=NTOT, CSIZE=csize, PSFILE=psfile, XTINT=xtint, $
              FINE=fine, INFIL=infil, LSIZE=lsize, YPS=yps, LTHICK=lthick, $
              VSPLT=vsplt

; lowzovi_prsdat -- Reads in DLA data to a structure

;if (N_params() LT 0) then begin 
;    print,'Syntax - ' + $
;      'fuse_velplt, strct, instr_list, vel_fil, NTOT=, CSIZE=, ' + $
;      'LSIZE=, PSFILE=, XTINT= (v1.1)' 
;    return
;endif 
;  
  if not keyword_set(PSFILE) then psfile = 'fig_fj0812_velp.ps'
  if not keyword_set( NTOT ) then ntot = 12L
  if not keyword_set( LSIZE ) then lsize = 1.4
  if not keyword_set( LTHICK ) then lthick = 6.
  if not keyword_set( INFIL ) then infil = 'Input/fig_fj0812_velp.inp'
  if not keyword_set( BLSZ ) then blsz = 1.7
  if not keyword_set( CSIZE ) then csize = 2.2

;  if not keyword_set( PSFIL ) then psfil = 'tmp.ps'

;  Read instrument file list
;readcol, instr_list, instr_fil, inst_dw, inst_w0, format='a,f,f'

  ;; Open vel_fil
  close, /all
  openr, 11, infil
  nlin = 0
  readf, 11, nlin

  datfil = '~/Keck/HIRES/RedData/FJ0812+32/FJ0812+32_f.fits'
  fx = x_readspec(datfil, sig=sig, wav=wave, /auto)
  vmnx = [-80, 45]
  zabs = 2.62633d

  npg = nlin/ntot + ((nlin MOD ntot) GT 0)

  for ff=0L,npg-1 do begin

      ;; Number of pages
      new_nlin = (nlin - ff*ntot) < ntot
      if new_nlin GT (ntot/2) then begin
          npx = 2 
          npy = new_nlin/2 + (new_nlin MOD 2)
      endif else begin
          npx = 1 
          npy = new_nlin
      endelse
      ;; PSFILE

      if npg EQ 1 then begin
          npsfile = psfile 
      endif else begin
          case ff of
              0: cps = 'a'
              1: cps = 'b'
              2: cps = 'c'
              3: cps = 'd'
              else: stop
          endcase
          pos = strpos(psfile,'.ps')
          npsfile = strmid(psfile,0,pos)+cps+'.ps'
      endelse

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
          spec_fil = ''
          readf, 11, wrest, flg
          readf, 11, ymin, ymax
          ymnx = [ymin, ymax]
          
          x_pixminmax, wave, wrest, zabs, vmnx[0], vmnx[1], PIXMIN=pmn, $
                       PIXMAX=pmx, VELO=velo

          readf, 11, nblnd
          
          ;; Y-axis
          if ymnx[1]-ymnx[0] GT 0.8 then ytint = 0.5 else begin
              if ymnx[1]-ymnx[0] GT 0.4 then ytint = 0.2 else ytint = 0.1
          endelse
          
          ;; Plot
          ysty=1
          if (ii MOD npy) NE (npy-1) and (ii NE new_nlin-1) then begin
              xspaces = replicate(' ',30) 
          endif else begin 
              if keyword_set(XSPACES) then delvarx, xspaces
          endelse
          if (ii GT npy-1) then yspaces = replicate(' ',30) 
          
          plot, velo, fx, xrange=vmnx, $
                yrange=ymnx, xtickn=xspaces, xmargin=[6,0], $
                ymargin=[0,0], NODATA=nblnd, $ ;ytickn=yspaces, $
                charsize=csize, psym=10, background=clr.white, color=clr.black, $
                xstyle=1, ystyle=ysty, thick=lthick, ytickinterval=ytint

          ;; Lines
          oplot, [0., 0.], ymnx, color=clr.blue, linestyle=2, thick=4
          oplot, [-10000., 10000.], [0.,0.], color=clr.red, linestyle=3,thick=4
          oplot, [-10000., 10000.], [1.,1.], color=clr.darkgreen, linestyle=1,thick=4
          
          ;; BLENDS
          nlow = pmn
          for jj=0L,nblnd-1 do begin
              readf, 11, vmin, vmax, FORMAT='(f,f)'
              mn = min(abs(velo-vmin),pixmin)
              mn = min(abs(velo-vmax),pixmax)
              
              ;; Plot good
              nlow = nlow < pixmin
              oplot, velo[nlow:pixmin], fx[nlow:pixmin], $
                     color=clr.black, psym=10, thick=lthick
              ;; Plot blend
              oplot, velo[pixmin:pixmax], fx[pixmin:pixmax], color=clr.orange, $
                     psym=10, linestyle=2, thick=1.
              ;; End
              if jj EQ nblnd-1 then begin
                  if pixmax LT pmx then $
                oplot, velo[pixmax:pmx], fx[pixmax:pmx], color=clr.black, $
                       psym=10, thick=lthick
              endif
              ;; Set nlow
              nlow = pixmax
          endfor
        
          ;; Labels
          getfnam, wrest, fv, nam 
          nam = strtrim(nam,2)
          if keyword_set(FINE) then begin
          abslin = x_setline(wrest)
          ;; SiII
          if strmid(nam,0,4) EQ 'SiII' and abslin.j EQ 0 then abslin.j = 0.5
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
          endif else begin
              if round(abslin.j) NE fix(abslin.j) then $
                nam2 = 'J='+strtrim(round(abslin.j*2),2)+'/2' $
              else nam2 = 'J='+strtrim(round(abslin.j),2)
          endelse
      endif

      case flg of 
          0: xyouts, 0.05*(vmnx[1]-vmnx[0])+vmnx[0], $
                     ymnx[0]+ (ymnx[1]-ymnx[0])*0.10, $
                     strtrim(nam,2), color=clr.black, charsize=LSIZE
          1: xyouts, 0.07*(vmnx[1]-vmnx[0])+vmnx[0], $
                     ymnx[0]+ (ymnx[1]-ymnx[0])*0.10, $
                     strtrim(nam,2), color=clr.black, charsize=LSIZE, align=0.
          2: xyouts, 0.05*(vmnx[1]-vmnx[0])+vmnx[0], $
                     ymnx[0]+ (ymnx[1]-ymnx[0])*0.75, $
                     strtrim(nam,2), color=clr.black, charsize=LSIZE
          3: xyouts, 0.72*(vmnx[1]-vmnx[0])+vmnx[0], $
                     ymnx[0]+ (ymnx[1]-ymnx[0])*0.14, $
                     strtrim(nam,2), color=clr.black, charsize=LSIZE
          4: xyouts, 0.03*(vmnx[1]-vmnx[0])+vmnx[0], $
                     ymnx[0]+ (ymnx[1]-ymnx[0])*0.12, $
                     strtrim(nam,2), color=clr.black, charsize=LSIZE
          5: xyouts, 0.40*(vmnx[1]-vmnx[0])+vmnx[0], $
                     ymnx[0]+ (ymnx[1]-ymnx[0])*0.64, $
                     strtrim(nam,2), color=clr.black, charsize=LSIZE
          6: xyouts, 0.37*(vmnx[1]-vmnx[0])+vmnx[0], $
                     ymnx[0]+ (ymnx[1]-ymnx[0])*0.12, $
                     strtrim(nam,2), color=clr.black, charsize=LSIZE
          7: xyouts, 0.2*(vmnx[1]-vmnx[0])+vmnx[0], $
                     ymnx[0]+ (ymnx[1]-ymnx[0])*0.15, $
                     strtrim(nam,2), color=clr.black, charsize=LSIZE
          8: xyouts, 0.70*(vmnx[1]-vmnx[0])+vmnx[0], $
                     ymnx[0]+ (ymnx[1]-ymnx[0])*0.75, $
                     strtrim(nam,2), color=clr.black, charsize=LSIZE
          else: stop
      endcase
      if keyword_set(FINE) then begin 
;          if abslin.E NE 0 or strmid(nam,0,4) EQ 'FeII' then begin
          case flg of 
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
  xyouts, 0.03, 0.55, 'Normalized Flux', $
          alignment=0.5, ORIENTATION=90., /normal, charsize=blsz
  if not keyword_set(YPS) then begin
      if ntot EQ 18 then yps = 0.03 else yps = 0.07
  endif
  xyouts, 0.57, yps, 'Relative Velocity  (km s!u-1!N)', $
          alignment=0.5, /normal, charsize=blsz

  if keyword_set( PSFILE ) then x_psclose
endfor
  !p.multi=[0,1,1]
  
  return
end
