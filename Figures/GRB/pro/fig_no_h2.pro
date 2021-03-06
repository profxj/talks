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
pro fig_no_h2, datfil, NTOT=NTOT, CSIZE=csize, PSFILE=psfile, XTINT=xtint, $
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
  if not keyword_set(PSFILE) then psfile = 'fig_no_h2.ps'
  if not keyword_set(INFIL) then infil = 'Input/fig_no_h2.inp'
  if not keyword_set( NTOT ) then ntot = 3L
  if not keyword_set( LSIZE ) then lsize = 1.5
  if not keyword_set( LTHICK ) then lthick = 7.
  if not keyword_set( BLSZ ) then blsz = 1.7
  if not keyword_set( CSIZE ) then csize = 3.1

;  if not keyword_set( PSFIL ) then psfil = 'tmp.ps'

;  Read instrument file list
;readcol, instr_list, instr_fil, inst_dw, inst_w0, format='a,f,f'


  h2lin = fuse_h2lin()

  npx = 2
  npy = 6
  nlin = 3
  if keyword_set( PSFILE ) then x_psopen, psfile, /maxs
  !p.multi=[npx*(npy+2),npx,npy+2,0,1]
  clr = getcolor(/load)

  ;; Open vel_fil
  close, /all
  openr, 11, infil
  datflg = 0
  datfil = ''
  flbl = ['(a)','(b)','(c)','(d)']

  for qq=0L,3 do begin

; PLOT
;if nlin LE ntot then begin

      readf, 11, datfil
      readf, 11, datflg
      fx = x_readspec(datfil, wav=wave, NPIX=npix, INFLG=datflg) 

;      if qq EQ 0 then airtovac, wave

      readf, 11, zabs, FORMAT='(f12.7)'
      readf, 11, vmin,vmax
      vmnx = [vmin,vmax]
      
      !p.multi=[npx*(npy+2) - (nlin+1)*(qq), npx,npy+2,0,1]

      ;;starting at p.multi[0], plot p.multi[1] rows and p.multi[2] columns
      ;;with p.multi[3] z dimensions and going from top to bottom (p.multi[4]

      ;; Vel strct
      wrest = 0.d
      svinst = -1  ;;because SiC 1B is binary 1, which makes instr list # = 0
      nblnd = 0
      flg = 0
      
      
      ;; LOOP
      for ii=0L,nlin-1 do begin
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
          if (ii MOD npy) NE (npy-1) and (ii NE nlin-1) then begin
              xspaces = replicate(' ',30) 
              if keyword_set(xtit) then delvarx, xtit
          endif else begin 
              if keyword_set(XSPACES) then delvarx, xspaces
              xtit = 'Relative Velocity (km s!u-1!N)'
          endelse
          if (ii GT npy-1) then yspaces = replicate(' ',30) 
          
          plot, velo[pmn:pmx], fx[pmn:pmx], xrange=vmnx, $
                yrange=ymnx, xtickn=xspaces, xmargin=[6,0], $
                ymargin=[0,0], NODATA=nblnd, $ ;ytickn=yspaces, $
                charsize=csize, psym=10, background=clr.black, color=clr.lightgray, $
                xstyle=1, ystyle=ysty, thick=lthick, ytickinterval=ytint, $
                xtitle=xtit 

          ;; Lines
;          oplot, [0., 0.], ymnx, color=clr.blue, linestyle=2, thick=1
          oplot, [-10000., 10000.], [0.,0.], color=clr.cyan, linestyle=3,thick=1
          oplot, [-10000., 10000.], [1.,1.], color=clr.green, linestyle=1,thick=1
          
          ;; BLENDS
          nlow = pmn
          for jj=0L,nblnd-1 do begin
              readf, 11, vmin, vmax, FORMAT='(f,f)'
              mn = min(abs(velo-vmin),pixmin)
              mn = min(abs(velo-vmax),pixmax)
              
              ;; Plot good
              nlow = nlow < pixmin
              oplot, velo[nlow:pixmin], fx[nlow:pixmin], $
                     color=clr.lightgray, psym=10, thick=lthick
              ;; Plot blend
              oplot, velo[pixmin:pixmax], fx[pixmin:pixmax], color=clr.orange, $
                     psym=10, linestyle=2, thick=1.
              ;; End
              if jj EQ nblnd-1 then begin
                  if pixmax LT pmx then $
                oplot, velo[pixmax:pmx], fx[pixmax:pmx], color=clr.lightgray, $
                       psym=10, thick=lthick
              endif
              ;; Set nlow
              nlow = pixmax
          endfor
        
          ;; Labels
          if ii EQ 0 then begin
              ;; H2
              mt = where(abs(h2lin.wrest-wrest) LT 0.001,nmt)
              if nmt NE 1 then stop
              nam = h2lin[mt].label
          endif else getfnam, wrest, fv, nam 

          nam = strtrim(nam,2)

          case flg of 
              0: xyouts, 0.05*(vmnx[1]-vmnx[0])+vmnx[0], $
                         ymnx[0]+ (ymnx[1]-ymnx[0])*0.10, $
                         strtrim(nam,2), color=clr.lightgray, charsize=LSIZE
              1: xyouts, 0.2*(vmnx[1]-vmnx[0])+vmnx[0], $
                         ymnx[0]+ (ymnx[1]-ymnx[0])*0.68, $
                         strtrim(nam,2), color=clr.lightgray, charsize=LSIZE
              2: xyouts, 0.05*(vmnx[1]-vmnx[0])+vmnx[0], $
                         ymnx[0]+ (ymnx[1]-ymnx[0])*0.31, $
                         strtrim(nam,2), color=clr.lightgray, charsize=LSIZE
              3: xyouts, 0.70*(vmnx[1]-vmnx[0])+vmnx[0], $
                         ymnx[0]+ (ymnx[1]-ymnx[0])*0.31, $
                         strtrim(nam,2), color=clr.lightgray, charsize=LSIZE
              4: xyouts, 0.70*(vmnx[1]-vmnx[0])+vmnx[0], $
                         ymnx[0]+ (ymnx[1]-ymnx[0])*0.14, $
                         strtrim(nam,2), color=clr.lightgray, charsize=LSIZE
              5: xyouts, 0.40*(vmnx[1]-vmnx[0])+vmnx[0], $
                         ymnx[0]+ (ymnx[1]-ymnx[0])*0.64, $
                         strtrim(nam,2), color=clr.lightgray, charsize=LSIZE
              6: xyouts, 0.37*(vmnx[1]-vmnx[0])+vmnx[0], $
                         ymnx[0]+ (ymnx[1]-ymnx[0])*0.12, $
                         strtrim(nam,2), color=clr.lightgray, charsize=LSIZE
              7: xyouts, 0.2*(vmnx[1]-vmnx[0])+vmnx[0], $
                         ymnx[0]+ (ymnx[1]-ymnx[0])*0.15, $
                         strtrim(nam,2), color=clr.lightgray, charsize=LSIZE
              8: xyouts, 0.70*(vmnx[1]-vmnx[0])+vmnx[0], $
                         ymnx[0]+ (ymnx[1]-ymnx[0])*0.75, $
                         strtrim(nam,2), color=clr.lightgray, charsize=LSIZE
              else: stop
          endcase

          if ii EQ 2 then $
            xyouts, vmnx[0]+0.9*(vmnx[1]-vmnx[0]), $
                    ymnx[0]+0.15*(ymnx[1]-ymnx[0]), flbl[qq], color=clr.lightgray, $
                    charsiz=lsize
      endfor 
  endfor

  xyouts, 0.03, 0.55, 'Normalized Flux', $
          alignment=0.5, color=clr.lightgray, ORIENTATION=90., /normal, charsize=blsz
  if not keyword_set(YPS) then begin
      if ntot EQ 18 then yps = 0.03 else yps = 0.05
  endif

  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]
  
  return
end
