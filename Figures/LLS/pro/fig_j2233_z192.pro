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
pro fig_j2233_z192, datfil, NTOT=NTOT, CSIZE=csize, PSFILE=psfile, XTINT=xtint, $
               FINE=fine, INFIL=infil, LSIZE=lsize, YPS=yps

; lowzovi_prsdat -- Reads in DLA data to a structure

;if (N_params() LT 0) then begin 
;    print,'Syntax - ' + $
;      'fuse_velplt, strct, instr_list, vel_fil, NTOT=, CSIZE=, ' + $
;      'LSIZE=, PSFILE=, XTINT= (v1.1)' 
;    return
;endif 
;  
  if not keyword_set(PSFILE) then psfile = 'fig_j2233_z192.ps'
  if not keyword_set( NTOT ) then ntot = 8L
  if not keyword_set( LSIZE ) then lsize = 1.5
  if not keyword_set( LTHICK ) then lthick = 3.
  if not keyword_set( INFIL ) then infil = 'Input/fig_j2233_z192.inp'
  if not keyword_set( DATFIL ) then datfil = '/u/xavier/HST/HDFS/f2233_f.fits'
;  if not keyword_set( PSFIL ) then psfil = 'tmp.ps'

;  Read instrument file list
;readcol, instr_list, instr_fil, inst_dw, inst_w0, format='a,f,f'



  ;; Open vel_fil
  close, /all
  openr, 11, infil
  readf, 11, zabs, FORMAT='(f12.7)'
  readf, 11, vmin,vmax
  vmnx = [vmin,vmax]
  nlin = 0
  readf, 11, nlin

; PLOT
;if nlin LE ntot then begin
  if nlin GT (ntot/2) then begin
      npx = 2 
      npy = nlin/2
  endif else begin
      npx = 1 
      npy = nlin
  endelse

  if not keyword_set( CSIZE ) then csize = 2.5

  ;; PSFILE
  if keyword_set( PSFILE ) then x_psopen, psfile, /maxs


  ;;starting at p.multi[0], plot p.multi[1] rows and p.multi[2] columns
  ;;with p.multi[3] z dimensions and going from top to bottom (p.multi[4]
  !p.multi=[npx*(npy+1),npx,npy+1,0,1]
  clr = getcolor(/load)

  ;; Vel strct
  wrest = 0.d
  svinst = -1  ;;because SiC 1B is binary 1, which makes instr list # = 0
  nblnd = 0
  flg = 0

  fx = x_readspec(datfil, wav=wave, NPIX=npix) 

  ;; LOOP
  for ii=0L,nlin-1 do begin
      if (ii mod npy) EQ 0 then begin
          !p.multi=[npx*(npy+1) - (npy+1)*(ii/npy),npx,npy+1,0,1]
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
      if ymnx[1]-ymnx[0] GT 0.6 then ytint = 0.5 else ytint = 0.2
      
      ;; Plot
      ysty=1
      if (ii MOD npy) NE (npy-1) and (ii NE nlin-1) then begin
          xspaces = replicate(' ',30) 
      endif else begin 
          if keyword_set(XSPACES) then delvarx, xspaces
      endelse
      if (ii GT npy-1) then yspaces = replicate(' ',30) 
      
      plot, velo[pmn:pmx], fx[pmn:pmx], xrange=vmnx, $
        yrange=ymnx, xtickn=xspaces, xmargin=[7,0], $
        ymargin=[0,0], NODATA=nblnd, $;ytickn=yspaces, $
        charsize=csize, psym=10, background=clr.white, color=clr.black, $
        xstyle=1, ystyle=ysty, thick=lthick, ytickinterval=ytint

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
          ;; FeII
          if strmid(nam,0,4) EQ 'FeII' then begin
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
          1: xyouts, 0.2*(vmnx[1]-vmnx[0])+vmnx[0], $
            ymnx[0]+ (ymnx[1]-ymnx[0])*0.68, $
            strtrim(nam,2), color=clr.black, charsize=LSIZE
          2: xyouts, 0.05*(vmnx[1]-vmnx[0])+vmnx[0], $
            ymnx[0]+ (ymnx[1]-ymnx[0])*0.31, $
            strtrim(nam,2), color=clr.black, charsize=LSIZE
          else: stop
      endcase
      if keyword_set(FINE) then begin 
          if abslin.E NE 0 or strmid(nam,0,4) EQ 'FeII' then $
            xyouts, 0.05*(vmnx[1]-vmnx[0])+vmnx[0], $
                    ymnx[0]+ (ymnx[1]-ymnx[0])*0.10, $
                    strtrim(nam2,2), color=clr.black, charsize=LSIZE
      endif
      
      
      ;; Lines
      oplot, [0., 0.], ymnx, color=clr.cyan, linestyle=2
      oplot, [-10000., 10000.], [0.,0.], color=clr.pink, linestyle=3
      oplot, [-10000., 10000.], [1.,1.], color=clr.green, linestyle=3
  endfor 
  xyouts, 0.04, 0.55, 'Normalized Flux', color=clr.black, $
          alignment=0.5, ORIENTATION=90., /normal, charsize=2.0 
  if not keyword_set(YPS) then begin
      if ntot EQ 18 then yps = 0.03 else yps = 0.05
  endif
  xyouts, 0.57, yps, 'Relative Velocity (km/s)', color=clr.black, $
          alignment=0.5, /normal, charsize=2.0 

  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]
  
  return
end
