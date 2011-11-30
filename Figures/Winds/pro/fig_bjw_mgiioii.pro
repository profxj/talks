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
pro fig_bjw_mgiioii, datfil, NTOT=NTOT, CSIZE=csize, PSFILE=psfile, XTINT=xtint, $
              FINE=fine, INFIL=infil, LSIZE=lsize, YPS=yps, LTHICK=lthick, $
              VSPLT=vsplt


  datfil = '~/Dropbox/Keck/DEIMOS/Ben_MgII/spec.coadd.all_mgii_inorm1.fits'

;  
  if not keyword_set(PSFILE) then psfile = 'fig_bjw_mgiioii.ps'
  if not keyword_set( NTOT ) then ntot = 8L
  if not keyword_set( LSIZE ) then lsize = 1.5
  if not keyword_set( LTHICK ) then lthick = 6.
;  if not keyword_set( INFIL ) then stop
  if not keyword_set( DATFIL ) then stop
  if not keyword_set( BLSZ ) then blsz = 1.7
  if not keyword_set( CSIZE ) then csize = 2.8

;  if not keyword_set( PSFIL ) then psfil = 'tmp.ps'

;  Read instrument file list
;readcol, instr_list, instr_fil, inst_dw, inst_w0, format='a,f,f'



  ;; Open vel_fil
  close, /all
;  openr, 11, infil
  zabs = 0.
;  readf, 11, zabs, FORMAT='(f12.7)'
;  readf, 11, vmin,vmax
;  vmnx = [vmin,vmax]
  vmnx = [-1700., 500]
  nlin = 0
;  readf, 11, nlin
  nlin = 2
  new_nlin = nlin

  npg = nlin/ntot + ((nlin MOD ntot) GT 0)

  npx = 1 
  npy = 3

  if keyword_set( PSFILE ) then x_psopen, psfile, /portrait


  ;;starting at p.multi[0], plot p.multi[1] rows and p.multi[2] columns
  ;;with p.multi[3] z dimensions and going from top to bottom (p.multi[4]
  !p.multi=[npx*(npy+1),npx,npy+1,0,1]
  clr = getcolor(/load)
  lclr = clr.white
;  lclr = clr.black

  ;; Vel strct
  wrest = 0.d
  svinst = -1 ;;because SiC 1B is binary 1, which makes instr list # = 0
  nblnd = 0
  flg = 0
      
  spec = xmrdfits(datfil,1)
  fx = spec.spec
  wave = spec.lambda
  npix = n_elements(fx)

  nblnd = 0
      
  ;; LOOP
  for ii=0L,new_nlin-1 do begin
     if (ii mod npy) EQ 0 then begin
        if ii GT 0 then offy = (new_nlin MOD 2) else offy = 0
        !p.multi=[npx*(npy+1) - (npy+1)*(ii/npy) - offy, $
                  npx,npy+1,0,1]
     endif
;      readf, 11, wrest, ymin, ymax, nblnd
     if ii EQ 0 then begin
        wrest = 3729.86 
        ymnx = [0., 550]
        flg = 0
        nam = '[OII]'
     endif else begin 
        wrest = 2803.531
        ymnx = [0., 45.]
        nam = 'MgII'
     endelse
          
     
     ;; Set vel array
     x_pixminmax, wave, wrest, zabs, vmnx[0], vmnx[1], PIXMIN=pmn, $
                  PIXMAX=pmx, VELO=velo
     
     ;; Y-axis
;     if ymnx[1]-ymnx[0] GT 0.8 then ytint = 0.5 else begin
;        if ymnx[1]-ymnx[0] GT 0.4 then ytint = 0.2 else ytint = 0.1
;     endelse
          
     ;; Plot
     ysty=1
     if (ii MOD npy) NE (npy-1) and (ii NE new_nlin-1) then begin
        xspaces = replicate(' ',30) 
     endif else begin 
        if keyword_set(XSPACES) then delvarx, xspaces
     endelse
     if (ii GT npy-1) then yspaces = replicate(' ',30) 
     
     plot, velo[pmn:pmx], fx[pmn:pmx], xrange=vmnx, $
           yrange=ymnx, xtickn=xspaces, xmargin=[6,1.5], $
           ymargin=[0,0], NODATA=nblnd, $ ;ytickn=yspaces, $
           charsize=csize, psym=10, background=clr.white, color=lclr, $
           xstyle=1, ystyle=ysty, thick=lthick, ytickinterval=ytint
     
     ;; Lines
     oplot, [0., 0.], ymnx, color=clr.cyan, thick=4
;     oplot, [-10000., 10000.], [0.,0.], color=clr.red, linestyle=3,thick=3
;     oplot, [-10000., 10000.], [1.,1.], color=clr.green, linestyle=1,thick=3

     if ii EQ 0 then begin
        dv = (3727.1-3729.86)/3728. * 3e5
        xyouts, -1550., 400., 'Weiner+09', color=clr.orange, charsi=2.0, align=0.
     endif else begin
        dv = (2796.352-2803.531)/2800. * 3e5
     endelse
     oplot, [dv, dv], ymnx, color=clr.cyan, thick=4
     
     ;; BLENDS
     nlow = pmn
     for jj=0L,nblnd-1 do begin
;        readf, 11, vmin, vmax, FORMAT='(f,f)'
        mn = min(abs(velo-vmin),pixmin)
        mn = min(abs(velo-vmax),pixmax)
        
        ;; Plot good
        nlow = nlow < pixmin
        oplot, velo[nlow:pixmin], fx[nlow:pixmin], $
               color=lclr, psym=10, thick=lthick
        ;; Plot blend
        oplot, velo[pixmin:pixmax], fx[pixmin:pixmax], color=clr.orange, $
               psym=10, linestyle=2, thick=1.
        ;; End
        if jj EQ nblnd-1 then begin
           if pixmax LT pmx then $
              oplot, velo[pixmax:pmx], fx[pixmax:pmx], color=lclr, $
                     psym=10, thick=lthick
        endif
        ;; Set nlow
        nlow = pixmax
     endfor
        
     ;; Labels
     case flg of 
        0: xyouts, 0.05*(vmnx[1]-vmnx[0])+vmnx[0], $
                   ymnx[0]+ (ymnx[1]-ymnx[0])*0.20, $
                   strtrim(nam,2), color=lclr, charsize=LSIZE
        1: xyouts, 0.5*(vmnx[1]-vmnx[0])+vmnx[0], $
                   ymnx[0]+ (ymnx[1]-ymnx[0])*0.68, $
                   strtrim(nam,2), color=clr.black, charsize=LSIZE
        2: xyouts, 0.70*(vmnx[1]-vmnx[0])+vmnx[0], $
                   ymnx[0]+ (ymnx[1]-ymnx[0])*0.10, $
                   strtrim(nam,2), color=clr.black, charsize=LSIZE
        3: xyouts, 0.50*(vmnx[1]-vmnx[0])+vmnx[0], $
                   ymnx[0]+ (ymnx[1]-ymnx[0])*0.10, $
                   strtrim(nam,2), color=clr.black, charsize=LSIZE
        4: xyouts, 0.40*(vmnx[1]-vmnx[0])+vmnx[0], $
                   ymnx[0]+ (ymnx[1]-ymnx[0])*0.10, $
                   strtrim(nam,2), color=clr.black, charsize=LSIZE
        5: xyouts, 0.20*(vmnx[1]-vmnx[0])+vmnx[0], $
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
      
  endfor 
  xyouts, 0.03, 0.75, 'Relative Flux', color=lclr, $
          alignment=0.5, ORIENTATION=90., /normal, charsize=blsz
  yps = 0.425
  xyouts, 0.57, yps, 'Relative Velocity  (km s!u-1!N)', $
          alignment=0.5, /normal, charsize=blsz, color=lclr

  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]
  
  return
end
