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
pro qpq5_boot_stack, NOPLOT=noplot, EW_VALS=ew_vals

  compile_opt strictarr
;  
  if not keyword_set(PSFILE) then psfile = 'qpq5_boot_stack.ps'
  if not keyword_set( NTOT ) then ntot = 12L
  if not keyword_set( LSIZE ) then lsize = 1.1
  if not keyword_set( LTHICK ) then lthick = 4.
  if not keyword_set( PTHICK ) then pthick = 2.
  if not keyword_set( BLSZ ) then blsz = 1.7
  if not keyword_set( CSIZE ) then csize = 2.5

  if not keyword_set( CONTI ) then conti = 0.8
  if not keyword_set( NBOOT ) then nboot = 100L
  if not keyword_set(VMNX) then vmnx = [-1,1]*3500. 
  if not keyword_set(YMNX) then ymnx = [0., 1.10]
  ;; Range for EW calculation
  ew_vmnx = 2000.* [-1, 1]
  
  ;; Compile
  cd, getenv('QSO_DIR')+'tex/QPQ6/Analysis/N_HI/pro/', curr=curr
  resolve_routine, 'qpq6_lyaew_from_stack', /compile_full_file, /is_function
  cd, curr

  ;;  Stack files
  stack_root = getenv('QSO_DIR')+'/tex/QPQ6/Analysis/Stacks/'
  boot_stack_fil = ['BOOT_NOESI_stack_QPQ6_Lya_0_300.fits', $
               'BOOT_NOESI_stack_QPQ6_Lya_0_100.fits', $
               'BOOT_NOESI_stack_QPQ6_Lya_100_200.fits', $
               'BOOT_NOESI_stack_QPQ6_Lya_200_300.fits']
  stack_fil = ['NOESI_stack_QPQ6_Lya_0_300.fits', $
               'NOESI_stack_QPQ6_Lya_0_100.fits', $
               'NOESI_stack_QPQ6_Lya_100_200.fits', $
               'NOESI_stack_QPQ6_Lya_200_300.fits']
  nfil = n_elements(stack_fil)

  npx = 1
  npy = nfil

  EW_vals = fltarr(3, nfil)

  if keyword_set( PSFILE ) and not keyword_set(NOPLOT) then begin 
     x_psopen, psfile, /portrait
     !p.multi=[npx*(npy+1),npx,npy+1,0,1]
     clr = getcolor(/load)
  endif
  ;lclr = clr.lightgray
  lclr = clr.white

  ;; LOOP
  for ii=0L,nfil-1 do begin
     if (ii mod npy) EQ 0 then begin
        if ii GT 0 then offy = (nfil MOD 2) else offy = 0
        if not keyword_set(NOPLOT) then $
           !p.multi=[npx*(npy+1) - (npy+1)*(ii/npy) - offy, $
                  npx,npy+1,0,1]
     endif

          
     ;; Y-axis
     if ymnx[1]-ymnx[0] GT 0.8 then ytint = 0.5 else begin
        if ymnx[1]-ymnx[0] GT 0.4 then ytint = 0.2 else ytint = 0.1
     endelse
          
     ;; Plot
     ysty=1
     if (ii MOD npy) NE (npy-1) and (ii NE nfil-1) then begin
        xspaces = replicate(' ',30) 
        xtitle = ''
     endif else begin 
        if keyword_set(XSPACES) then delvarx, xspaces
        xtitle='Relative Velocity (km s!u-1!N)'
     endelse
     if (ii GT npy-1) then yspaces = replicate(' ',30) 
          
     if not keyword_set(NOPLOT) then begin
        plot, [0], [0], $
              xrange=vmnx, $
              yrange=ymnx, xtickn=xspaces, xmargin=[7,20], $
              ymargin=[0,0], NODATA=nblnd, $ ;ytickn=yspaces, $
              charsize=csize, psym=10, background=clr.white, color=lclr, $
              xstyle=1, ystyle=ysty, thick=pthick, ytickinterval=ytint,$
              xtitle=xtitle, xtickint=1000.

        ;; Lines
        oplot, [0., 0.], ymnx, color=clr.cyan, linestyle=2, thick=4
                                ;oplot, [-10000., 10000.], [0.,0.], color=clr.red, linestyle=3,thick=1
        oplot, [-10000., 10000.], [1.,1.], color=clr.green,  linestyle=2,thick=4
                                ;oplot, [-10000., 10000.], replicate(CONTI, 2), color=clr.green,  linestyle=2,thick=3
     endif
        
     ;; Plot Bootstrap
     velb = xmrdfits( stack_root + boot_stack_fil[ii], 0, /silen)
     npixb =n_elements(velb)
     strct = xmrdfits( stack_root + boot_stack_fil[ii], 1, /silen)
     sv_boot = fltarr(npixb, nboot)
     for jj=0L,nboot-1 do begin
        fx = xmrdfits( stack_root + boot_stack_fil[ii], 2+jj, /silen)
        if not keyword_set(NOPLOT) then oplot, velb, fx, color=clr.darkgray, thick=2
        sv_boot[*,jj] = fx
     endfor

     ;; Stats on Boot
     rms = fltarr(npixb)
     for kk=0L,npixb-1 do begin
        djs_iterstat, sv_boot[kk,*], sigrej=3., sigma=trms
        rms[kk] = trms
     endfor
     if not keyword_set(NOPLOT) then oplot, velb, rms, color=clr.red, thick=2
          
     ;; Read standard stack
     fx = xmrdfits( stack_root + stack_fil[ii], 0+keyword_set(MEDIAN), /silen)
     vel = xmrdfits( stack_root + stack_fil[ii], 3, /silen)
     npix = n_elements(vel)
     if npix NE npixb then stop
     strct = xmrdfits( stack_root + stack_fil[ii], 4, /silen)

     if not keyword_set(NOPLOT) then oplot, vel, fx, color=lclr, psym=10, thick=5

     ;; Label
      ;xyouts, 0.05*(vmnx[1]-vmnx[0])+vmnx[0], $
      ;        ymnx[0]+ (ymnx[1]-ymnx[0])*0.10, $
      ;        stack_fil[ii],  color=clr.black, charsize=LSIZE

     ;; <sigma>
     if not keyword_set(NOPLOT) then $
      xyouts, 0.05*(vmnx[1]-vmnx[0])+vmnx[0], $
              ymnx[0]+ (ymnx[1]-ymnx[0])*0.12, $
              '<!9s!X!dz!N> = '+string(round(mean(strct.z_fsig)), format='(i3)')+' km/s', $
              color=lclr, charsize=LSIZE

      ;; Nspec
      xtwo = 0.65
      if not keyword_set(NOPLOT) then $
         xyouts, xtwo*(vmnx[1]-vmnx[0])+vmnx[0], $
                 ymnx[0]+ (ymnx[1]-ymnx[0])*0.30, $
                 'N!dspec!N = '+string(n_elements(strct), format='(i2)'), $
                 color=lclr, charsize=LSIZE, align=0.
      ;; R range
      prs = strsplit(stack_fil[ii], '_', /extract)
      nprs = n_elements(prs)
      r1 = prs[nprs-2]
      p2 = strpos(prs[nprs-1], '.fits')
      r2 = strmid(prs[nprs-1], 0, p2)
      
      if not keyword_set(NOPLOT) then $
         xyouts, xtwo*(vmnx[1]-vmnx[0])+vmnx[0], $
                 ymnx[0]+ (ymnx[1]-ymnx[0])*0.15, $
                 'R!dphys!N = ['+r1+','+r2+'] kpc', $
                 color=lclr, charsize=LSIZE, align=0.

      ;; EW
      EW = qpq6_lyaew_from_stack(vel, fx, ew_vmnx, CLIM=[-4000.,4000]) 
      ;; Crude error from RMS
      vpix = where(velb GT ew_vmnx[0] and velb LT ew_vmnx[1], new)
      dv = abs(median(velb-shift(velb,1))) ; km/s
      dwv = 1215.6710 * dv/3e5
      sigEW = sqrt(total(rms[vpix]^2)) * dwv

      EW_vals[*,ii] = [EW,sigEW,round(mean(strct.R))]

      if not keyword_set(NOPLOT) then $
         xyouts, 0.05*(vmnx[1]-vmnx[0])+vmnx[0], $
                 ymnx[0]+ (ymnx[1]-ymnx[0])*0.25, $
                 'W!S!u2000!N!R!dLy!9a!X!N = '+string(EW, format='(f4.2)')+$
                 ' +/- '+string(sigEW,format='(f4.2)')+' A', $
                 color=lclr, charsize=LSIZE

      ;; <Rphys>
      if not keyword_set(NOPLOT) then $
         xyouts, 0.05*(vmnx[1]-vmnx[0])+vmnx[0], $
                 ymnx[0]+ (ymnx[1]-ymnx[0])*0.40, $
                 '<R!dphys!N> = '+strtrim(round(mean(strct.R)),2)+' kpc', $
                 color=lclr, charsize=LSIZE
   endfor 

  if not keyword_set(NOPLOT) then $
     xyouts, 0.03, 0.55, 'Normalized Flux', color=lclr, $
             alignment=0.5, ORIENTATION=90., /normal, charsize=blsz

  if keyword_set( PSFILE ) and not keyword_set(NOPLOT) then begin
     x_psclose
     !p.multi=[0,1,1]
  endif
  
  return
end
