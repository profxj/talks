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
;; fig_lya_stacks, /NEW
pro tpe_full_stack, NOPLOT=noplot, CSIZE=csize, LSIZE=lsize, XMRG=xmrg, $
                        MED=med

  compile_opt strictarr

  cd, getenv('QSO_DIR')+'tex/QPQ6/Analysis/N_HI/pro/', curr=curr
  resolve_routine, 'qpq6_lyaew_from_stack', /compile_full_file, /is_function
  resolve_routine, 'qpq6_fit_gauss_stack', /compile_full_file, /is_function
  cd, curr

;  
  if not keyword_set(PSFILE) then psfile = 'tpe_full_stack.ps'
  if keyword_set(MED) then psfile = 'tpe_full_stack_median.ps'
  if not keyword_set( NTOT ) then ntot = 12L
  if not keyword_set( LSIZE ) then lsize = 1.2
  if not keyword_set( LTHICK ) then lthick = 4.
  if not keyword_set( PTHICK ) then pthick = 2.
  if not keyword_set( BLSZ ) then blsz = 1.5
  if not keyword_set( CSIZE ) then csize = 2.5

  if not keyword_set( CONTI ) then conti = 0.8
  if not keyword_set(VMNX) then vmnx = [-5000., 5000]

  ew_vmnx = 1350.* [-1, 1]
  if not keyword_set(XMRG) then xmrg = [7,20]

  ;;  Stack files
  stack_root = getenv('QSO_DIR')+'tex/TPE/Letter/Analysis/Stacks/'
  stack_fil = ['stack_TPE_Lya_0_3000.fits', $
               'stack_TPE_Lya_0_300.fits', $
               'stack_TPE_Lya_300_600.fits', $
               'stack_TPE_Lya_600_1200.fits', $
               'stack_TPE_Lya_1200_2400.fits', $
               'stack_TPE_Lya_2400_3000.fits']
;               'stack_TPE_Lya_0_300.fits', $
;               'stack_TPE_Lya_300_1000.fits', $
;               'stack_TPE_Lya_1000_3000.fits']
  nfil = n_elements(stack_fil)
  pref = ''

  npx = 1
  npy = nfil
  if keyword_set( PSFILE )  and not keyword_set(NOPLOT) then begin 
     x_psopen, psfile, /portrait
     !p.multi=[npx*(npy+1),npx,npy+1,0,1]
  endif
  clr = getcolor(/load)
  fclr = clr.white

  ;; LOOP
  for ss=0,0 do begin

     for ii=0L,nfil-1 do begin

     ;; Read
     fx = xmrdfits( stack_root + pref+ stack_fil[ii], 0+keyword_set(MED), /silen)
     vel = xmrdfits( stack_root + pref+ stack_fil[ii], 3, /silen)
     strct = xmrdfits( stack_root + pref+ stack_fil[ii], 4, /silen)
      
     ymnx = [0.4, 1.05]
          
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
          
     plot, [0], [0], $
           xrange=vmnx, $
           yrange=ymnx, xtickn=xspaces, xmargin=xmrg, $
           ymargin=[0,0], NODATA=nblnd, $ ;ytickn=yspaces, $
           charsize=csize, psym=10, background=clr.white, color=fclr, $
           xstyle=1, ystyle=ysty, thick=pthick, ytickinterval=ytint,$
           xtitle=xtitle, xtickint=2000.

     ;; Lines
     oplot, [0., 0.], ymnx, color=clr.lightgray, linestyle=2, thick=3
     ;oplot, [-10000., 10000.], [0.,0.], color=clr.red, linestyle=3,thick=1
     ;oplot, [-10000., 10000.], [1.,1.], color=clr.green,  linestyle=1,thick=1
     ;oplot, [-10000., 10000.], replicate(CONTI, 2), color=clr.green,  linestyle=2,thick=3

     ;if ii GT 0 then oplot, sv_vel, sv_fx, color=clr.gray, psym=10, thick=2
     oplot, vel, fx, color=fclr, psym=10, thick=5

     ;; Label
     ;if ii EQ 0 then xyouts, 2000., 0.90, 'Mean', charsi=lsize, color=clr.black $ 
     ;else xyouts, 2000., 0.90, 'Median', charsi=lsize, color=clr.black 
     

     if ss EQ 0 then begin
        ;; Label
        ;xyouts, 0.05*(vmnx[1]-vmnx[0])+vmnx[0], $
        ;        ymnx[0]+ (ymnx[1]-ymnx[0])*0.10, $
        ;        pref+ stack_fil[ii],  color=clr.black, charsize=LSIZE
                                ;xyouts, 0.05*(vmnx[1]-vmnx[0])+vmnx[0], $
                                ;        ymnx[0]+ (ymnx[1]-ymnx[0])*0.20, $
                                ;        '<!9s!X!dz!N> = '+string(round(mean(strct.z_fsig)), format='(i3)')+' km/s', $
                                ;        color=clr.black, charsize=LSIZE
        avgR = mean(strct.R)
        thisletter = byte(94)
        perpletter = '!9' + string(thisletter) + '!X'
        xyouts, 0.70*(vmnx[1]-vmnx[0])+vmnx[0], $
                ymnx[0]+ (ymnx[1]-ymnx[0])*0.40, $
                'N!dspec!N = '+string(n_elements(strct), format='(i4)'), $
                color=fclr, charsize=LSIZE
        xyouts, 0.70*(vmnx[1]-vmnx[0])+vmnx[0], $
                ymnx[0]+ (ymnx[1]-ymnx[0])*0.20, $
                '<R!d'+perpletter+'!N> = '+string(round(avgR),'(i4)')+' kpc', $
                color=fclr, charsize=LSIZE
        
        ;; EW
                                ;EW = qpq6_lyaew_from_stack(vel, fx, ew_vmnx, FIN_CONTI=fin_conti, CLIM=[-4000, 4000]) 
        EW = qpq6_fit_gauss_stack(vel, fx, FIT=fit, CONTI=conti, ACOEFF=acoeff)
        oplot, vel, conti*(1-fit), color=clr.cyan, thick=3
        xyouts, 0.05*(vmnx[1]-vmnx[0])+vmnx[0], $
                ymnx[0]+ (ymnx[1]-ymnx[0])*0.40, $
                'W!d1350!N = '+string(EW, format='(f4.2)')+'A', $
                color=fclr, charsize=LSIZE
        
        xyouts, 0.05*(vmnx[1]-vmnx[0])+vmnx[0], $
                ymnx[0]+ (ymnx[1]-ymnx[0])*0.20, $
                '!9s!X = '+string(round(acoeff[2]), format='(i3)')+' km/s', $
                color=fclr, charsize=LSIZE
     endif else begin ;;

        if ii EQ 1000 then begin
           ;; Plot Bootstrap
           velb = xmrdfits( stack_root + boot_fil[ii], 0, /silen)
           npixb =n_elements(velb)
           strct = xmrdfits( stack_root + boot_fil[ii], 1, /silen)
           sv_boot = fltarr(npixb, nboot)
           sv_bootEW = fltarr(npixb, nboot)
           for jj=0L,nboot-1 do begin
              bootfx = xmrdfits( stack_root + boot_fil[ii], 2+jj, /silen)
              if not keyword_set(NOPLOT) then oplot, velb, bootfx, color=clr.lightgray, thick=2
              bootEW = qpq6_fit_gauss_stack(velb, bootfx, FIT=fit, CONTI=conti, ACOEFF=acoeff)
              sv_boot[*,jj] = bootfx
              sv_bootEW[*,jj] = bootEW
           endfor
        endif

        ;; Random
        rand_fx = xmrdfits( stack_root + rand_fil[ii], ii, /silen)
        rand_vel = xmrdfits( stack_root + rand_fil[ii], 3, /silen)
        if not keyword_set(NOPLOT) then oplot, rand_vel, rand_fx, color=clr.cyan, psym=10, thick=2

        ;; Replot data
        oplot, vel, fx, color=fclr, psym=10, thick=5
     endelse
   endfor 
  endfor

  xyouts, 0.03, 0.55, 'Normalized Flux', $
          alignment=0.5, ORIENTATION=90., /normal, charsize=blsz, color=fclr
  if not keyword_set(YPS) then begin
      if ntot EQ 18 then yps = 0.03 else yps = 0.09
  endif
;  xyouts, 0.57, yps, 'Relative Velocity  (km s!u-1!N)', $
;          alignment=0.5, /normal, charsize=blsz

  if keyword_set( PSFILE ) and not keyword_set(NOPLOT) then begin
     x_psclose
     !p.multi=[0,1,1]
  endif
  
  return
end
