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
pro qpq6_lya_stack2, ORIG=orig, RANDOM=random, MEDIAN=median, NODLA=nodla, $
                    PROX=prox, ZLYA=zlya, NEW=new

  compile_opt strictarr

  cd, getenv('QSO_DIR')+'tex/QPQ6/Analysis/N_HI/pro/', curr=curr
  resolve_routine, 'qpq6_lyaew_from_stack', /compile_full_file, /is_function
  resolve_routine, 'qpq6_fit_gauss_stack', /compile_full_file, /is_function
  cd, curr

;  
  if not keyword_set(PSFILE) then psfile = 'qpq6_lya_stack2.ps'
  if not keyword_set( NTOT ) then ntot = 12L
  if not keyword_set( LSIZE ) then lsize = 1.4
  if not keyword_set( LTHICK ) then lthick = 4.
  if not keyword_set( PTHICK ) then pthick = 2.
  if not keyword_set( BLSZ ) then blsz = 1.7
  if not keyword_set( CSIZE ) then csize = 2.8

  if not keyword_set( CONTI ) then conti = 0.8
  if not keyword_set(VMNX) then vmnx = [-5000., 5000]

  ;ew_vmnx = 1350.* [-1, 1]
  NEW = 1

  ;;  Stack files
  stack_root = getenv('QSO_DIR')+'tex/QPQ6/Analysis/Stacks/'
  stack_fil = ['stack_QPQ6_Lya_0_300.fits', $
               'stack_QPQ6_Lya_0_100.fits', $
               'stack_QPQ6_Lya_100_200.fits', $
               'stack_QPQ6_Lya_200_300.fits']
  pref = ''
  if keyword_set(ORIG) then begin
     pref = 'ORIGz_' 
     psfile = 'fig_lya_stacks_orig.ps'
  endif
  if keyword_set(RANDOM) then begin
     pref = 'RAND_' 
     psfile = 'fig_lya_stacks_rand.ps'
  endif
  if keyword_set(NODLA) then begin
     pref = 'NODLA_' 
     psfile = 'fig_lya_stacks_nodla.ps'
  endif
  ;if keyword_set(MEDIAN) then begin
  ;   psfile = 'fig_lya_stacks_median.ps'
  ;   if keyword_set(RANDOM) then psfile = 'fig_lya_stacks_rand_median.ps'
  ;   conti = 0.9
  ;endif
  if keyword_set(PROX) then begin
     pref = 'PROX_' 
     psfile = 'fig_lya_stacks_prox.ps'
  endif
  if keyword_set(ZLYA) then begin
     pref = 'ZLYA_' 
     psfile = 'fig_lya_stacks_zlya.ps'
     vmnx = [-4500, 4500]
  endif
  if keyword_set(NEW) then begin
     stack_fil = ['stack_QPQ6_Lya_0_100.fits', $
                  'stack_QPQ6_Lya_100_200.fits', $
                  'stack_QPQ6_Lya_200_500.fits', $
                  'stack_QPQ6_Lya_500_1000.fits']
     if keyword_set(MEDIAN) then psfile = 'fig_lya_stacks_new_median.ps'
     vmnx = [-4500, 4500]
  endif
  nfil = n_elements(stack_fil)

  npx = 1
  npy = nfil
  if keyword_set( PSFILE ) then x_psopen, psfile, /maxs

  ;;starting at p.multi[0], plot p.multi[1] rows and p.multi[2] columns
  ;;with p.multi[3] z dimensions and going from top to bottom (p.multi[4]
  npy = 2
  !p.multi=[6,2,npy+1,0,1]
  clr = getcolor(/load)
  fclr = clr.lightgray 
  pclr = clr.white 

  ;; LOOP
  for ii=0L,nfil-1 do begin
     if ii EQ 2 then !p.multi=[3,2,npy+1,0,1]
     ;if (ii mod npy) EQ 0 then begin
     ;   if ii GT 0 then offy = (nfil MOD 2) else offy = 0
     ;   !p.multi=[npx*(npy+1) - (npy+1)*(ii/npy) - offy, $
     ;             npx,npy+1,0,1]
     ;endif

     ;; Read
     fx = xmrdfits( stack_root + pref+ stack_fil[ii], 0+keyword_set(MEDIAN), /silen)
     vel = xmrdfits( stack_root + pref+ stack_fil[ii], 3, /silen)
     strct = xmrdfits( stack_root + pref+ stack_fil[ii], 4, /silen)
     if ii EQ 0 then begin
        sv_fx = fx
        sv_vel = vel
     endif
      
     ymnx = [0., 1.05]
          
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
           yrange=ymnx, xtickn=xspaces, xmargin=[7,0], $
           ymargin=[0,0], NODATA=nblnd, $ ;ytickn=yspaces, $
           charsize=csize, psym=10, background=clr.white, color=fclr, $
           xstyle=1, ystyle=ysty, thick=pthick, ytickinterval=ytint,$
           xtitle=xtitle, xtickint=2000.

     ;; Lines
     oplot, [0., 0.], ymnx, color=clr.gray, linestyle=2, thick=1
     ;oplot, [-10000., 10000.], [0.,0.], color=clr.red, linestyle=3,thick=1
     ;oplot, [-10000., 10000.], [1.,1.], color=clr.green,  linestyle=1,thick=1
     ;oplot, [-10000., 10000.], replicate(CONTI, 2), color=clr.green,  linestyle=2,thick=3

     ;if ii GT 0 then oplot, sv_vel, sv_fx, color=clr.gray, psym=10, thick=2
     oplot, vel, fx, color=pclr, psym=10, thick=5
          
     ;; Label
     prs = strsplit(stack_fil[ii], '_', /extract)
     nprs = n_elements(prs)
     R0 = strtrim(long(prs[nprs-2]),2)
     ipos = strpos(prs[nprs-1], '.fits')
     R1 = strmid(prs[nprs-1], 0, ipos)
      xyouts, 0.05*(vmnx[1]-vmnx[0])+vmnx[0], $
              ymnx[0]+ (ymnx[1]-ymnx[0])*0.30, $
              'R=['+R0+','+R1+'] kpc',  color=fclr, charsize=LSIZE
              ;pref+ stack_fil[ii],  color=clr.black, charsize=LSIZE

      ;xyouts, 0.05*(vmnx[1]-vmnx[0])+vmnx[0], $
      ;        ymnx[0]+ (ymnx[1]-ymnx[0])*0.20, $
      ;        '<!9s!X!dz!N> = '+string(round(mean(strct.z_fsig)), format='(i3)')+' km/s', $
      ;        color=clr.black, charsize=LSIZE
      xyouts, 0.05*(vmnx[1]-vmnx[0])+vmnx[0], $
              ymnx[0]+ (ymnx[1]-ymnx[0])*0.15, $
              'N!dspec!N = '+string(n_elements(strct), format='(i3)'), $
              color=fclr, charsize=LSIZE

      ;; EW
      ;EW = qpq6_lyaew_from_stack(vel, fx, ew_vmnx, FIN_CONTI=fin_conti, CLIM=[-4000, 4000]) 
      EW = qpq6_fit_gauss_stack(vel, fx, FIT=fit, CONTI=conti, ACOEFF=acoeff)
      oplot, vel, conti*(1-fit), color=clr.cyan, thick=3
      xpos = 0.65
      xyouts, xpos*(vmnx[1]-vmnx[0])+vmnx[0], $
              ymnx[0]+ (ymnx[1]-ymnx[0])*0.30, $
              'W!dLy!9a!X!N = '+string(EW, format='(f4.2)')+'A', $
              color=fclr, charsize=LSIZE

      xyouts, xpos*(vmnx[1]-vmnx[0])+vmnx[0], $
              ymnx[0]+ (ymnx[1]-ymnx[0])*0.15, $
              '!9s!X = '+string(round(acoeff[2]), format='(i3)')+' km/s', $
              color=fclr, charsize=LSIZE

      ;; Plot the continuum
      ;if not keyword_set(NOPLOT) then oplot, ew_vmnx, fin_conti, color=clr.cyan, thick=2
   endfor 

  ylbl = 0.65
  xyouts, 0.03, ylbl, 'Normalized Flux', $
          alignment=0.5, ORIENTATION=90., /normal, charsize=blsz, color=fclr
  xyouts, 0.53, ylbl, 'Normalized Flux', $
          alignment=0.5, ORIENTATION=90., /normal, charsize=blsz, color=fclr
  if not keyword_set(YPS) then begin
      if ntot EQ 18 then yps = 0.03 else yps = 0.09
  endif
;  xyouts, 0.57, yps, 'Relative Velocity  (km s!u-1!N)', $
;          alignment=0.5, /normal, charsize=blsz

  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]
  
  return
end
