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
pro qpq6_all_stack, ORIG=orig, RANDOM=random, MEDIAN=median, NODLA=nodla, $
                    PROX=prox, ZLYA=zlya, NEW=new

  compile_opt strictarr

  cd, getenv('QSO_DIR')+'tex/QPQ6/Analysis/N_HI/pro/', curr=curr
  resolve_routine, 'qpq6_lyaew_from_stack', /compile_full_file, /is_function
  resolve_routine, 'qpq6_fit_gauss_stack', /compile_full_file, /is_function
  cd, curr

;  
  if not keyword_set(PSFILE) then psfile = 'qpq6_all_stack.ps'
  if not keyword_set( NTOT ) then ntot = 12L
  if not keyword_set( LSIZE ) then lsize = 1.5
  if not keyword_set( LTHICK ) then lthick = 4.
  if not keyword_set( PTHICK ) then pthick = 2.
  if not keyword_set( BLSZ ) then blsz = 1.7
  if not keyword_set( CSIZE ) then csize = 2.5

  if not keyword_set( CONTI ) then conti = 0.8
  if not keyword_set(VMNX) then vmnx = [-4500., 4500]

  ew_vmnx = 1350.* [-1, 1]

  ;;  Stack files
  stack_root = getenv('QSO_DIR')+'tex/QPQ6/Analysis/Stacks/'
  stack_fil = ['stack_QPQ6_Lya_0_1000.fits', $
               'stack_QPQ6_Lya_0_1000.fits']
  if not keyword_set( NBOOT ) then nboot = 100L
  boot_fil = ['BOOT_stack_QPQ6_Lya_0_1000.fits', $
               'BOOT_stack_QPQ6_Lya_0_1000.fits']
  rand_fil = ['RAND_stack_QPQ6_Lya_0_1000.fits', $
               'RAND_stack_QPQ6_Lya_0_1000.fits']
  nfil = n_elements(stack_fil)
  pref = ''

  npx = 2
  npy = nfil
  if keyword_set( PSFILE ) then x_psopen, psfile, /maxs

  ;;starting at p.multi[0], plot p.multi[1] rows and p.multi[2] columns
  ;;with p.multi[3] z dimensions and going from top to bottom (p.multi[4]
  !p.multi=[npx*(npy+1),npx,npy+1,0,1]
  clr = getcolor(/load)
  fclr = clr.white

  ;; LOOP on Pages
  for qq=0,4 do begin
     case qq of
        0: begin
           nss = 0
           nfil = 1
        end
        1: begin
           nss = 0
           nfil = 2
        end
        2:
        3: nss=1
        4:
        else: stop
     endcase

     for ss=0,nss do begin
        
        for ii=0L,nfil-1 do begin
           
           if ii EQ 0 then begin
              if ss EQ 0 then !p.multi=[0, npx,npy+1,0,1] $
              else !p.multi=[(npy+1), npx,npy+1,0,1] 
           endif
           
           ;; Read
           fx = xmrdfits( stack_root + pref+ stack_fil[ii], ii, /silen)
           vel = xmrdfits( stack_root + pref+ stack_fil[ii], 3, /silen)
           strct = xmrdfits( stack_root + pref+ stack_fil[ii], 4, /silen)
           
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
                 yrange=ymnx, xtickn=xspaces, xmargin=[7,2], $
                 ymargin=[0,0], NODATA=nblnd, $ ;ytickn=yspaces, $
                 charsize=csize, psym=10, background=clr.white, color=fclr, $
                 xstyle=1, ystyle=ysty, thick=pthick, ytickinterval=ytint,$
                 xtitle=xtitle, xtickint=2000.
           
           ;; Lines
           oplot, [0., 0.], ymnx, color=clr.cyan, linestyle=2, thick=3
                                ;oplot, [-10000., 10000.], [0.,0.], color=clr.red, linestyle=3,thick=1
           oplot, [-10000., 10000.], [1.,1.], color=clr.green,  linestyle=2,thick=4
     ;oplot, [-10000., 10000.], replicate(CONTI, 2), color=clr.green,  linestyle=2,thick=3
           
                                ;if ii GT 0 then oplot, sv_vel, sv_fx, color=clr.gray, psym=10, thick=2
           oplot, vel, fx, color=fclr, psym=10, thick=5
           
           ;; Label
           xpos = 1500.
           xyouts, xpos, 0.30, 'R!dphys!N < 1 Mpc', charsi=lsize, color=fclr, align=0. 
           ypos = 0.15
           if ii EQ 0 then xyouts, xpos, ypos, 'Mean', charsi=lsize, color=fclr, align=0. $ 
           else xyouts, xpos, ypos, 'Median', charsi=lsize, color=fclr, align=0. 
     

           if ss EQ 0 and qq GT 1 then begin
              EW = qpq6_fit_gauss_stack(vel, fx, FIT=fit, CONTI=conti, ACOEFF=acoeff)
              oplot, vel, conti*(1-fit), color=clr.red, thick=3
              xyouts, 0.05*(vmnx[1]-vmnx[0])+vmnx[0], $
                      ymnx[0]+ (ymnx[1]-ymnx[0])*0.30, $
                      'W!dLy!9a!X!N = '+string(EW, format='(f4.2)')+'A', $
                      color=fclr, charsize=LSIZE
              
              xyouts, 0.05*(vmnx[1]-vmnx[0])+vmnx[0], $
                      ymnx[0]+ (ymnx[1]-ymnx[0])*0.15, $
                      '!9s!X = '+string(round(acoeff[2]), format='(i3)')+' km/s', $
                      color=fclr, charsize=LSIZE
           endif else begin ;;
              
              if ii EQ 0 and qq GT 3 then begin
                 ;; Plot Bootstrap
                 velb = xmrdfits( stack_root + boot_fil[ii], 0, /silen)
                 npixb =n_elements(velb)
                 strct = xmrdfits( stack_root + boot_fil[ii], 1, /silen)
                 sv_boot = fltarr(npixb, nboot)
                 sv_bootEW = fltarr(nboot)
                 for jj=0L,nboot-1 do begin
                    bootfx = xmrdfits( stack_root + boot_fil[ii], 2+jj, /silen)
                    if not keyword_set(NOPLOT) then oplot, velb, bootfx, color=clr.darkgray, thick=2
                    bootEW = qpq6_fit_gauss_stack(velb, bootfx, FIT=fit, CONTI=conti, ACOEFF=acoeff)
                    sv_boot[*,jj] = bootfx
                    sv_bootEW[jj] = bootEW
                 endfor
                 ;; Error
                 djs_iterstat, sv_bootEW, sigma=sig 
                 xyouts, 0.05*(vmnx[1]-vmnx[0])+vmnx[0], $
                         ymnx[0]+ (ymnx[1]-ymnx[0])*0.30, $
                         '!9s!X(W!dLy!9a!X!N) = '+string(sig, format='(f4.2)')+'A', $
                         color=fclr, charsize=LSIZE
              endif

              if qq GT 2 then begin
                 ;; Random
                 rand_fx = xmrdfits( stack_root + rand_fil[ii], ii, /silen)
                 rand_vel = xmrdfits( stack_root + rand_fil[ii], 3, /silen)
                 if not keyword_set(NOPLOT) then oplot, rand_vel, rand_fx, color=clr.cyan, psym=10, thick=4
              endif

              ;; Replot data
              oplot, vel, fx, color=fclr, psym=10, thick=5
           endelse
        endfor 
     endfor
  endfor

  xyouts, 0.03, 0.65, 'Normalized Flux', $
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
