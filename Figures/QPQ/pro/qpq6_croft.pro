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
pro qpq6_croft, ORIG=orig, RANDOM=random, MEDIAN=median, NODLA=nodla, $
                      PROX=prox, ZLYA=zlya, NEW=new, KLUDGE=kludge

  compile_opt strictarr

  cd, getenv('QSO_DIR')+'tex/QPQ6/Analysis/N_HI/pro/', curr=curr
  resolve_routine, 'qpq6_lyaew_from_stack', /compile_full_file, /is_function
  resolve_routine, 'qpq6_fit_gauss_stack', /compile_full_file, /is_function
  cd, curr

;  
  if not keyword_set(PSFILE) then psfile = 'qpq6_croft.ps'
  if not keyword_set( NTOT ) then ntot = 12L
  if not keyword_set( LSIZE ) then lsize = 1.2
  if not keyword_set( LTHICK ) then lthick = 4.
  if not keyword_set( PTHICK ) then pthick = 2.
  if not keyword_set( BLSZ ) then blsz = 1.7
  if not keyword_set( CSIZE ) then csize = 2.9
  if not keyword_set( KLUDGE ) then kludge = 1.

  if not keyword_set( CONTI ) then conti = 0.8
  if not keyword_set(LIT_H) then LIT_H = 0.70  ;; (Hubble)

  ew_vmnx = 1350.* [-1, 1]

  ;;  Stack files
  stack_root = getenv('QSO_DIR')+'tex/QPQ6/Analysis/Stacks/'
  ;stack_fil = ['cMpc_stack_QPQ6_Lya_1_2.fits', $
  ;             'cMpc_stack_QPQ6_Lya_2_3.fits', $
  ;             'cMpc_stack_QPQ6_Lya_3_4.fits']
  stack_fil = ['cMpc_stack_QPQ6_Lya_1.0_1.5.fits', $
               'cMpc_stack_QPQ6_Lya_1.5_2.0.fits', $
               'cMpc_stack_QPQ6_Lya_2.0_2.5.fits']
  pref = ''
  vmnx = [-4000, 4000]
  xmnx = [-50, 50.] ; cMpc
  nfil = n_elements(stack_fil)

  npx = 1
  npy = nfil
  if keyword_set( PSFILE ) then x_psopen, psfile, /portrait

  ;;starting at p.multi[0], plot p.multi[1] rows and p.multi[2] columns
  ;;with p.multi[3] z dimensions and going from top to bottom (p.multi[4]
  !p.multi=[npx*(npy+1),npx,npy+1,0,1]
  clr = getcolor(/load)
  pclr = clr.white

  thisletter = byte(94)
  perpletter = '!9' + string(thisletter) + '!X'

  ;; LOOP
  for ii=0L,nfil-1 do begin
     if (ii mod npy) EQ 0 then begin
        if ii GT 0 then offy = (nfil MOD 2) else offy = 0
        !p.multi=[npx*(npy+1) - (npy+1)*(ii/npy) - offy, $
                  npx,npy+1,0,1]
     endif

     ;; Read
     fx = xmrdfits( stack_root + pref+ stack_fil[ii], 0+keyword_set(MEDIAN), /silen)
     vel = xmrdfits( stack_root + pref+ stack_fil[ii], 3, /silen)
     strct = xmrdfits( stack_root + pref+ stack_fil[ii], 4, /silen)
     npix = n_elements(vel)
     npair = n_elements(strct)
      
     ymnx = [0.5, 0.85]
          
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
        ;xtitle='s (cMpc h!u-1!N)'
        xtitle='Relative Velocity (km/s)'
     endelse
     if (ii GT npy-1) then yspaces = replicate(' ',30) 
          
     plot, [0], [0], $
           xrange=vmnx, $
           yrange=ymnx, xtickn=xspaces, xmargin=[7,19], $
           ymargin=[0,0], NODATA=nblnd, $ ;ytickn=yspaces, $
           charsize=csize, psym=10, background=clr.white, color=pclr, $
           xstyle=1, ystyle=ysty, thick=pthick, ytickinterval=ytint,$
           xtitle=xtitle;, xtickint=1000.

     ;; Lines
     oplot, [0., 0.], ymnx, color=clr.gray, linestyle=2, thick=1
     ;oplot, [-10000., 10000.], [0.,0.], color=clr.red, linestyle=3,thick=1
     ;oplot, [-10000., 10000.], [1.,1.], color=clr.green,  linestyle=1,thick=1
     ;oplot, [-10000., 10000.], replicate(CONTI, 2), color=clr.green,  linestyle=2,thick=3

     ;; Label
      ;xyouts, 0.05*(xmnx[1]-xmnx[0])+xmnx[0], $
      ;        ymnx[0]+ (ymnx[1]-ymnx[0])*0.10, $
      ;        pref+ stack_fil[ii],  color=clr.black, charsize=LSIZE


     MODEL=1
     if keyword_set(MODEL) then begin
        dVel = 10. ;; km/s
        nstep = (vmnx[1]-vmnx[0])/dVel + 1;; every 10km/s
        vmodel = vmnx[0] + dVel*findgen(nstep)

        all_model = fltarr(nstep, npair)
        all_smooth = fltarr(nstep, npair)
        all_smooth2 = fltarr(nstep, npair) ;; Low Mass

        for kk=0L,npair-1 do begin

           ;; Generate model in r space
           Rcom = (strct[kk].zfg+1)*strct[kk].R / 1e3 ;; cMpc 
           zval = strct[kk].zfg + vmodel * (1+strct[kk].zfg) / 3e5
           Dval = cosm_dist(zval, /w05map) ; cMpc
           mn = min(abs(zval-strct[kk].zfg), imn)
           rval = sqrt( Rcom^2 + (Dval-Dval[imn])^2 )  ; cMpc

           teff = taueff_evo(strct[kk].zfg)

           ;; Kim & Croft (12.85)
           ;r0 = 3.33
           ;gamma = 1.51
        
           ;; Eyeball of KC08 for 12.5
           r0 = 2.9
           gamma = 1.43

           r0_2 = 2.0
           gamma_2 = 1.2
        
           Fr = exp(-1*teff*(1+(rval/r0)^(-1*gamma)))
           Fr2 = exp(-1*teff*(1+(rval/r0_2)^(-1*gamma_2)))
           all_model[*,kk] = Fr
        
           ;; Smooth
           kernel = gauss_kernel(strct[kk].z_fsig / dVel)
           smooth_Fr = convol(Fr, kernel, /edge_wrap)
           smooth_Fr2 = convol(Fr2, kernel, /edge_wrap)
           all_smooth[*,kk] = smooth_Fr
           all_smooth2[*,kk] = smooth_Fr2

        endfor
        ;; Sum
        final_Fr = total(all_model, 2) / float(npair)
        final_smooth = total(all_smooth, 2) / float(npair)
        final_smooth2 = total(all_smooth2, 2) / float(npair)

        ;; Interpolate onto Stacks
        ;inter_Fr = interpol(final_Fr, vmodel, vel)
        x_specrebin, vmodel, final_Fr, vel, inter_Fr, /flamb 
        ;inter_sFr = interpol(final_smooth, vmodel, vel)
        x_specrebin, vmodel, final_smooth, vel, inter_sFr, /flamb 
        ;inter_sFr2 = interpol(final_smooth2, vmodel, vel)
        x_specrebin, vmodel, final_smooth2, vel, inter_sFr2, /flamb 

        ;x_psclose
        ;!p.multi=[0,1,1]
        ;stop

        ;; Scale
        case ii of
           0: SCALE = 0.95
           1: SCALE = 0.95
           2: SCALE = 0.97
           else: stop
        endcase

        ;; Modify the edges
        nset = 14L
        inter_Fr[0:nset-1] = inter_Fr[nset]
        inter_sFr[0:nset-1] = inter_sFr[nset]
        inter_sFr2[0:nset-1] = inter_sFr2[nset]
        nset = 86L
        inter_Fr[nset+1:*] = inter_Fr[nset]
        inter_sFr[nset+1:*] = inter_sFr[nset]
        inter_sFr2[nset+1:*] = inter_sFr2[nset]

        ;; Plot
        if keyword_set(SHOW_NOSMOOTH) then $
           oplot, vel, inter_Fr*SCALE, color=clr.cyan, linesty=1
        oplot, vel, inter_sFr*SCALE, color=clr.cyan
        oplot, vel, inter_sFr2*SCALE, color=clr.yellow, linesty=2
     ;;
     endif

     ;; Data
     oplot, vel, fx, color=pclr, psym=10, thick=3
     Rcom = mean( strct.R*(strct.zfg+1) ) * LIT_H/ 1e3
     xyouts, 0.02*(vmnx[1]-vmnx[0])+vmnx[0], $
             ymnx[0]+ (ymnx[1]-ymnx[0])*0.10, $
             '<R!S!d'+perpletter+'!R!ucom!N> = '+ $
             string(Rcom, format='(f4.2)')+' cMpc h!u-1!N', $
             color=pclr, charsize=LSIZE

      ;; Plot the continuum
      ;if not keyword_set(NOPLOT) then oplot, ew_vmnx, fin_conti, color=clr.cyan, thick=2
   endfor 

  xyouts, 0.04, 0.60, 'Normalized Flux', $
          alignment=0.5, ORIENTATION=90., /normal, charsize=csize*0.5, color=pclr
  if not keyword_set(YPS) then begin
      if ntot EQ 18 then yps = 0.03 else yps = 0.09
  endif
;  xyouts, 0.57, yps, 'Relative Velocity  (km s!u-1!N)', $
;          alignment=0.5, /normal, charsize=blsz

  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]
  
  return
end
