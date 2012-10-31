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
pro qpq6_r0gamma, ORIG=orig, RANDOM=random, MEDIAN=median, NODLA=nodla, $
                      PROX=prox, ZLYA=zlya, NEW=new, KLUDGE=kludge

  compile_opt strictarr

  cd, getenv('QSO_DIR')+'tex/QPQ6/Analysis/N_HI/pro/', curr=curr
  resolve_routine, 'qpq6_lyaew_from_stack', /compile_full_file, /is_function
  resolve_routine, 'qpq6_fit_gauss_stack', /compile_full_file, /is_function
  cd, curr

;  
  if not keyword_set(PSFILE) then psfile = 'qpq6_r0gamma.ps'
  if not keyword_set( NTOT ) then ntot = 12L
  if not keyword_set( LSIZE ) then lsize = 1.2
  if not keyword_set( LTHICK ) then lthick = 4.
  if not keyword_set( PTHICK ) then pthick = 2.
  if not keyword_set( BLSZ ) then blsz = 1.7
  if not keyword_set( CSIZE ) then csize = 2.5
  if not keyword_set( KLUDGE ) then kludge = 0.015

  if not keyword_set( CONTI ) then conti = 0.8

  ew_vmnx = 1350.* [-1, 1]

  ;;  Stack files
  stack_root = getenv('QSO_DIR')+'tex/QPQ6/Analysis/Stacks/'
  stack_fil = ['stack_QPQ6_Lya_200_500.fits', $
               'stack_QPQ6_Lya_500_1000.fits']
  pref = ''
  vmnx = [-5000, 5000]
  xmnx = [-50, 50.] ; cMpc
  nfil = n_elements(stack_fil)

  npx = 1
  npy = nfil
  if keyword_set( PSFILE ) then x_psopen, psfile, /portrait

  ;;starting at p.multi[0], plot p.multi[1] rows and p.multi[2] columns
  ;;with p.multi[3] z dimensions and going from top to bottom (p.multi[4]
  !p.multi=[npx*(npy+1),npx,npy+1,0,1]
  clr = getcolor(/load)
  fclr = clr.white

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
        xtitle='s (cMpc h!u-1!N)'
     endelse
     if (ii GT npy-1) then yspaces = replicate(' ',30) 
          
     plot, [0], [0], $
           xrange=xmnx, $
           yrange=ymnx, xtickn=xspaces, xmargin=[7,10], $
           ymargin=[0,0], NODATA=nblnd, $ ;ytickn=yspaces, $
           charsize=csize, psym=10, background=clr.white, color=fclr, $
           xstyle=1, ystyle=ysty, thick=pthick, ytickinterval=ytint,$
           xtitle=xtitle;, xtickint=1000.

     ;; Lines
     oplot, [0., 0.], ymnx, color=clr.yellow, linestyle=2, thick=3
     ;oplot, [-10000., 10000.], [0.,0.], color=clr.red, linestyle=3,thick=1
     ;oplot, [-10000., 10000.], [1.,1.], color=clr.green,  linestyle=1,thick=1
     ;oplot, [-10000., 10000.], replicate(CONTI, 2), color=clr.green,  linestyle=2,thick=3

     ;; Label
      ;xyouts, 0.05*(xmnx[1]-xmnx[0])+xmnx[0], $
      ;        ymnx[0]+ (ymnx[1]-ymnx[0])*0.10, $
      ;        pref+ stack_fil[ii],  color=clr.black, charsize=LSIZE
      ;xyouts, 0.05*(xmnx[1]-xmnx[0])+xmnx[0], $
      ;        ymnx[0]+ (ymnx[1]-ymnx[0])*0.20, $
      ;        '<!9s!X!dz!N> = '+string(round(mean(strct.z_fsig)), format='(i3)')+' km/s', $
      ;        color=clr.black, charsize=LSIZE
      ;xyouts, 0.05*(xmnx[1]-xmnx[0])+xmnx[0], $
      ;        ymnx[0]+ (ymnx[1]-ymnx[0])*0.30, $
      ;        'N!dspec!N = '+string(n_elements(strct), format='(i3)'), $
      ;        color=clr.black, charsize=LSIZE

      ;; Get <R>
     avgR = mean(strct.R)  ;; Physical
     avgz = mean(strct.zfg)  
     Rcom = avgR * (1+avgz) /1e3 ;; comoving Mpc

      ;; Map v -> R
      nplt = 20001L
      vplt = -10000. + dindgen(nplt)
      zplt = avgz + vplt / 3e5 * (1+avgz)
      Dplt = cosm_dist(zplt, /w05map)  ; cMpc
      Dplt = Dplt - Dplt[(nplt-1)/2] ;; 
      rplt = sqrt( Rcom^2 + Dplt^2 ) ; cMpc

     ;; Stack
     Dstack = interpol(Dplt, vplt, vel)
     oplot, Dstack, fx, color=fclr, psym=10, thick=5
          
     ;; Continuum
     if not keyword_set(CLIM) then clim = [-3000., 3000] ; km/s
     px = where( abs( abs(Dstack)-40 ) LT 10)
     ;low_c = where(vel LT CLIM[0])
     ;low_c = low_c[1:*]  ;; Trim the first pixel
     ;hi_c = where(vel GT CLIM[1])
     ;hi_c = hi_c[0:n_elements(hi_c)-2] ;; Trim the last pixel
     ;teff = -1*alog(median([fx[low_c], fx[hi_c]]))
     teff = -1*alog(median(fx[px] + KLUDGE*(ii EQ 0))) ;; KLUDGE!

     ;; Model

     ;; Kim & Croft (12.85)
     r0 = 3.33
     gamma = 1.51

     ;; Eyeball of KC08 for 12.5
     r0 = 2.9
     gamma = 1.43

     Fr = exp(-1*teff*(1+(rplt/r0)^(-1*gamma)))
     ;Fr = exp(-1*teff*(1+(rplt/3.3)^(-1.51)))

     oplot, Dplt, Fr, color=clr.cyan, linesty=1

     ;; Smooth
     Dsmooth = interpol(Dplt, vplt, 500.)
     dDplt = median( abs(Dplt - shift(Dplt,1) ) )
     kernel = gauss_kernel(Dsmooth / dDplt)
     smooth_Fr = convol(Fr, kernel)
     oplot, Dplt, smooth_Fr, color=clr.cyan
     ;print, total( exp(-1.*teff) - Fr ), total(exp(-1.*teff)-smooth_Fr)
     ;stop

     ;;
      xyouts, 0.05*(xmnx[1]-xmnx[0])+xmnx[0], $
              ymnx[0]+ (ymnx[1]-ymnx[0])*0.20, $
              '<R> = '+string(Rcom, format='(f4.2)')+' cMpc h!u-1!N', $
              color=fclr, charsize=LSIZE

      ;; Plot the continuum
      ;if not keyword_set(NOPLOT) then oplot, ew_vmnx, fin_conti, color=clr.cyan, thick=2
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
