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
pro qpq6_ewstack_vs_rho, PROX=prox, ZLYA=zlya, NEW=new, FSTRCT=fstrct

  compile_opt strictarr

  cd, getenv('QSO_DIR')+'tex/QPQ6/Analysis/N_HI/pro/', curr=curr
  resolve_routine, 'qpq6_lyaew_from_stack', /compile_full_file, /is_function
  resolve_routine, 'qpq6_fit_gauss_stack', /compile_full_file, /is_function
  cd, curr

;  
  if not keyword_set(PSFILE) then psfile = 'qpq6_ewstack_vs_rho.ps'
  if not keyword_set( NTOT ) then ntot = 12L
  if not keyword_set( LSIZE ) then lsize = 1.4
  if not keyword_set( LTHICK ) then lthick = 4.
  if not keyword_set( PTHICK ) then pthick = 2.
  if not keyword_set( BLSZ ) then blsz = 1.7
  if not keyword_set( CSIZE ) then csize = 2.8

  if not keyword_set( CONTI ) then conti = 0.8
  if not keyword_set(VMNX) then vmnx = [-5000., 5000]

  ;ew_vmnx = 1350.* [-1, 1]

  ;;  Stack files
  stack_root = getenv('QSO_DIR')+'tex/QPQ6/Analysis/Stacks/'
  stack_fil = ['stack_QPQ6_Lya_0_100.fits', $
               'stack_QPQ6_Lya_100_200.fits', $
               'stack_QPQ6_Lya_200_500.fits', $
               'stack_QPQ6_Lya_500_1000.fits']
  if not keyword_set( NBOOT ) then nboot = 100L
  boot_fil = 'BOOT_'+stack_fil 
  pref = ''
  vmnx = [-4500, 4500]
  nfil = n_elements(stack_fil)

  ;; Do the measurements
  if not keyword_set(FSTRCT) then begin
     EWval = fltarr(nfil)
     Rval = fltarr(nfil)
     Rmnx = fltarr(nfil,2)
     sigEW = fltarr(nfil)
     for ii=0L,nfil-1 do begin
        ;; Read
        fx = xmrdfits( stack_root + pref+ stack_fil[ii], 0+keyword_set(MEDIAN), /silen)
        vel = xmrdfits( stack_root + pref+ stack_fil[ii], 3, /silen)
        strct = xmrdfits( stack_root + pref+ stack_fil[ii], 4, /silen)
        
        ;; R
        Rval[ii] = median(strct.R)
        mn = min(strct.R, max=mx)
        Rmnx[ii,*] = [mn,mx]
        
        ;; EW
        EWval[ii] = qpq6_fit_gauss_stack(vel, fx, FIT=fit, CONTI=conti, ACOEFF=acoeff)
        
        ;; Bootstrap
        velb = xmrdfits( stack_root + boot_fil[ii], 0, /silen)
        npixb =n_elements(velb)
        strct = xmrdfits( stack_root + boot_fil[ii], 1, /silen)
        sv_boot = fltarr(npixb, nboot)
        sv_bootEW = fltarr(nboot)
        for jj=0L,nboot-1 do begin
           bootfx = xmrdfits( stack_root + boot_fil[ii], 2+jj, /silen)
           bootEW = qpq6_fit_gauss_stack(velb, bootfx, FIT=fit, CONTI=conti, ACOEFF=acoeff)
           sv_boot[*,jj] = bootfx
           sv_bootEW[jj] = bootEW
        endfor
        
        ;; Stats
        djs_iterstat, sv_bootEW, sigma=sig 
        sigEW[ii] = sig
     endfor
     
     fstrct = { $
              Rval: Rval, $
              Rmnx: Rmnx, $
              EWval: EWval, $
              sigEW: sigEW $
              }
  endif else begin
     Rval = fstrct.rval
     Rmnx = fstrct.Rmnx
     EWval = fstrct.EWval
     sigEW = fstrct.sigEW
  endelse

  if keyword_set( PSFILE ) then x_psopen, psfile, /maxs

  ;;starting at p.multi[0], plot p.multi[1] rows and p.multi[2] columns
  ;;with p.multi[3] z dimensions and going from top to bottom (p.multi[4]
  clr = getcolor(/load)
  fclr = clr.white
  pclr = clr.cyan

  xtitle='R!dphys!N (kpc)'
  ytitle='W!dLya!N (Ang)'
  ymnx = [0., 3.3]
          
  plot, [0], [0], $
        xrange=[0., 1000], $
        yrange=ymnx, xtickn=xspaces, xmargin=[8,11], $
        ymargin=[4,1], /NODATA, $ ;ytickn=yspaces, $
        charsize=csize, psym=10, background=clr.white, color=fclr, $
        xstyle=1, ystyle=1, thick=pthick, ytickinterval=ytint, xtitle=xtitle, ytitle=ytitle

  ;; Label
  ;xpos = 0.65
  ;xyouts, xpos*(vmnx[1]-vmnx[0])+vmnx[0], $
  ;        ymnx[0]+ (ymnx[1]-ymnx[0])*0.30, $
  ;        'W!dLy!9a!X!N = '+string(EW, format='(f4.2)')+'A', $
  ;        color=clr.black, charsize=LSIZE

  ;; Plot
  oploterror, Rval, EWval, (Rval-Rmnx[*,0]), sigEW, color=pclr, errcol=pclr, psym=1, /lob
  oploterror, Rval, EWval, abs(Rval-Rmnx[*,1]), sigEW, color=pclr, errcol=pclr, psym=1, /hib

  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]
  
  return
end
