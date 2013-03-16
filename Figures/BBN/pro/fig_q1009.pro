;+ 
; NAME:
; fig_lyseries
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
pro fig_q1009, datfil, NTOT=NTOT, CSIZE=csize, PSFILE=psfile, XTINT=xtint, $
              FINE=fine, INFIL=infil, LSIZE=lsize, YPS=yps, LTHICK=lthick

; lowzovi_prsdat -- Reads in DLA data to a structure

;if (N_params() LT 0) then begin 
;    print,'Syntax - ' + $
;      'fuse_velplt, strct, instr_list, vel_fil, NTOT=, CSIZE=, ' + $
;      'LSIZE=, PSFILE=, XTINT= (v1.1)' 
;    return
;endif 
;  
  if not keyword_set(PSFILE) then psfile = 'fig_q1009.ps'
  if not keyword_set( NTOT ) then ntot = 16L
  if not keyword_set( LSIZE ) then lsize = 1.8
  if not keyword_set( LTHICK ) then lthick = 6.
  if not keyword_set( BLSZ ) then blsz = 1.4
  if not keyword_set( CSIZE ) then csize = 1.7
  if not keyword_set( CSZ2 ) then csz2 = 1.4

;  if not keyword_set( PSFIL ) then psfil = 'tmp.ps'

;  Read instrument file list
;readcol, instr_list, instr_fil, inst_dw, inst_w0, format='a,f,f'



  ;; Open vel_fil
  close, /all
  datfil = '~/Keck/HIRES/RedData/Q1009+2956/Q1009+2956a_f.fits'
  FWHM = 4.
  zhabs = 2.50369d
  fx = x_readspec(datfil, wav=wave, npix=npix) 

  if keyword_set( PSFILE ) then x_psopen, psfile, /maxs
  clr = getcolor(/load)
  sclr = clr.gray

  ;; Lya first
  x_pixminmax, wave, 1215.6701d, zhabs, -8000, 8000, PIXMIN=pmn, $
               PIXMAX=pmx, VELO=velo

  xr = [4256.5, 4263]
  plot, wave[pmn:pmx], fx[pmn:pmx], xrange=xr, $
        yrange=[-0.05,1.1], xtickn=xspaces, xmargin=[10,30], $
        ymargin=[10,14], xtitl='Wavelength (Ang)', ytit='Normalized Flux', $
        charsize=csz2, psym=10, background=clr.white, color=clr.black, $
        xstyle=1, ystyle=1, thick=1, ytickinterval=ytint, $
        /nodata

  ;; PLOT
  oplot, [-10000., 10000.], [0.,0.], color=clr.gray, linestyle=3,thick=3
  oplot, [-10000., 10000.], [1.,1.], color=clr.gray, linestyle=2,thick=3
  
  ;lya = x_setline(1215.6701)
  ;lya.N = NHI
  ;lya.b = bHI
  ;lya.zabs = zhabs
  ;vfx[pmn:pmx] = x_voigt(wave[pmn:pmx], lya, FWHM=FWHM)

  ;; Data
  oplot, wave[pmn:pmx], fx[pmn:pmx], color=clr.black, thick=5, psym=10

  ;; Label
  xyouts, 4259.5, 0.25, 'HI Ly!9a!X', color=clr.blue, charsize=LSIZE, align=0.5
  xyouts, 4258.15, 0.85, 'DI Ly!9a!X', color=clr.darkgreen, charsize=LSIZE, align=0.2

  plotsym, 1, 2.5, thick=6
  oplot, [4258.15], [0.8], psym=8, color=clr.darkgreen

  ;oplot, wave[pmn:pmx],vfx[pmn:pmx]*conti[pmn:pmx], color=clr.cyan, thick=3

  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]
  
  return
end
