;+ 
; NAME:
; fig_velplt
;  V1.1
;
; PURPOSE:
;    Velcoity plot of FUSE transitions for a given absorption 
;    system
; CALLING SEQUENCE:
;   
;   fig_velplt, strct_fil, instr_list, vel_fil, NTOT=, CSIZE=,
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
;
; OPTIONAL OUTPUTS:
;  PSFILE -- Postscript filename
;
; COMMENTS:
;
; EXAMPLES:
;   fig_calcewn, struct, fil_instr
;
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   12-Sep-2003 Written by JXP
;-
;------------------------------------------------------------------------------
pro fig_q2000, dat_fil, NTOT=NTOT, CSIZE=csize, $
                LSIZE=lsize, PSFILE=psfile, XTINT=xtint, $
                SCREEN=screen, TWO_PANEL=two_panel

; lowzovi_prsdat -- Reads in DLA data to a structure

  if not keyword_set(CFILB) then $
    cfilb = '/u/xavier/LCO/MIKE/RedData/Q2000-330/Q2000_conti_b.fits'
  if not keyword_set(CFILR) then $
    cfilr = '/u/xavier/LCO/MIKE/RedData/Q2000-330/Q2000_conti_r.fits'
  if not keyword_set(DATFILB) then $
    datfilb = '/u/xavier/LCO/MIKE/RedData/Q2000-330/Q2000-330a_b_F.fits'
  if not keyword_set(DATFILR) then $
    datfilr = '/u/xavier/LCO/MIKE/RedData/Q2000-330/Q2000-330a_r_F.fits'
;  if not keyword_set( DAT_FIL ) then $
;    dat_fil = '/u/xavier/Keck/HIRES/RedData/TiII/HD195965/HD195965_f.fits'
  if not keyword_set( NTOT ) then ntot = 16L
  if not keyword_set( LSIZE ) then lsize = 1.6
  if not keyword_set( LSZ ) then LSZ = 2.4
  if not keyword_set( LTHICK ) then lthick = 3
  if not keyword_set( PSFILE ) and not keyword_set( SCREEN ) $
    then psfile = 'fig_q2000.ps'
  if not keyword_set( CSIZE ) then csize = 2.4
  LSIZE = 2.0

  ;; Read data
  fxb = x_readspec(datfilb, wav=waveb, sig=sigb, NPIX=npix)
  contib = xmrdfits(cfilb)
  fxr = x_readspec(datfilr, wav=waver, sig=sigr, NPIX=npix)
  contir = xmrdfits(cfilr)

  wvcut = 4775.
  sclr = 4./3 * 0.985
  gdb = where(waveb LT wvcut)
  gdr = where(waver GT wvcut)
  fx = [fxb[gdb], fxr[gdr]*sclr]
  wave = [waveb[gdb], waver[gdr]]
  conti = [contib[gdb], contir[gdr]*sclr]

; Open vel_fil
  close, /all

  xmrg = [6,0.5]
  ymrg = [3.5,0.5]
  xr = [4200., 5900]
  yr = [0., 7.5]
  xspaces = replicate(' ',30) 

  ;; PSFILE
  if keyword_set( PSFILE ) then x_psopen, psfile, /maxs
  !p.multi=[0,1,2,0,1]
  clr = getcolor(/load)

; LOOP

  for ss=1,0,-1 do begin
      ;; 
      !p.multi=[2-ss,1,2,0,1]
      plot, wave, fx, xrange=xr, yrange=yr, $
            xmargin=xmrg, ymargin=ymrg, /NODATA, $
            charsize=csize, psym=10, background=clr.white, $
            xstyle=1, ystyle=1, thick=5, $ ;ytitle='Brightness', $
            color=clr.black, xtickn=xspaces, ytickn=xspaces
  
      if ss EQ 1 then oplot, wave, smooth(fx,3), color=clr.blu5, psym=10, thick=3 $
      else oplot, wave, conti, color=clr.blu5, thick=7

      ;; Label
      asz = 3.
      xv = 4255.
      if ss EQ 1 then lbl = 'Cosmic Web' else lbl= 'No Cosmic Web'
          ;; LLS
      xyouts, xv, 6.3, lbl, color=clr.black, charsiz=lsz, alignment=0.

  endfor
  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]

  return
end
