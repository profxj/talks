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
pro fig_1558_dh, datfil, NTOT=NTOT, CSIZE=csize, PSFILE=psfile, XTINT=xtint, $
              FINE=fine, INFIL=infil, LSIZE=lsize, YPS=yps, LTHICK=lthick

; lowzovi_prsdat -- Reads in DLA data to a structure

;if (N_params() LT 0) then begin 
;    print,'Syntax - ' + $
;      'fuse_velplt, strct, instr_list, vel_fil, NTOT=, CSIZE=, ' + $
;      'LSIZE=, PSFILE=, XTINT= (v1.1)' 
;    return
;endif 
;  
  if not keyword_set(PSFILE) then psfile = 'fig_1558_dh.ps'
  if not keyword_set( NTOT ) then ntot = 16L
  if not keyword_set( LSIZE ) then lsize = 1.5
  if not keyword_set( LTHICK ) then lthick = 4.
  if not keyword_set( INFIL ) then infil = 'Input/fig_1558_dh.inp'
  if not keyword_set( BLSZ ) then blsz = 1.4
  if not keyword_set( CSIZE ) then csize = 2.2
  if not keyword_set( CSZ2 ) then csz2 = 1.2
  lybfil = '../Data/J1558-0031a_94.fits'
  lybfit = '../Data/fit_lybHIRES.idl'

;  if not keyword_set( PSFIL ) then psfil = 'tmp.ps'

;  Read instrument file list
;readcol, instr_list, instr_fil, inst_dw, inst_w0, format='a,f,f'



  ;; Open vel_fil
  close, /all
  openr, 11, infil
  datfil = ''
  readf, 11, datfil
  readf, 11, FWHM
  readf, 11, zhabs, FORMAT='(f12.8)'
  readf, 11, zdabs, FORMAT='(f12.8)'
  readf, 11, NHI, FORMAT='(f5.2)'
  readf, 11, sigNHI, FORMAT='(f5.2)'
  readf, 11, bHI, FORMAT='(f)'
  readf, 11, Dval, FORMAT='(d)'
  readf, 11, sigD, FORMAT='(d)'
  readf, 11, bDH, FORMAT='(f)'

; PLOT
;if nlin LE ntot then begin

  fx = x_readspec(mikefil, wav=wave, NPIX=npix) 
  vfx = replicate(1.,npix)
  ufx = replicate(1.,npix)
  lfx = replicate(1.,npix)

  if keyword_set( PSFILE ) then x_psopen, psfile, /portrait
  clr = getcolor(/load)
  sclr = clr.gold

  ;; Lyb next
  fx = x_readspec(lybfil, wav=wave, NPIX=npix, inflg=2) 
  vfx = replicate(1.,npix)
  x_pixminmax, wave, 1025.7223d, zhabs, -1500, 1500, PIXMIN=pmn, $
               PIXMAX=pmx, VELO=velo

  plot, wave[pmn:pmx], fx[pmn:pmx], xrange=[3785, 3810], /noerase, $
        yrange=[-0.025,0.64], xtickn=xspaces, xmargin=[8,2], $
        ymargin=[0,0], /NODATA, ytitl='Intensity', $ ;ytickn=yspaces, $
        charsize=csz2, psym=10, background=clr.white, color=clr.black, $
        xstyle=1, ystyle=1, thick=lthick, ytickinterval=0.1, $
        pos=[0.1,0.24,0.88,0.80], xtitle='Wavelength (A)'
  
  lyb = x_setline(1025.7223)
  lyb.N = NHI
  lyb.b = bHI
  lyb.zabs = zhabs

  restore, lybfit
  vfx[pmn:pmx] = x_voigt(wave[pmn:pmx], lyb, FWHM=FWHM)

  ;; Shaded region
  lyb.N = NHI + sigNHI
  ufx[pmn:pmx] = x_voigt(wave[pmn:pmx], lyb, FWHM=FWHM)
  x_curvefill, wave[pmn:pmx], vfx[pmn:pmx]*conti[pmn:pmx], $
               ufx[pmn:pmx]*conti[pmn:pmx], color=sclr
  lyb.N = NHI - sigNHI
  lfx[pmn:pmx] = x_voigt(wave[pmn:pmx], lyb, FWHM=FWHM)
  x_curvefill, wave[pmn:pmx], vfx[pmn:pmx]*conti[pmn:pmx], $
               lfx[pmn:pmx]*conti[pmn:pmx], color=sclr

  ;; Data
  oplot, wave[pmn:pmx], fx[pmn:pmx], color=clr.black, thick=1, psym=10

  ;; Label
  xyouts, 3786., 0.52, 'HI 1025', color=clr.black, charsize=LSIZE

  ;; PLOT
  oplot, [-10000., 10000.], [0.,0.], color=clr.gray, linestyle=3,thick=1
  oplot, wave, conti, color=clr.purple, linestyle=2, thick=2
  oplot, wave[pmn:pmx],vfx[pmn:pmx]*conti[pmn:pmx], color=clr.green, thick=2

  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]

  return

  ;;starting at p.multi[0], plot p.multi[1] rows and p.multi[2] columns
  ;;with p.multi[3] z dimensions and going from top to bottom (p.multi[4]
  npy = 7
  npx = 2

  ;; Data
  fx = x_readspec(datfil, wav=wave, NPIX=npix) 
  vfx = replicate(1.,npix)

  ;; Vel strct
  wrest = 0.d
  dwrest = 0.d
  svinst = -1  ;;because SiC 1B is binary 1, which makes instr list # = 0
  nblnd = 0
  flg = 0
      
  ;; Lyman series
  readf, 11, vmin,vmax
  vmnx = [vmin,vmax]
  nlin = 0
  readf, 11, nlin
  
  ;; LOOP
  for ii=0L,nlin-1 do begin
      if ii LT 3 then offy = (npy+4) else offy = (npy)
      !p.multi=[offy-ii, npx,npy,0,1]
;      readf, 11, wrest, ymin, ymax, nblnd
          readf, 11, wrest, flg
          readf, 11, dwrest
          readf, 11, ymin, ymax
          ymnx = [ymin, ymax]
          
          readf, 11, nblnd
          
          
          ;; Set vel array
          x_pixminmax, wave, wrest, zhabs, vmnx[0], vmnx[1], PIXMIN=pmn, $
                       PIXMAX=pmx, VELO=velo
          
          ;; Y-axis
          if ymnx[1]-ymnx[0] GT 0.8 then ytint = 0.5 else begin
              if ymnx[1]-ymnx[0] GT 0.4 then ytint = 0.2 else ytint = 0.1
          endelse
          
          ;; Plot
          ysty=1
          if ((ii+1) MOD (nlin/2)) NE 0 then begin
              xspaces = replicate(' ',30) 
          endif else begin 
              if keyword_set(XSPACES) then delvarx, xspaces
          endelse
          if (ii GT npy-1) then yspaces = replicate(' ',30) 
          
          plot, velo[pmn:pmx], fx[pmn:pmx], xrange=vmnx, $
                yrange=ymnx, xtickn=xspaces, xmargin=[6,0], $
                ymargin=[0,0], NODATA=nblnd, $ ;ytickn=yspaces, $
                charsize=csize, psym=10, background=clr.white, color=clr.black, $
                xstyle=1, ystyle=ysty, thick=lthick, ytickinterval=ytint

          ;; Lines
          oplot, [0., 0.], ymnx, color=clr.blue, linestyle=2, thick=1
          oplot, [-82., -82.], ymnx, color=clr.gray, linestyle=2, thick=1
          oplot, [-10000., 10000.], [0.,0.], color=clr.gray, linestyle=3,thick=1
          oplot, [-10000., 10000.], [1.,1.], color=clr.gray, linestyle=1,thick=1
          
          ;; Labels
          getfnam, wrest, fv, nam 
          nam = strtrim(nam,2)
          if keyword_set(FINE) then begin
          abslin = x_setline(wrest)
          ;; SiII
          if strmid(nam,0,4) EQ 'SiII' and abslin.j EQ 0 then abslin.j = 0.5
          ;; CII
          if strmid(nam,0,3) EQ 'CII' and abslin.j EQ 0 then abslin.j = 0.5
          ;; OI
          if strmid(nam,0,2) EQ 'OI' and abslin.j EQ 0 $
            and abs(wrest-1306.0286) GT 0.1 then abslin.j = 2.0

          ;; FeII
          if strmid(nam,0,4) EQ 'FeII' then begin
              if abslin.j EQ 0 then abslin.j = 4.5
              nam2 = 'J='+strtrim(round(abslin.j*2),2)+'/2' 
          endif else begin
              if round(abslin.j) NE fix(abslin.j) then $
                nam2 = 'J='+strtrim(round(abslin.j*2),2)+'/2' $
              else nam2 = 'J='+strtrim(round(abslin.j),2)
          endelse
      endif

      case flg of 
          0: xyouts, 0.05*(vmnx[1]-vmnx[0])+vmnx[0], $
                     ymnx[0]+ (ymnx[1]-ymnx[0])*0.13, $
                     strtrim(nam,2), color=clr.black, charsize=LSIZE
          1: xyouts, 0.2*(vmnx[1]-vmnx[0])+vmnx[0], $
                     ymnx[0]+ (ymnx[1]-ymnx[0])*0.68, $
                     strtrim(nam,2), color=clr.black, charsize=LSIZE
          2: xyouts, 0.05*(vmnx[1]-vmnx[0])+vmnx[0], $
                     ymnx[0]+ (ymnx[1]-ymnx[0])*0.61, $
                     strtrim(nam,2), color=clr.black, charsize=LSIZE
          3: xyouts, 0.70*(vmnx[1]-vmnx[0])+vmnx[0], $
                     ymnx[0]+ (ymnx[1]-ymnx[0])*0.31, $
                     strtrim(nam,2), color=clr.black, charsize=LSIZE
          4: xyouts, 0.70*(vmnx[1]-vmnx[0])+vmnx[0], $
                     ymnx[0]+ (ymnx[1]-ymnx[0])*0.14, $
                     strtrim(nam,2), color=clr.black, charsize=LSIZE
          5: xyouts, 0.40*(vmnx[1]-vmnx[0])+vmnx[0], $
                     ymnx[0]+ (ymnx[1]-ymnx[0])*0.64, $
                     strtrim(nam,2), color=clr.black, charsize=LSIZE
          else: stop
      endcase
      
      fitlin = [x_setline(wrest), x_setline(dwrest)]
      fitlin[0].N = NHI
      fitlin[0].b = bHI
      fitlin[1].b = bDH
      fitlin.zabs = zhabs

      fitlin[1].N = Dval
      fitlin[1].zabs = zdabs
      vfx[pmn:pmx] = x_voigt(wave[pmn:pmx], fitlin, FWHM=FWHM)

      ;; Shaded region
      fitlin[1].N = Dval - sigD
      ufx[pmn:pmx] = x_voigt(wave[pmn:pmx], fitlin, FWHM=FWHM)
      x_curvefill, velo[pmn:pmx], vfx[pmn:pmx], $
                   ufx[pmn:pmx], color=sclr

      fitlin[1].N = Dval + sigD
      lfx[pmn:pmx] = x_voigt(wave[pmn:pmx], fitlin, FWHM=FWHM)
      x_curvefill, velo[pmn:pmx], vfx[pmn:pmx], $
                   lfx[pmn:pmx], color=sclr

      ;; Data
      if NBLND EQ 0 then $
        oplot, velo[pmn:pmx], fx[pmn:pmx], color=clr.black, $
               thick=lthick, psym=10

      ;; BLENDS
      nlow = pmn
      for jj=0L,nblnd-1 do begin
          readf, 11, vmin, vmax, FORMAT='(f,f)'
          mn = min(abs(velo-vmin),pixmin)
          mn = min(abs(velo-vmax),pixmax)
          
          ;; Plot good
          nlow = nlow < pixmin
          oplot, velo[nlow:pixmin], fx[nlow:pixmin], $
                 color=clr.black, psym=10, thick=lthick
          ;; Plot blend
          oplot, velo[pixmin:pixmax], fx[pixmin:pixmax], color=clr.orange, $
                 psym=10, linestyle=2, thick=1.
          ;; End
          if jj EQ nblnd-1 then begin
              if pixmax LT pmx then $
                oplot, velo[pixmax:pmx], fx[pixmax:pmx], color=clr.black, $
                       psym=10, thick=lthick
          endif
          ;; Set nlow
          nlow = pixmax
      endfor
        

      ;; Best Fit
      oplot, velo[pmn:pmx],vfx[pmn:pmx], color=clr.green, thick=1

  endfor 

  xyouts, 0.03, 0.80, 'Relative Flux', $
          alignment=0.5, ORIENTATION=90., /normal, charsize=blsz
  xyouts, 0.03, 0.35, 'Normalized Flux', $
          alignment=0.5, ORIENTATION=90., /normal, charsize=blsz
  yps = 0.08
  xyouts, 0.57, yps, 'Relative Velocity  (km s!u-1!N)', $
          alignment=0.5, /normal, charsize=blsz

  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]
  
  return
end
