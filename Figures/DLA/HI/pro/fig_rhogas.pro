;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro fig_rhogas, PSFILE=psfile, GZFIL=gzfil, DLALST=dlalst, PLTRAO=pltrao, $
                 TOTON=toton, NOPRX=noprx, NO_PS=no_ps, $
                MODELS=models, STRCT=strct, BINS=bins, LBLSZ=lblsz, $
                CSZ=csz, XRNG=xrng, z0=z0, NDAEME=ndaeme, $
                DRK=drk

  ;; Get structure if necessary
  if not keyword_set( GZFIL ) then gzfil = '~/SDSS/DR5_QSO/dr5_dlagz_s2n4.fits'
  if not keyword_set( XZFIL ) then xzfil = $
     getenv('PSDSS')+'/DR5/Analysis/xz_val_L70_M30.fits'
  if not keyword_set( PSFILE ) and not keyword_set(NO_PS) then $
    psfile = 'fig_rhogas.ps'
  if not keyword_set(CSZ) then csz = 2.2

  c = x_constants()
  Hub = 72.  ; km/s/Mpc
  if not keyword_set(XRNG) then xrng = [2.0, 5.]
  yrng = [0.0, 1.5]

  ;; PLOT
  if keyword_set( PSFILE ) then x_psopen, psfile, /maxs
  clr = getcolor(/load)
  if keyword_set(DRK) then begin
     lclr = clr.white
     zclr = clr.yellow
  endif
     lclr = clr.black
     zclr = x_fsc_color('blu5')
  endelse

  if keyword_set(NO_PS) then begin
      xspaces = replicate(' ',30) 
      xtit='' 
      ymrg = [0.,0]
  endif else begin
      ymrg = [3.5,3]
      xtit='z' 
  endelse

  if not keyword_set( NOPLT ) then $
    plot, [0.], [0.], color=lclr, background=clr.white, charsize=csz,$
          xmargin=[8,4], ymargin=ymrg, xtitle=xtit, xtickn=xspaces, $
          ytitle='!9r!X!dHI!N (10!u8!N M!dSun!N Mpc!u-3!N h!d'+ $
          string(round(Hub),format='(i2)')+'!N)', $
          /nodata, xrange=xrng, ystyle=1, yrange=yrng, xstyle=1

;  xspaces = replicate(' ',30) 
;  axis, xaxis=1, charsize=csz, xstyle=1, xtickn=xspaces, xminor=1, xthick=5

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; TOTAL
;  bins = [ [1.7,2.2], [2.2, 2.5], [2.5, 3.0], [3.0, 3.5], [3.5,5.3]]
  if not keyword_set( BINS ) then $
    bins = [ [2.2, 2.4], [2.4,2.7], [2.7, 3.0], [3.0, 3.5], $
             [3.5,4.5]];, [4.,5.5]]
  if n_elements(bins) EQ 2 then sz = [2,1] else sz = size(bins,/dimensions)
  nbin = sz[1]

  ;; Zwaan
  oHI = 3.5e-4
  sigHI = 0.4e-4
  
  NCST = (c.rhoc * (Hub/100.) * (Hub/100.) )  / c.msun * (c.Mpc)^3 / 1e8 $
         * (75./Hub)
  rhoHI = oHI * NCST
  sig_rhoHI = sigHI * NCST

  x_curvefill, [0., 10.], replicate(rhoHI+sig_rhoHI,2), $
               replicate(rhoHI-sig_rhoHI,2), color=zclr, /line_fill, $
               orientation=45.
  lsz= 2.4
  xyouts, 4.6, 0.6, 'z=0', color=zclr, charsi=lsz


  ;; DR5
  if not keyword_set(z0) then begin
      sdss_omgdla, GZFIL=GZFIL, STRCT=STRCT, BINS=bins, /ALL, CST=cst, XZFIL=xzfil
      NCST = (hub*1e5/c.Mpc)*c.mp/c.c / CST / c.msun * (c.Mpc)^3
      strct.omega = strct.omega * NCST ;  Msun Mpc^-3 h
      strct.somg = strct.somg * NCST ;  Msun Mpc^-3 h

      if keyword_set(NDAEME) and strct[0].medz LT 2.4 then begin
          ;; Take mean of two values and error for the two
          orig = strct[0].omega
          strct[0].omega = mean(strct[0:1].omega)
          strct[0].somg[0] = strct[1].omega
          strct[0].somg[1] = strct[0].somg[1] + 0.05
      endif
      
      for ii=0L,nbin-1 do begin 
          oploterror, [strct[ii].medz], [strct[ii].omega/1e8], $
                      [strct[ii].medz-strct[ii].bins[0]], $
                      (strct[ii].omega-strct[ii].somg[1])/1e8, $
                      color=lclr, /lobar, errcolor=lclr, psym=psym
          oploterror, [strct[ii].medz], [strct[ii].omega/1e8], $
                      [strct[ii].bins[1]-strct[ii].medz], $
                      (strct[ii].omega-strct[ii].somg[0])/1e8, $
                      color=lclr, /hibar, errcolor=lclr, psym=psym
      endfor
      printcol, [strct.medz], [strct.omega/1e8], $
                [strct.medz-strct.bins[0]], $
                (strct.omega-strct.somg[1])/1e8 
  endif


  ;; Close Ps
  if keyword_set( PSFILE ) then x_psclose
  print, 'fig_rhogas: All done!'

  return
end
      
      
