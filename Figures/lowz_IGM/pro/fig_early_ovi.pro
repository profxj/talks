;fig_1671 
pro fig_early_ovi, PSFILE=psfile, XR=xr, YR=yr

  ;; keywords
  if not keyword_set( LSZ ) then lsz = 2.5
  if not keyword_set( LBLS ) then lbls = 1.5
  if not keyword_set( XR ) then xr = [3501., 6800.]
  if not keyword_set( ROOT ) then $
    root = '/u/xavier/LCO/OVI/FUSE/data/PKS1302-102/'
  if not keyword_set(PSFILE) then psfile='fig_early_ovi.ps'
                                           
  ;; PLOT
  if keyword_set( PSFILE ) then x_psopen, psfile, /maxs
  clr = getcolor(/load)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Early
  wfccd_wrfspec, galstr, root+'GSpec/PKS1302-102_wfccd.fits', /read

  yr = [0., 0.8]
  plot, [0.], [0.], xrange=xr, yrange=yr, color=clr.black, $
    background=clr.white, charsize=1.5, thick=7, $
    xmargin=[9,0.5], ymargin=[4.5,0.1], ytitle='F!d!9l!N!X (10!u-16!n cgs)', $
    xtitle='Rest Wavelength ('+string("305B)+')',  /nodata, xstyle=1, ystyle=1

  xyouts, 3985., 0.4, 'CaH+K', orientation=90, charsize=lbls, color=clr.black
  xyouts, 4325., 0.4, 'Gband', orientation=90, charsize=lbls, color=clr.black
  xyouts, 6564.7, 0.53, 'H!9a!x', orientation=90, charsize=lbls, color=clr.black
  xyouts, 4882.7, 0.53, 'H!9b!x', orientation=90, charsize=lbls, color=clr.black
  xyouts,5096.,0.5,'[OIII] 5006',orientation=90, charsize=lbls, color=clr.black
  xyouts, 4132., 0.55, 'H!9d!x', orientation=90, charsize=lbls, color=clr.black
  xyouts, 4400., 0.55, 'H!9g!x', orientation=90, charsize=lbls, color=clr.black
  xyouts, 3748., 0.43, '[OII]', orientation=90, charsize=lbls, color=clr.black
  xyouts, 5200., 0.5, 'MgI 5175', orientation=90,  charsize=lbls, color=clr.black
  xyouts, 5920., 0.5, 'NaI 5894', orientation=90, charsize=lbls, color=clr.black

  idx = where(galstr.slit_id EQ 2033L)

  ;; bad pix
;  bd = where(galstr[idx].var LE 0 and galstr[idx].wave LT 4668 and $
;             galstr[idx].wave GT 0.,nbd)
;  if nbd NE 0 then galstr[idx].fx[bd] = galstr[idx].fx[bd[0]-1]

  ;; Data
  oplot, galstr[idx].wave / (galstr[idx].zans.z+1.), $
    galstr[idx].fx*1e16, color=clr.black, thick=7

  oplot, galstr[idx].wave / (galstr[idx].zans.z+1.), $
    sqrt(galstr[idx].var)*1e16, color=clr.red, thick=4

  ;; Label
  xyouts, 0.14, 0.90, 'PKS1302-102: ID 2033', color=clr.black, charsize=lsz, /norm
  xyouts, 0.14, 0.85, $
    'z = '+string(galstr[idx].zans.z,format='(f6.4)'), $
    color=clr.black, charsize=lsz, /norm

  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1,1]

  return
end
      
      
