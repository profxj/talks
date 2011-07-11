pro fig_fnlya_lowvshighz, infil, PSFILE=psfile

  ;; Get structure if necessary
  if not keyword_set( PSFILE ) then psfile = 'fig_fnlya_lowvshighz.ps'
  if not keyword_set( DATFIL ) then datfil = 'fig_fnall.dat'

  if not keyword_set( NMIN ) then nmin = 12.0
  if not keyword_set( NMAX ) then nmax = 22.0
  if not keyword_set( LSZ ) then lsz = 1.8

  if not keyword_set( GZFIL ) then $
    gzfil = '~/SDSS/DR3_QSO/dr3_dlagz_s2n4.fits'


  lyabin =  findgen(10)*0.36+11.70
  kt_spots = [489., 498., 499., 485., 468., 447., 419., 383., 333., 315.]
  log_lyafn = -10. + 4. * (kt_spots-512.)/(512.-316) 
  log_lyao3n = log_lyafn + lyabin + alog10(alog(10.)) $
                -alog10(sqrt(3.7/0.3)) + alog10(1.2)
  log_lyafn = alog10(10.^log_lyafn * (0.201/0.361) )   ; WMAP06 vs SCDM with H=100
  printcol, lyabin, log_lyafn

  close, /all
  openw, 11, datfil

  ;; Plot

  x_psopen, psfile, /portrait
  !p.multi=[0,1,1,0,0]
  clr = getcolor(/load)
  lclr = clr.white

  lbl = 'Comb: N!dmin!N='+string(nmin,format='(f4.1)')
  XTIT='log N!dHI!N'

  yrng = [-17.5,-9.5]
  xrng = [11.5, 17.]
  csz = 2.2
  ymrg = [3.5,0.5]
  plot, [0], [0], color=lclr, $
    background=clr.white, charsize=csz,$
    xmargin=[8,1.2], ymargin=ymrg, xtitle=XTIT, $
    ytitle='log f(N!dHI!N, X)', yrange=yrng, thick=4, $
    xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata, xtickn=xtck

  ;; High z Lya forest
  oplot, lyabin, log_lyafn, color=clr.pink, psym=2, symsiz=1.5
  xplt = 12.
  oplot, [xplt], [-16.0], color=clr.pink, psym=2, symsiz=1.5
  xyouts, xplt+0.2, -16.1, 'z~2.5; Kirkman+02', color=clr.pink, charsiz=lsz
;  writecol, 'dum', lyabin, log_lyafn, filnum=11

  ;; Low z
  readcol, '../lowz_IGM/penton_fN.dat', HI_clm, fN, sigHI, sigfN

  lya = where(HI_clm LT 14.5)
  oploterror, HI_clm[lya], fN[lya], $;sigHI[lya], $
              sigfN[lya], $
              color=clr.cyan, errcolor=clr.cyan, thick=6, psym=1

  ;; Penton (14.5-17.5)
  beta = -1.33
  norm = 5.2

;  nplt = 100L
;  HI_clm = 14.5 + 3.0*findgen(nplt)/(nplt-1)
;  f_z = 10.^norm * (10.^HI_clm)^beta
;  f_x = f_z / dxdz
;  oplot, HI_clm, alog10(f_x), color=clr.black, linestyle=1

  high = where(HI_clm GT 14.5)
  oploterror, HI_clm[high], fN[high], sigHI[high], sigfN[high], $
              color=clr.cyan, errcolor=clr.cyan, thick=6, psym=3, errstyle=1
  oplot, [xplt], [-16.5], color=clr.cyan, psym=1, symsiz=1.5
  xyouts, xplt+0.2, -16.6, 'z~0; Penton+04', color=clr.cyan, charsiz=lsz

  ;; Close Ps
  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]
  close, /all
  
  return

end
