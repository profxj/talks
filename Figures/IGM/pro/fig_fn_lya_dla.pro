pro fig_fn_lya_dla, infil, PSFILE=psfile

  ;; Get structure if necessary
  if not keyword_set( PSFILE ) then psfile = 'fig_fn_lya_dla.ps'
  if not keyword_set( DATFIL ) then datfil = 'fig_fnall.dat'

  if not keyword_set( NMIN ) then nmin = 12.0
  if not keyword_set( NMAX ) then nmax = 22.0

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


  x_psopen, psfile, /portrait
  !p.multi=[0,1,1,0,0]
  clr = getcolor(/load)
  lclr = clr.white

  lbl = 'Comb: N!dmin!N='+string(nmin,format='(f4.1)')
  XTIT='log N!dHI!N'

  yrng = [-26.,-10]
  xrng = [nmin, nmax]
  csz = 2.2
  ymrg = [3.5,0.5]
  plot, [0], [0], color=lclr, $
    background=clr.white, charsize=csz,$
    xmargin=[8,1.2], ymargin=ymrg, xtitle=XTIT, $
    ytitle='log f(N!dHI!N, X)', yrange=yrng, thick=4, $
    xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata, xtickn=xtck

  ;; Lya forest
  oplot, lyabin, log_lyafn, color=clr.red, psym=2
  writecol, 'dum', lyabin, log_lyafn, filnum=11

  ;; DLA
  if not keyword_set( GZFIL ) then gzfil = '~/SDSS/DR5_QSO/dr5_dlagz_s2n4.fits'
  sdss_fndla, GZFIL=GZFIL, STRCT=fnstr,  STP=stp, $
              /SALL, ALLDR=alldr, CL=0.95, NPLT=nplt, /NODBL
  ;; Shortcut
  xval = fnstr.xval
  yplt = fnstr.yplt
  bins = fnstr.bins
  yerr1 = fnstr.yerr1
  yerr2 = fnstr.yerr2

  val = where(yplt GT -99., complement=lim, ncomplement=nlim)
  oploterror, [xval[val]], [yplt[val]], [xval[val]-bins[val]],  $
              [yplt[val]-yerr2[val]], psym=3, thick=6, $
              color=clr.cyan, /lobar, errcolor=clr.cyan
  oploterror, [xval[val]], [yplt[val]], [bins[val]+stp-xval[val]],  $
              [yerr1[val]-yplt[val]], psym=3, thick=6, $
              color=clr.cyan, /hibar, errcolor=clr.cyan

  writecol, 'dum', xval[val], yplt[val], (yerr1[val]-yerr2[val])/2., filnum=11

  ;; Upper limits
  if nlim NE 0 then begin
      plotsym,1, 3., thick=5
      oploterror, [xval[lim]], [yerr1[lim]], [xval[lim]-bins[lim]], [0.], $
                  color=clr.cyan, psym=8, /lobar, errcolor=clr.cyan, thick=6
      oploterror, [xval[lim]], [yerr1[lim]], [bins[lim]+stp-xval[lim]], [0.], $
                  color=clr.cyan, psym=8, /hibar, errcolor=clr.cyan, thick=6
  endif

  ;; LLS
  rws_parm = [-0.677, -0.227, 21.63]  ;; dashed-dotted line in Fig 9 of O'Meara et al.
  NHIval = 16 + 3*findgen(1001)/1000.
  dNHIval = 10.d^NHIval - 10.d^shift(NHIval,1)
  dNHIval[0] = dNHIval[1]

  ;; RWS
  logON = rws_parm[0] + rws_parm[1]*(NHIval - 19.) - 10.^(NHIval - rws_parm[2])
  logfN  = logON - NHIval - alog10(alog(10.))
;  oplot, NHIval, logfN, color=clr.black, linesty=2

  ;; 15 to 16 spline
  xval = [lyabin, NHIval[0:5]]
  yval = [log_lyafn, logfN[0:5]]

  splin = spl_init(xval, yval, /doub)
  xeval = 15. + findgen(101)/100.
  yeval = spl_interp(xval, yval, splin, xeval,/doub)
  
;  oplot, xeval, yeval, color=clr.black, linesty=1
   
  ;; Close Ps
  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]
  close, /all
  
  return

end
