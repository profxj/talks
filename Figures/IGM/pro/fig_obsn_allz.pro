;; Plotting f(N,X)
;; Normalize to O_m = 0.3, O_L = 0.7
pro fig_obsn_allz, infil, PSFILE=psfile

  ;; Get structure if necessary
  if not keyword_set( PSFILE ) then psfile = 'fig_obsn_allz.ps'
  if not keyword_set( DATFIL ) then datfil = 'fig_fnall.dat'

  if not keyword_set( NMIN ) then nmin = 12.0
  if not keyword_set( NMAX ) then nmax = 22.0
  if not keyword_set( LSZ ) then lsz = 1.8

  if not keyword_set( GZFIL ) then $
    gzfil = '~/SDSS/DR3_QSO/dr3_dlagz_s2n4.fits'


  ;; Kirkman & Tytler 1997
  ;;  Not sure the dX is for 0.3, 0.7 cosmology [JXP, 15-07-2011]
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
  lclr = clr.black
  z0clr = clr.red
  z3clr = clr.blue
  z4clr = clr.tomato

  lbl = 'Comb: N!dmin!N='+string(nmin,format='(f4.1)')
  XTIT='log N!dHI!N'
  YTIT='O(N!dHI!N)'

  yrng = [0.0001, 1000]
  xrng = [12, 22.]
  csz = 2.0
  ymrg = [3.5,0.5]

  for qq=0,1 do begin

  plot, [0], [0], color=lclr, $
    background=clr.white, charsize=csz,$
    xmargin=[8,1.2], ymargin=ymrg, xtitle=XTIT, $
    ytitle=YTIT, yrange=yrng, thick=4, $
    xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata, xtickn=xtck, /ylog

  ;; High z Lya forest
  gdp = where(lyabin GT 12.4)
  a10 = alog(10.)
  oplot, lyabin[gdp], 10.^(log_lyafn[gdp]+lyabin[gdp])*a10, $
         color=z3clr, psym=2, symsiz=1.5
  xplt = 12.
;  oplot, [xplt], [-16.0], color=z2clr, psym=2, symsiz=1.5
;  xyouts, xplt+0.2, -16.1, 'z~2.5; Kirkman+02', color=clr.pink, charsiz=lsz
;  writecol, 'dum', lyabin, log_lyafn, filnum=11

  ;; Low z
  readcol, '../lowz_IGM/penton_fN.dat', HI_clm, fN, sigHI, sigfN

  lya = where(HI_clm LT 14.5)
  if qq GE 1 then begin
     oploterror, HI_clm[lya], 10.d^(HI_clm[lya]+fN[lya])*a10, $ ;sigHI[lya], $
                 sigfN[lya]*10.^(HI_clm[lya]+fN[lya]) * a10, $
                 color=z0clr, errcolor=z0clr, thick=6, psym=1

  ;; Penton (14.5-17.5)
  beta = -1.33
  norm = 5.2

;  nplt = 100L
;  HI_clm = 14.5 + 3.0*findgen(nplt)/(nplt-1)
;  f_z = 10.^norm * (10.^HI_clm)^beta
;  f_x = f_z / dxdz
;  oplot, HI_clm, alog10(f_x), color=clr.black, linestyle=1

  high = where(HI_clm GT 14.5)
  oploterror, HI_clm[high], 10.d^(fN[high]+HI_clm[high])*a10, sigHI[high], $
              sigfN[high]*10.d^(HI_clm[high]+fn[high])*a10, $
              color=z0clr, errcolor=z0clr, thick=6, psym=3, errstyle=1
;  oplot, [xplt], [-16.5], color=z0clr, psym=1, symsiz=1.5
;  xyouts, xplt+0.2, -16.6, 'z~0; Penton+04', color=z0clr, charsiz=lsz
  endif


  ;; ;;;;;;;;;;;
  ;; SLLS (z=3)
  N_SLLS = [19.15, 19.45, 19.75, 20.1]
  bin_SLLS = [0.15, 0.15, 0.15, 0.2]
  fN_SLLS = [-20.17, -20.7, -21.04, -21.48]
  sigfN_SLLS = [0.1, 0.14, 0.15, 0.15]
  oploterror, N_SLLS, 10.d^(fN_SLLS+N_SLLS)*a10, bin_SLLS, $
              sigfN_SLLS*10.^(N_SLLS+fN_SLLS)*a10, $
              color=z3clr, errcol=z3clr, thick=6, psym=1

  ;; ;;;;
  ;; LLS (beta=1)
  N_LLS = 17.0 + findgen(10)*2./9.
  fN_LLS = fN_SLLS[0] - 1.*(N_LLS-19.15)
  oplot, N_LLS, 10.d^(fN_LLS+N_LLS)*a10, color=z3clr, linesty=2

  ;; ;;;;
  ;; pLLS (inferred)
  xval= [N_LLS[0], max(lyabin)]
  oplot, xval, 10.d^([fN_LLS[0], min(log_lyafn)] + xval)*a10, $
         color=z3clr, linesty=1

  ;; ;;;;;;;;;;;
  ;; DLAs (<z> = 3)
  if qq EQ 0 then begin
     if not keyword_set( GZFIL ) then gzfil = '~/SDSS/DR5_QSO/dr5_dlagz_s2n4.fits'
     sdss_fndla, GZFIL=GZFIL, STRCT=fnstr,  STP=stp, $
                 /SALL, ALLDR=alldr, CL=0.95, NPLT=nplt, /NODBL
  endif
  ;; Shortcut
  xval = fnstr.xval
  yplt = 10.^(fnstr.yplt + xval)*a10
  bins = fnstr.bins
  yerr1 = 10.^(fnstr.yerr1 + xval)*a10
  yerr2 = 10.d^(fnstr.yerr2 + xval)*a10

  val = where(yplt GT -99., complement=lim, ncomplement=nlim)
  oploterror, [xval[val]], [yplt[val]], [xval[val]-bins[val]],  $
              [yplt[val]-yerr2[val]], psym=3, thick=6, $
              color=z3clr, /lobar, errcolor=z3clr
  oploterror, [xval[val]], [yplt[val]], [bins[val]+stp-xval[val]],  $
              [yerr1[val]-yplt[val]], psym=3, thick=6, $
              color=z3clr, /hibar, errcolor=z3clr

  writecol, 'dum', xval[val], yplt[val], (yerr1[val]-yerr2[val])/2., filnum=11

  ;; Upper limits
  if nlim NE 0 then begin
     plotsym,1, 3., thick=5
     oploterror, [xval[lim]], [yerr1[lim]], [xval[lim]-bins[lim]], [0.], $
                 color=z3clr, psym=8, /lobar, errcolor=z3clr, thick=6
     oploterror, [xval[lim]], [yerr1[lim]], [bins[lim]+stp-xval[lim]], [0.], $
                 color=z3clr, psym=8, /hibar, errcolor=z3clr, thick=6
  endif
  
  ;; DLA (Zwaan et al. 2005)
  if qq GE 1 then begin
     readcol, getenv('PSDSS')+'/DR3/Figures/Data/zwaan_fn.dat', $
              HI_clm, fN, sigfN, /sile
     gd = where(HI_clm GT 19.8 and HI_clm LE 22., ngd)
     oploterror, HI_clm[gd], 10.^(fn[gd]+HI_clm[gd])*a10, $ ;replicate(0.05, ngd), $
                 sigfn[gd]* 10.d^(fn[gd]+HI_clm[gd])*a10, $
                 color=z0clr, errcolor=z0clr, thick=6, $
                 psym=3
     
     xyouts, 18., -11., 'z~0', color=z0clr, charsiz=lsz
  endif

  xyouts, 18., -12., 'z~3', color=z3clr, charsiz=lsz

endfor

  ;; Close Ps
  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]
  close, /all
  
  return

end
