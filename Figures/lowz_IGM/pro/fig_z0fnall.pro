pro fig_z0fnall, infil, PSFILE=psfile

  ;; Get structure if necessary
  if not keyword_set( PSFILE ) then psfile = 'fig_z0fnall.ps'
  if not keyword_set( DATFIL ) then datfil = 'fig_fnall.dat'

  if not keyword_set( NMIN ) then nmin = 12.0
  if not keyword_set( NMAX ) then nmax = 22.2

  x_psopen, psfile, /maxs
  !p.multi=[0,1,1,0,0]
  clr = getcolor(/load)

  cc = clr.black
  lbl = 'Comb: N!dmin!N='+string(nmin,format='(f4.1)')
  XTIT='log N!dHI!N'

  yrng = [-31.,-9]
  xrng = [nmin, nmax]
  csz = 2.2
  ymrg = [3.5,0.5]
  plot, [0], [0], color=clr.black, $
    background=clr.white, charsize=csz,$
    xmargin=[8,1.2], ymargin=ymrg, xtitle=XTIT, $
    ytitle='log f(N!dHI!N, X)', yrange=yrng, thick=4, $
    xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata, xtickn=xtck

  ;; Penton et al. 2004
  beta = -1.65
  norm = 10.3 

  readcol, 'penton_fN.dat', HI_clm, fN, sigHI, sigfN

;  dxdz = abs(cosm_xz(0.01,/wmap06) -  cosm_xz(0.05,/wmap06))/0.04
  dxdz = 1.  ;; This has a 3% error (low)

;  nplt = 100L
;  HI_clm = 12.3 + 2.2*findgen(nplt)/(nplt-1)
;  f_z = 10.^norm * (10.^HI_clm)^beta
;  f_x = f_z / dxdz

  ;; Lya
  lya = where(HI_clm LT 14.5)
  oploterror, HI_clm[lya], fN[lya], $;sigHI[lya], $
              sigfN[lya], $
              color=clr.red, errcolor=clr.red, thick=6, psym=3

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
              color=clr.red, errcolor=clr.red, thick=6, psym=3, errstyle=1

  lsz = 2.
  xyouts, 14.5, -12., 'Ly!9a!X Forest', color=clr.red, charsize=lsz

  ;; DLA (Zwaan et al. 2005)
  
  readcol, getenv('PSDSS')+'/DR3/Figures/Data/zwaan_fn.dat', $
           HI_clm, fN, sigfN
  gd = where(HI_clm GT 19.8 and HI_clm LE 22., ngd)
  oploterror, HI_clm[gd], fn[gd], $;replicate(0.05, ngd), $
              sigfn[gd], color=clr.blue, errcolor=clr.blue, thick=6, $
              psym=3
  xyouts, 20.5, -20., '21cm', color=clr.blue, charsize=lsz

  ;; H2
  nplt = 100L
  HI_clm = 22.0 + 2.0*dindgen(nplt)/(nplt-1)
  fstar = 1.1d-25
  sigma = 0.65
  mu = 20.6
  fN = fstar * exp(-0.5*((alog10(10.d^HI_clm)-mu)/sigma)^2)
  
;  oplot, HI_clm, alog10(fN), color=clr.black, linestyle=1
;  xyouts, 23.0, -26., 'CO', color=clr.black, charsize=lsz, alignment=0.5
;  printcol, HI_clm, alog10(fN)
  
  


  ;; Close Ps
  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]
  close, /all
  
  return

end
