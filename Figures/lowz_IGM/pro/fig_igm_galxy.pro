pro fig_igm_galxy, infil, PSFILE=psfile

  ;; Get structure if necessary
  if not keyword_set( PSFILE ) then psfile = 'fig_igm_galxy.ps'
  if not keyword_set( DATFIL ) then datfil = 'fig_fnall.dat'

  if not keyword_set( NMIN ) then nmin = 12.0
  if not keyword_set( NMAX ) then nmax = 22.2

  x_psopen, psfile, /maxs
  !p.multi=[0,1,1,0,0]
  clr = getcolor(/load)
  lclr = clr.white

  lbl = 'Comb: N!dmin!N='+string(nmin,format='(f4.1)')
  XTIT='log N!dHI!N'

  yrng = [-30.,-9]
  xrng = [nmin, nmax]
  if not keyword_set(WIDE) then ymrg = [3.5,3.5] else ymrg = [6.5, 5]
  csz = 2.2

  for qq=0,4 do begin
     plot, [0], [0], color=lclr, $
           background=clr.white, charsize=csz,$
           xmargin=[8,1.2], ymargin=ymrg, xtitle=XTIT, $
           ytitle='log f(N!dHI!N, X)', yrange=yrng, thick=4, $
           xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata, xtickn=xtck

     ;; Penton et al. 2004
     beta = -1.65
     norm = 10.3 
     
     readcol, 'penton_fN.dat', HI_clm, fN, sigHI, sigfN

;  dxdz = abs(cosm_xz(0.01,/wmap06) -  cosm_xz(0.05,/wmap06))/0.04
     dxdz = 1. ;; This has a 3% error (low)

;  nplt = 100L
;  HI_clm = 12.3 + 2.2*findgen(nplt)/(nplt-1)
;  f_z = 10.^norm * (10.^HI_clm)^beta
;  f_x = f_z / dxdz

     ;; Lya
     lya = where(HI_clm LT 14.5)
     oploterror, HI_clm[lya], fN[lya], $ ;sigHI[lya], $
                 sigfN[lya], $
                 color=clr.cyan, errcolor=clr.cyan, thick=6, psym=3

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
     
     lsz = 2.
     xyouts, 12.3, -14.5, 'Ly!9a!X Forest', color=clr.cyan, charsize=lsz

     if qq GT 0 and qq LT 4 then begin
        oplot, replicate(14.5,2), yrng, color=lclr, linesty=2
        xyouts, 13.3, -28., 'IGM', color=clr.cyan, charsiz=lsz, align=0.5
     endif

     ;; LLS
     xyouts, 18.5, -18., 'LLS', color=clr.orange, charsize=lsz, align=0.5

     if qq GT 1 and qq Lt 4 then begin
        oplot, replicate(19.8,2), yrng, color=lclr, linesty=2
        xyouts, 17.15, -28., 'CGM', color=clr.orange, charsiz=lsz, align=0.5
     endif
     
     ;; DLA (Zwaan et al. 2005)
     readcol, getenv('PSDSS')+'/DR3/Figures/Data/zwaan_fn.dat', $
              HI_clm, fN, sigfN, /sile
     gd = where(HI_clm GT 19.8 and HI_clm LE 22., ngd)
     oploterror, HI_clm[gd], fn[gd], $ ;replicate(0.05, ngd), $
                 sigfn[gd], color=clr.yellow, errcolor=clr.yellow, thick=6, $
                 psym=3
     xyouts, 20.0, -20.5, 'DLA (21cm)', color=clr.yellow, charsize=lsz
     
     if qq GT 2 then begin
        oplot, replicate(19.8,2), yrng, color=lclr, linesty=2
        xyouts, 21.0, -28., 'ISM', color=clr.yellow, charsiz=lsz, align=0.5
     endif
  
;  xyouts, 12.5, -25.5, 'Penton+04', color=lclr, charsize=lsz, align=0
;  xyouts, 12.5, -27., 'Lehner+07', color=lclr, charsize=lsz, align=0
;  xyouts, 12.5, -28.5, 'Zwaan+05', color=lclr, charsize=lsz, align=0

  endfor


  ;; Close Ps
  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]
  close, /all
  
  return

end
