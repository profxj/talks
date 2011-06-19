pro fig_z0_lofx, infil, PSFILE=psfile

  ;; Get structure if necessary
  if not keyword_set( PSFILE ) then psfile = 'fig_z0_lofx.ps'

  if not keyword_set( NMIN ) then nmin = 12.0
  if not keyword_set( NMAX ) then nmax = 22.2
  lsz = 2.

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

  for qq=0,3 do begin

     plot, [0], [0], color=clr.black, $
           background=clr.white, charsize=csz,$
           xmargin=[8,1.2], ymargin=ymrg, xtitle=XTIT, $
           ytitle='log f(N!dHI!N, X)', yrange=yrng, thick=4, $
           xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata, xtickn=xtck
     
     ;; Penton et al. 2004
     beta = -1.65
     norm = 10.3 
     
     readcol, 'penton_fN.dat', HI_clm, fN, sigHI, sigfN, /silen

     dxdz = 1. ;; This has a 3% error (low)

     ;; Lya
     lya = where(HI_clm LT 14.5, nlya)
     oploterror, HI_clm[lya], fN[lya], $ ;sigHI[lya], $
                 sigfN[lya], $
                 color=clr.red, errcolor=clr.red, thick=6, psym=3

     if qq GT 0 then begin
        x_curvefill, HI_clm[lya], replicate(yrng[0], nlya), fN[lya], $ 
                     color=clr.red 
        dN = 10.^HI_clm[lya] - shift(10.^HI_clm[lya],1)
        dN[0] = dN[1]
        lofX_lya = total(10.^fN[lya] * dN)
        xyouts, 13.25, -11., 'l(X) = '+string(lofX_lya,format='(i3)'), $
                color=clr.red, charsiz=lsz
     endif

     ;; Penton (14.5-17.5)
     beta = -1.33
     norm = 5.2

     high = where(HI_clm GT 14.5, nhigh)
     oploterror, HI_clm[high], fN[high], sigHI[high], sigfN[high], $
                 color=clr.darkgreen, errcolor=clr.darkgreen, thick=6, psym=3, errstyle=1
     if qq GT 1 then begin
        x_curvefill, HI_clm[high], replicate(yrng[0], nhigh), fN[high], $ 
                     color=clr.darkgreen 
        dN = 10.^HI_clm[high] - shift(10.^HI_clm[high],1)
        dN[0] = dN[1]
        lofX_high = total(10.^fN[high] * dN)
        xyouts, 16.0, -15., 'l(X) = '+string(lofX_high,format='(f3.1)'), $
                color=clr.darkgreen, charsiz=lsz
     endif

     if qq EQ 0 then xyouts, 14.5, -12., 'Ly!9a!X Forest', color=clr.red, charsize=lsz

     ;; DLA (Zwaan et al. 2005)
     
     readcol, getenv('PSDSS')+'/DR3/Figures/Data/zwaan_fn.dat', $
              HI_clm, fN, sigfN, /sile
     gd = where(HI_clm GT 19.8 and HI_clm LE 22., ngd)
     oploterror, HI_clm[gd], fn[gd], $ ;replicate(0.05, ngd), $
                 sigfn[gd], color=clr.blue, errcolor=clr.blue, thick=6, $
                 psym=3
     if qq GT 2 then begin
        x_curvefill, HI_clm[gd], replicate(yrng[0], ngd), fN[gd], $ 
                     color=clr.blue 
        dN = 10.^HI_clm[gd] - shift(10.^HI_clm[gd],1)
        dN[0] = dN[1]
        lofX_gd = total(10.^fN[gd] * dN)
        xyouts, 20.0, -21., 'l(X) = '+string(lofX_gd,format='(f4.2)'), $
                color=clr.blue, charsiz=lsz
     endif
     if qq EQ 0 then xyouts, 20.5, -20., '21cm', color=clr.blue, charsize=lsz

  endfor

  ;; Close Ps
  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]
  close, /all
  
  return

end
