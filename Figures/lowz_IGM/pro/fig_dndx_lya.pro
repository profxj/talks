;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro fig_dndx_lya, summ_fil, ALL_GAL=all_gal

  if not keyword_set(psfile) then psfile = 'fig_dndx_lya.ps'
  if not keyword_set(binsz) then binsz = 0.01
  if not keyword_set(lsz) then lsz = 1.5 
  if not keyword_set(csz) then csz = 2.0 
  if not keyword_set(dlim) then dlim = 5. ;arcmin
  ;; Blanton lum function
  if not keyword_set(phi_str) then phi_str = 1.49  ; 10^-2 h^3 Mpc^-3
  if not keyword_set(alpha) then alpha = -1.05
  if not keyword_set(Mstar) then Mstar = -20.44  ; M* - 5 log h

  ;; Initialize
  root = '~/paper/OVI/Galaxies/'
  cd, root+'/Analysis/pro/', curr=curr
  galabs_initparm, init
  cd, curr

  if not keyword_set(hubb) then hubb = init.H0 / 100.

  ;; Phi(M) = 0.4 log(10) Phi* 10^(-0.4 [M-M*][alpha+1]) exp(-
                    ; 10^(-0.4[M-M*] ) )
  ;; Phi(L) = Phi* (L/L*)^alpha exp(-L/L*)

  c = x_constants()
  phi_str_cgs = (phi_str * 1e-2 * hubb^3) / (c.mpc^3)
  dndx_const = c.c / (hubb * 100 * 1e5 / c.Mpc) 


  ;; Plot
  if keyword_set(psfile) then x_psopen,psfile,/maxs
  !p.multi=[0,1,1]
  clr = getcolor(/load)
  lclr = clr.black

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; dN/dX
  yrng = [0.1, 300]
  xrng = [0.001, 10] 
  xmrg = [8,1.0]
  ylbl = -1
  
;  xspaces = replicate(' ',30) 
  plot,[0.],[0.],xrange=xrng,yrange=yrng,color=lclr,$
       background=clr.white, charsize=csz, $
       xmargin=xmrg, ymargin=[4.5,2],xtitle='Luminosity (!8L*!X)', $
       ytitle='',/nodata,xstyle=1,ystyle=1, /xlog, /ylog, $ ;, xtickn=xspaces
       xtickformat='x_logticks'

  xlbl = 0.05
  xyouts, 0.05, 0.35, 'Cumulative', alignment=0.0, color=lclr,$
    charsize=2.3, /normal, orientation=90
  !p.font = -1
  xyouts, xlbl, 0.6, '!12l!X', alignment=0.0, color=lclr,$
    charsize=2.8, /normal, orientation=90
  !p.font = 0
  xyouts, xlbl, 0.61, '!dLy!9a!X!N(X)', alignment=0.0, color=lclr,$
    charsize=2.4, /normal, orientation=90

  ;; Shade out faint-end
  x_curvefill, [0.001, 0.01], replicate(yrng[1],2), replicate(yrng[0],2), $
               color=clr.gray, /line_fill, orient=135.
  
 
  ;; Observations
  ;; (Penton et al. 2004)  Assumed to be z=0.05
  xv = cosm_xz([0.03, 0.08])
  dx = abs(xv[0]-xv[1])/0.05
  xlbl3 = 2.

  ;; 300mA  [CHECK!]
  x_curvefill, xrng, [15.0, 15.0]/dx, [9., 9.]/dx, color=clr.green, /line_fill, orient=45.
  xyouts, xlbl3, 15./dx, 'W!uLy!9a!X!N > 300 mA', color=lclr, charsi=lsz

  ;; 100mA
  x_curvefill, xrng, [70.0, 70.0]/dx, [50., 50.]/dx, color=clr.orange, /line_fill, orient=45.
  xyouts, xlbl3, 70./dx, 'W!uLy!9a!X!N > 100 mA', color=lclr, charsi=lsz

  ;; 30mA
  x_curvefill, xrng, [180.0, 180.0]/dx, [160., 160.]/dx, color=clr.red, /line_fill, orient=45.
  xyouts, xlbl3, 180./dx, 'W!uLy!9a!X!N > 30 mA', color=lclr, charsi=lsz

  ;;;;;;;;;;;;;;;
  ;;  r_p = r_vir
  ;;  r_vir = 250 kpc * (L/L*)^beta   ;; Using beta = 0.2

  Lval = xrng[0] * 10.^(alog10(xrng[1]/xrng[0]) * findgen(101)/100.)
  rvir_0 = init.rvir[2] * c.kpc ; kpc
  beta = init.beta
  x = alpha + 1 + beta
  dndx_rvir = dndx_const * phi_str_cgs * (!pi * rvir_0^2) * $
              gamma(x) * ( igamma(x,xrng[1]) - igamma(x,Lval))

  oplot, Lval, dndx_rvir, color=lclr

  ;;;;;;;;;;;;;;;
  ;;  r_p = r_cgm
  ;;  r_cgm = 300 kpc 

  rcgm_0 = init.rcgm * c.kpc ; kpc
  beta = 0.
  x = alpha + 1 + beta
  dndx_rcgm = dndx_const * phi_str_cgs * (!pi * rcgm_0^2) * $
              gamma(x) * ( igamma(x,xrng[1]) - igamma(x,Lval))

  oplot, Lval, dndx_rcgm, color=clr.blue, linesty=1

;  x1 = alpha + 1 
;  dndx_subls1 = dndx_const * phi_str_cgs * (!pi * (100.*c.kpc)^2) * $
;               gamma(x1) * ( igamma(x1,1.) - igamma(x1,0.1))
;  dndx_subls2 = dndx_const * phi_str_cgs * (!pi * rvir_0^2) * $
;              gamma(x) * ( igamma(x,1.) - igamma(x,0.1))
;  print, 'SubL*', dndx_subls1, dndx_subls2

  xlbl2 = [0.02, 0.03]
  ylbl2 = 0.3
  oplot, xlbl2, replicate(ylbl2,2), color=clr.blue, linest=1
  xyouts, xlbl2[1]*1.3, ylbl2*0.9, 'A!dp!N = !9p!X r!S!dCGM!R!N!u2!N', $
          color=clr.blue, charsiz=lsz

  ylbl2 = ylbl2/1.5
  oplot, xlbl2, replicate(ylbl2,2), color=lclr
  xyouts, xlbl2[1]*1.3, ylbl2*0.9, 'A!dp!N = !9p!X r!S!dvir!R!N!u2!N', $
          color=lclr, charsiz=lsz
              
  x_psclose
  !p.multi=[0,1,1]

  return
end
