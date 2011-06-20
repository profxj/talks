;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro fig_dndx_ovi, summ_fil, ALL_GAL=all_gal

  if not keyword_set(psfile) then psfile = 'fig_dndx_ovi.ps'
  if not keyword_set(binsz) then binsz = 0.01
  if not keyword_set(lsz) then lsz = 1.9 
  if not keyword_set(csz) then csz = 2.0 
  if not keyword_set(dlim) then dlim = 5. ;arcmin
  ;; Blanton lum function
  if not keyword_set(phi_str) then phi_str = 1.49  ; 10^-2 h^3 Mpc^-3
  if not keyword_set(alpha) then alpha = -1.05
  if not keyword_set(Mstar) then Mstar = -20.44  ; M* - 5 log h

  ;; Initialize
  root = '~/paper/OVI/Galaxies/'
  cd, root+'Analysis/pro/', curr=curr
  galabs_initparm, init
  cd, curr

  if not keyword_set(hubb) then hubb = init.H0 / 100.

  c = x_constants()


  ;; Plot
  if keyword_set(psfile) then x_psopen,psfile,/maxs
  !p.multi=[0,1,1]
  clr = getcolor(/load)
  lclr = clr.white

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; dN/dX
  yrng = [0., 15]
  xrng = [0, 4] 
  xmrg = [8,0.5]
  ylbl = -1
  
  xspaces = replicate(' ',30) 
  plot,[0.],[0.],xrange=xrng,yrange=yrng,color=lclr,$
       background=clr.white, charsize=csz, xminor=1, $
       xmargin=xmrg, ymargin=[4.5,1],$;xtitle='Luminosity', $
       ytitle='',/nodata,xstyle=1,ystyle=1, xtickn=xspaces, thick=9

  !p.font = -1
  xyouts, 0.06, 0.5, '!12l!X', alignment=0.0, color=lclr,$
    charsize=2.8, /normal, orientation=90
  !p.font = 0
  xyouts, 0.06, 0.51, '!dOVI!N(X)', alignment=0.0, color=lclr,$
    charsize=2.4, /normal, orientation=90
 
  ;; Observations
  ;; (Tripp et al. 2008) -- Evaluated at z~0.2
  xv = cosm_xz([0.2, 0.25])
  dx = abs(xv[0]-xv[1])/0.05
  ;; 70mA
  x_curvefill, xrng, [10.9, 10.9]/dx, [7.1, 7.1]/dx, color=clr.darkgray, /line_fill, orient=45.
  xyouts, 0.5, 8.8/dx, 'W!u1031!N > 70 mA', color=clr.white, charsi=lsz
  ;; 100mA
  x_curvefill, xrng, [3.3, 3.3]/dx, [6.0, 6.0]/dx, color=clr.darkgray, /line_fill, orient=45.
  xyouts, 0.5, 4.5/dx, 'W!u1031!N > 100 mA', color=clr.white, charsi=lsz
  ;; 30mA
  x_curvefill, xrng, [18.5, 18.5]/dx, [13.2, 13.2]/dx, color=clr.darkgray, /line_fill, orient=45.
  xyouts, 0.5, 15.6/dx, 'W!u1031!N > 30 mA', color=clr.white, charsi=lsz

  ;;;;;;;;;;;;;;;;;;;;
  ;; Dwarfs
  mrang = -2.5 * alog10( [0.01, 0.1] )  ;; 5 to 2.5 mag fainter
  dm = 0.01
  nval = abs(mrang[1]-mrang[0])/dm
  mval = min(mrang) + dm * findgen(nval+1)

  phiM =  0.4 * alog(10.) * phi_str * 10.^(-0.4 * mval*(alpha+1)) * $
          exp(-10.^(-0.4*mval))
  ndwarf = total(phiM*dm) * 1e-2  * hubb^3  ; Mpc^-3
  ;; Check 
  x = alpha+1
  ndwarf2 = phi_str * gamma(x) * ( igamma(x,0.1) - igamma(x,0.01)) * hubb^3 * 1e-2
  nall = phi_str * gamma(x) * ( igamma(x,10.) - igamma(x,0.01)) * hubb^3 * 1e-2
  print, 'ndwarf = ', ndwarf, ndwarf2
  print, 'ngal (L>0.01L*) = ', nall
  dwarf_area = !pi * (init.rvir[0]*c.kpc)^2 ; Upper limit
  dwarf_dndx =  c.c / (hubb * 100 * 1e5 / c.Mpc) * (ndwarf/ c.Mpc^3) * dwarf_area

  x_curvefill, [0.75,1.25], [0., 0.], replicate(dwarf_dndx,2), color=clr.yellow
  xyouts, 1., ylbl, 'Dwarfs', color=clr.yellow, charsi=lsz, align=0.5

  ;;;;;;;;;;;;;;;;;;;;
  ;; sub-L*
  mrang = -2.5 * alog10( [0.1, 1.0] )  ;; 5 to 2.5 mag fainter
  dm = 0.01
  nval = abs(mrang[1]-mrang[0])/dm
  mval = min(mrang) + dm * findgen(nval+1)

  phiM =  0.4 * alog(10.) * phi_str * 10.^(-0.4 * mval*(alpha+1)) * $
          exp(-10.^(-0.4*mval))
  nsubls = total(phiM*dm) * 1e-2  * hubb^3  ; Mpc^-3
  nsubls2 = phi_str * gamma(x) * ( igamma(x,1) - igamma(x,0.1)) * hubb^3 * 1e-2
  print, 'nsubls = ', nsubls, nsubls2
  subls_area = !pi * (init.rvir[1]*c.kpc)^2 ; r_vir
  subls_area2 = !pi * (init.rcgm*c.kpc)^2 ;  CGM
  subls_dndx =  c.c / (hubb * 100 * 1e5 / c.Mpc) * (nsubls/ c.Mpc^3) * subls_area
  subls_dndx2 =  c.c / (hubb * 100 * 1e5 / c.Mpc) * (nsubls/ c.Mpc^3) * subls_area2

  x_curvefill, [1.75,2.25], [0., 0.], replicate(subls_dndx2,2), color=clr.cyan, /line_fil,$
               orient=-45.
  x_curvefill, [1.75,2.25], [0., 0.], replicate(subls_dndx,2), color=clr.cyan
  xyouts, 2., ylbl, 'Sub-L*', color=clr.cyan, charsi=lsz, align=0.5

  ;;;;;;;;;;;;;;;;;;;;
  ;; L*
  mrang = -2.5 * alog10( [1.0,10.] )  ;; 5 to 2.5 mag fainter
  dm = 0.01
  nval = abs(mrang[1]-mrang[0])/dm
  mval = min(mrang) + dm * findgen(nval+1)

  phiM =  0.4 * alog(10.) * phi_str * 10.^(-0.4 * mval*(alpha+1)) * $
          exp(-10.^(-0.4*mval))
  nlstar = total(phiM*dm) * 1e-2  * hubb^3  ; Mpc^-3
  nlstar2 = phi_str * gamma(x) * ( 1. - igamma(x,1.)) * hubb^3 * 1e-2
  print, 'nlstar = ', nlstar, nlstar2
  lstar_area = !pi * (init.rvir[2]*c.kpc)^2 ; Upper limit
  lstar_dndx =  c.c / (hubb * 100 * 1e5 / c.Mpc) * (nlstar/ c.Mpc^3) * lstar_area

  x_curvefill, [2.75,3.25], [0., 0.], replicate(lstar_dndx,2), color=clr.tomato
  xyouts, 3., ylbl, 'L*', color=clr.tomato, charsi=lsz, align=0.5

  x_psclose
  !p.multi=[0,1,1]

  return
end
