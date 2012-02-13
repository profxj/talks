;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro cgm_baryfrac, summ_fil, ALL_GAL=all_gal

  if not keyword_set(psfile) then psfile = 'cgm_baryfrac.ps'
  if not keyword_set(binsz) then binsz = 0.01
  if not keyword_set(lsz) then lsz = 2.2 
  if not keyword_set(csz) then csz = 2.2 
  if not keyword_set(dlim) then dlim = 5. ;arcmin
  ;; Blanton lum function
  if not keyword_set(phi_str) then phi_str = 1.49  ; 10^-2 h^3 Mpc^-3
  if not keyword_set(alpha) then alpha = -1.05
  if not keyword_set(Mstar) then Mstar = -20.44  ; M* - 5 log h

  ;; Initialize
  c = x_constants()

  ;; Plot
  if keyword_set(psfile) then x_psopen,psfile,/maxs
  !p.multi=[0,1,1]
  clr = getcolor(/load)
;  lclr = clr.white
  lclr = clr.black

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; dN/dX
  yrng = [1e8, 1e13]
  xrng = [0.001, 10] 
  xmrg = [8,1.0]
  ylbl = -1

  
;  xspaces = replicate(' ',30) 
  plot,[0.],[0.],xrange=xrng,yrange=yrng,color=lclr,$
       background=clr.white, charsize=csz, $
       xmargin=xmrg, ymargin=[4.5,2],xtitle='Luminosity (!8L*!X)', $
       ytitle='Mass (M!dSun!N)',/nodata,xstyle=1,ystyle=1, /xlog, /ylog, $ ;, xtickn=xspaces
       xtickformat='x_logticks'
  
  L_val = xrng[0] * 10.^(alog10(xrng[1]/xrng[0]) * findgen(101)/100)

  ;; DM Halo
  rvir_dumb = 250. * (L_val)^0.2  ;; Prochaska et al. 2011 kludge
  Mvir_dumb = 200. * c.rhoc * 4 * !pi * (rvir_dumb*c.kpc)^3 / 3. / c.msun

  oplot, L_val, Mvir_dumb, color=lclr
  oplot, L_val, Mvir_dumb*0.17, color=lclr, linesty=1

  ;; CGM
  cgm_NH = [19.d, 20] ;; Total H
  ncgm = n_elements(cgm_NH)
  cgm_clr = [clr.red, clr.blue]
  M_300 = c.mp * 10.^cgm_NH * !pi * (300.*c.kpc)^2 * 1.3 / c.msun

  for ii=0L,ncgm-1 do begin
     ;; 300kpc
     oplot, xrng, replicate(M_300[ii],2), color=cgm_clr[ii], linesty=2

     ;; rvir
     M_cgm = c.mp * 10.^cgm_NH[ii] * !pi * (rvir_dumb*c.kpc)^2 * 1.3 / c.msun
     oplot, L_val, M_cgm, color=cgm_clr[ii], linesty=3
  endfor

  ;; Label
  xlbl = 0.7
  xyouts, xlbl, 1e9, 'DM Halo', color=lclr, charsiz=lsz
  xyouts, xlbl, 1e9/2., 'N!dH!N = 10!u20!N cm!u-2!N', color=clr.blue, charsiz=lsz
  xyouts, xlbl, 1e9/4., 'N!dH!N = 10!u19!N cm!u-2!N', color=clr.red, charsiz=lsz


  x_psclose
  !p.multi=[0,1,1]

  return
end
