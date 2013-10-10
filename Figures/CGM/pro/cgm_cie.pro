;;   fig_cldycoll, [2e4, 2e6], psfile='Figures/fig_cldycoll.ps'
pro cgm_cie, tlim, CLDY=cldy, MTL=mtl, NHI=NHI, ION_FIL=ion_fil, $
                CMPION=cmpion, PSFILE=psfile, CHK=chk

  ;; Get structure if necessary
  if not keyword_set( TLIM ) then tlim = [1.7e4, 2e6]
  if not keyword_set( LSZ ) then lsz = 1.9
  if not keyword_set( ymnx ) then ymnx=[13., 17.3]
  if not keyword_set( CLDY ) then $
    cldy='/u/xavier/idl/xidl/Cloudy/gnat_CIE.fits'

  ;; Open
  coll = xmrdfits(cldy,1)

; PLOT
  psfile = 'cgm_cie.ps'
  if keyword_set( PSFILE ) then x_psopen, psfile, /maxs
  clr = getcolor(/load)
  fclr = clr.white


  xrng = [1e5, 1e7]
  xmarg=[8,15]
  yrng = [0., 1]
  plot, [0], [0], color=fclr, background=clr.white, charsize=2.2,$
        xmargin=xmarg, ymargin=[5,1], xtitle='T (K)', ytitle='Ionic fraction', $
        xrang=xrng, yrang=yrng, $
        /nodata, xstyle=1, ystyle=1, /xlog


  ;; Let's just do some Oxygen
  oplot, coll.T, coll.O[4], color=fclr, thick=5, linesty=1
  oplot, coll.T, coll.O[5], color=fclr, thick=5 
  oplot, coll.T, coll.O[6], color=fclr, thick=5, linesty=2

  ;; Label
  xyouts, 2e5, 0.55, 'OV', color=clr.cyan, charsi=lsz
  xyouts, 3.5e5, 0.15, 'OVI', color=clr.cyan, charsi=lsz
  xyouts, 8.5e5, 0.85, 'OVII', color=clr.cyan, charsi=lsz

  if keyword_set( PSFILE ) then x_psclose

  return
end
      
      
