pro fig_magestack, RREAL=rreal

  if not keyword_set( PSFILE ) then psfile = 'fig_magestack.ps'
  if not keyword_set(PAD_FRAC) then pad_frac = 0.1
  if not keyword_set(CSZ) then csz = 1.9
  if not keyword_set(CSZ2) then csz2 = 2.3
  if not keyword_set(lSZ) then lsz = 1.9
  if not keyword_set(lSZ2) then lsz2 = 1.7
  if not keyword_set(XNCOLORS) then xncolors=200L
  if not keyword_set(YSTP) then ystp = 0.07

  thk = 9
  xlbl = 0.3
  xlbl2 = 0.05
  ylbl = 0.90

  close, /all
  fx = x_readspec('../../../OCIW_2010/mage_stack.fits', inf=2, wav=wave)

  ;;; BEGIN PLOTS
  x_psopen, psfile, /maxs
  !p.multi = [0,1,1]
  clr = getcolor(/load)

  xmrg = [7,1]
  ymrg = [4.0,0.5]

  ;; Plot FeII
  yrng=[-0.5, 6.0]
  xrng=[830., 1400]
  xcut = 2605.
  off = 1.

  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, ytitle='Relative Flux', $
        xtitle='Rest Wavelength (Ang)', yrange=yrng, thick=4, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata, /noerase, $
        xthi=thk, ythi=thk
  
  oplot, wave, fx, color=clr.black, thick=5, psym=10
  xyouts, 900., 4., 'MagE Stack', color=clr.red, charsiz=lsz

  if keyword_set( PSFILE ) then x_psclose
  !p.multi = [0,1,1]
  close, /all

  return

end
