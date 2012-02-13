pro fig_grb_excite, conti, J0_NH2=J0_NH2, J1_NH2=J1_NH2, AV=AV, $
            NOSMC=nosmc, NOPS=nops, SUB=sub, NOFIT=nofit, $
            NOERR=noerr, NOLY=noly, REVIEW=review, NOLBL=nolbl

  ;; Get structure if necessary
  if not keyword_set( CSZ ) then csz = 1.4
  if not keyword_set( LSZ ) then lsz = 1.5
  if not keyword_set( YSTP ) then ystp = 0.
  if not keyword_set( SMTH ) then SMTH = 1
  if not keyword_set(NOPS) then psfile = 'fig_grb_excite.ps'
  if not keyword_set(CONTI) then conti = 0.5
  if not keyword_set(AV) then AV = 0.2
  if not keyword_set(NOSMC) then SMC = 1
  if not keyword_set(C0) then c0 = 1.6 ;; lambda = 1340Ang
  if not keyword_set(NHISIG) then nhisig = 0.15

  npix = 10000L
  wave = 1400. + 0.05 * dindgen(npix)
  fx = findgen(npix)
 
  ;; Unexcited lines
  atomic_lines = [x_setline([1526.7066, 1548.195, 1550.770, 1608.4511]), $
                  x_setline([1526.7066, 1548.195, 1550.770, 1608.4511])] 
  atomic_lines[4:*].zabs  = 25. / 3e5
  atomic_lines.b  = 8. 
  atomic_lines.N  = 16.

  CI_lines = x_setline([1560.309])
  CI_lines.b = 5.
  CI_lines.N = 13.
  atomic_lines = [atomic_lines, CI_lines]
  fx_atoms = x_voigt(wave, atomic_lines, FWHM=3.)

  ;; CO lines
  Tex = 10. ; K
  N_CO = 16.5
  colines = read_colines()

  COwv = [1322.153, 1344.183, 1367.623, 1392.525, 1419.044d, 1477.565, 1509.748, 1544.302]
  COlbl = [$
          'CO 8-0', $
          'CO 7-0', $
          'CO 6-0', $
          'CO 5-0', $
          'CO 4-0', $
          'CO 2-0', $
          'CO 1-0', $
          'CO 0-0']
  c=x_constants()
  xrng = [1460., 1650]

  gdCO = where( colines.c_isotope EQ 12 AND $
                colines.wrest GT xrng[0] and colines.wrest LT xrng[1], nco) 
;              print, 'Parsing ', np, nco, tex
              
              ;; Set up lines
  colin = replicate({newabslinstrct}, nco)
  colin.f = colines[gdCO].f
  colin.gamma = colines[gdCO].gamma
  colin.wrest = colines[gdCO].wrest
  colin.zabs = 0.
  colin.b = 2.0

  ;; Set columns
  gJ = 2*colines[gdCO].Jpp + 1
  u_Jpp = colines[gdCO[uniq(colines[gdCO].Jpp, $
                            sort(colines[gdCO].Jpp))]].Jpp
  nJpp = n_elements(u_Jpp)
  dE = dblarr(nJpp)
  for ff=0L,nJpp-1 do $
     dE[ff] = abs(co_xlevel(0,u_Jpp[ff])-co_xlevel(0,0)) ; cm^-1
  dE = dE * c.h * c.c  ;; ergs
  nJ = gJ * exp(-1. * dE / c.k / Tex)
  norm = 10.^N_CO / total(nJ)
  nJ = nJ * norm
  for ff=0L,nJpp-1 do begin
     idx = where(colines[gdCO].Jpp EQ u_Jpp[ff], nidx)
     if nidx EQ 0 then stop
     colin[idx].N = alog10( nJ[ff] )
  endfor
  fx_CO = x_voigt(wave, colin, FWHM=4.)
  
  ;; PLOT
  if keyword_set( PSFILE ) then begin
     x_psopen, psfile, /maxs
      !p.multi=[0,1,2]
   endif
  
  clr = getcolor(/load)
  c = x_constants()

  ;; Idealized, no excitation
  yrng = [-0.05, 1.05]
  ymrg = [4.,0.5]
  plot, [0], [0],  color=clr.black, background=clr.white, charsize=csz,$
        xmargin=[8,2], ymargin=ymrg, pos=pos, $
        xrange=xrng, yrange=yrng, $
        ytitle='Normalized Flux', $ 
        xtitle='Rest Wavelength (Ang)', $ 
        /nodata, xstyle=1, ystyle=1

  oplot, wave, fx_CO*fx_atoms, color=clr.black, psym=10, thick=4

  ;; Label CO
  zgrb = 0.
  xpt = [replicate(1477.565,2), replicate(1509.748,3), replicate(1544.302,2)]
  hrz = -0.03
  ypt = [0.1, hrz, hrz, 0.1, hrz, hrz, 0.1]
  oplot, xpt*(1+zgrb), ypt, color=clr.darkgreen
  xyouts, 1480., 0, 'CO (T!dex!N=10K)', color=clr.darkgreen, charsize=lsz, alignment=0.0
  xyouts, 1530., 0.1, 'SiII', color=clr.red, charsize=lsz, alignment=0.0
  xyouts, 1610., 0.1, 'FeII', color=clr.red, charsize=lsz, alignment=0.0

  ;; Actual

  r400 = '~/GRB/data/080607/LRIS/GRB080607_R400N.fits'
  fx = x_readspec(r400, wav=wv, sig=sig, inflg=2)

  plot, [0], [0],  color=clr.black, background=clr.white, charsize=csz,$
        xmargin=[8,2], ymargin=ymrg, pos=pos, $
        xrange=xrng, yrange=yrng, $
        ytitle='Normalized Flux', $ 
        xtitle='Rest Wavelength (Ang)', $ 
        /nodata, xstyle=1, ystyle=1

  oplot, wv/(1+3.03626), fx, color=clr.black, psym=10, thick=4

  ;; Label CO
  zgrb = 0.
  xpt = [replicate(1477.565,2), replicate(1509.748,3), replicate(1544.302,2)]
  hrz = -0.03
  ypt = [0.1, hrz, hrz, 0.1, hrz, hrz, 0.1]
  oplot, xpt*(1+zgrb), ypt, color=clr.darkgreen
  xyouts, 1480., 0, 'CO (T!dex!N=10K)', color=clr.darkgreen, charsize=lsz, alignment=0.0
  xyouts, 1530., 0.1, 'SiII', color=clr.red, charsize=lsz, alignment=0.0
  xyouts, 1610., 0.1, 'FeII', color=clr.red, charsize=lsz, alignment=0.0

  
  ;; H2
  ;xyouts, 3900., 0.70, 'H!d2!N Lyman-Werner', color=clr.blue, charsiz=2.
  ;oplot, [912., 912., 1115., 1115]*(1+zgrb), 0.2+[0.3,0.4,0.4,0.3], $
  ;       color=clr.blue, thick=2
  
  if keyword_set( PSFILE ) then begin
      x_psclose
      !p.multi=[0,1,1]
  endif

  return
end
      
      
