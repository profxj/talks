pro fig_h2_080607, conti, J0_NH2=J0_NH2, J1_NH2=J1_NH2, AV=AV, $
                   NOSMC=nosmc, NOPS=nops, SUB=sub, NOFIT=nofit, $
                   NOERR=noerr, NOLY=noly, REVIEW=review

  ;; Get structure if necessary
  if not keyword_set( CSZ ) then csz = 1.9
  if not keyword_set( LSZ ) then lsz = 1.1
  if not keyword_set( YSTP ) then ystp = 0.
  if not keyword_set( SMTH ) then SMTH = 1
  if not keyword_set(NOPS) then psfile = 'fig_h2_080607.ps'
  if not keyword_set(CONTI) then conti = 0.5
  if not keyword_set(AV) then AV = 0.2
  if not keyword_set(NOSMC) then SMC = 1
  if not keyword_set(C0) then c0 = 1.6 ;; lambda = 1340Ang
  if not keyword_set(NHISIG) then nhisig = 0.15

  b600 = '~/GRB/data/080607/LRIS/GRB080607_B600X.fits'
  lyafit = '~/GRB/data/080607/LRIS/GRB080607_lya.idl'
  restore, lyafit
  lyalines = lines

  model = { $
        Tex: 0., $
        b: 0., $
        fwhm: 0., $
        NH2: 0. }
  if not keyword_set(DEF_TEX) then def_Tex = 200.  ;; K
  if not keyword_set(DEF_B) then def_b = 3. ;; km/s
  if not keyword_set(DEF_NH2) then def_NH2 = 21.2 ;; log
  if not keyword_set(DEF_FWHM) then def_FWHM = 5. ;; log
  model.NH2 = DEF_NH2
  model.b = DEF_b
  model.fwhm = DEF_FWHM
  model.Tex = DEF_TEX

  ;; Dust values
  AV = 3.206
  R_V = 4.008
  c3 = 1.33
  c4 = 0.28

  ;; PLOT
  if keyword_set( PSFILE ) then begin
      x_psopen, psfile, /maxs
      !p.multi=[0,1,1]
  endif
  clr = getcolor(/load)
  c = x_constants()

  zgrb = 3.03626

  xrest = [950., 1400]
  xrng = xrest*(1+zgrb)
  yrng = [-0.1, 1.1]
  ymrg = [4,1]

  ;; Loop
;  xtit = 'Rest Wavelength (Ang)'
  xtit = 'Wavelength (Ang)'
  fx = x_readspec(b600, wav=wv, inflg=2)
  rwv = wv / (1+zgrb)
  rnu = c.c/rwv ;; Goofy units

  ;; Convert to fnu
  fx = fx * wv * (wv/1d8) / c.c * 1e11 ;; erg/s/cm^2/Hz 1e23

  ;; Correct flux for Galactic extinction
  AV_MW = 0.07
  Alambda_MW = AV_MW * x_fitzpalav(wv)
  fx = fx * 10.^(0.4*Alambda_MW)

  ;; Now correct N_HI
  lya_abs = x_voigt(rwv, lyalines, FWHM=6)

  ;; Dust
  Al = AV * x_fitzpalav(rwv, R_V=R_V, c3=c3, c4=c4)

  ;; Tie to 7140Ang (~1750Ang rest)
  pix = where(wv GT 5326.8 and wv LT 5331.75)
  medfx = median(fx[pix])

  ;; Take out Lya opacity
  medfx = medfx / (lya_abs[round(mean(pix))])

  medwv = median(wv[pix]) / (1+zgrb)
  mednu = c.c/medwv

  Almed = AV * x_fitzpalav(medwv, R_V=R_V, c3=c3, c4=c4)

  trufx0 = medfx * 10.^(0.4*Almed)
  trufx = (rnu/mednu)^(-0.5) * trufx0  ;; Only good for fnu fluxes

  ;; Reddened continuum
  red_conti = trufx * 10.^(-0.4*Al)

  ;; Smooth
  fx = smooth(fx,SMTH)

  ;; Plot
  pos = [0.08, 0.45, 0.98, 0.98]
  plot, [0], [0],  color=clr.lightgray, background=clr.black, charsize=csz,$
        xmargin=[8,2], ymargin=ymrg, xtitle=xtit, pos=pos, $
        xrange=xrng, yrange=yrng, xtickn=xspaces, $
        ytitle='Normalized Flux', $ 
;        ytitle='Relative Flux (erg/s/Ang/cm!u2!N)', $
        /nodata, xstyle=1, ystyle=1

  oplot, xrng, [0., 0.], color=clr.orange, thick=5, linesty=2

  ;; Turn back to f_lambda
  conv = 1. / wv / (wv/1d8) * c.c / 1e11 ;; erg/s/cm^2/Hz 1e23
  flam = fx * conv
  red_conti =  red_conti*conv

  ;; Continuum
;  oplot, wv, conti, color=clr.blue, thick=1, linesty=1
;  oplot, wv, red_conti, color=clr.orange, thick=5, linesty=2

  N_H2 = model.NH2
  
  nJpp = 14L
  Jval = lindgen(nJpp)
  gJ = 2*Jval + 1
  dE = dblarr(nJpp)
  for ff=0L,nJpp-1 do $
         dE[ff] = h2_xlevel(Jval[ff]) ; cm^-1
  
  dE = dE * c.h * c.c  ;; ergs
  nJ = gJ * exp(-1. * dE / c.k / model.Tex)
  norm = 10.^N_H2 / total(nJ)
  nJ = nJ * norm
  
  ;; Truncate
  keep = where(nJ GT 1e10, nJpp)
  maxJ = Jval[nJpp-1]
  gJ = gJ[keep]
  dE = dE[keep]
  print, 'Max J = ', maxJ
  
  nJ = gJ * exp(-1. * dE / c.k / model.Tex)
  norm = 10.^N_H2 / total(nJ)
  nJ = nJ * norm
  
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; H2
  h2lin = read_h2lin()
  gdh2 = where(h2lin.Jpp LE maxJ and $
               h2lin.wrest GT xrest[0])
  h2lin = h2lin[gdh2]
  
  for ff=0L,nJpp-1 do begin
      idx = where(h2lin.Jpp EQ ff, nidx)
      if nidx EQ 0 then stop
      h2lin[idx].N = alog10( nJ[ff] )
  endfor
;      print, 'NH2 = ', alog10(total(nJ))
  
  h2lin.zabs = 0.  ;; Rest-frame
  h2lin.b = model.b
  
  ;; Lyman
  lyman = x_setline([949.7431,972.5368, 1025.7223d,1215.6701d])
  nly = n_elements(lyman)
  lyman.N = lyalines.N
  lyman.b = 30.
  tmp = h2lin[0:nly-1]
  tmp.N = lyman.N
  tmp.b = 30.
  tmp.f = lyman.f
  tmp.wrest = lyman.wrest
  tmp.gamma = lyman.gamma
  lines = [h2lin,tmp]

  ;; Plot Lyman error
  if not keyword_set(NOERR) then begin
      lyman.N = lyalines.N - NHIsig
      modellow = x_voigt(rwv, lyman, FWHM=model.fwhm)
      lyman.N = lyalines.N + NHIsig
      modelhi = x_voigt(rwv, lyman, FWHM=model.fwhm)
      x_curvefill, wv, modellow, modelhi, color=clr.yellow
;      x_curvefill, wv, modellow*red_conti, modelhi*red_conti, color=clr.yellow
  endif
      
  ;; Data
  oplot, wv, flam/red_conti, color=clr.lightgray, psym=10, thick=5

  
  ;; Total model (plot)
  model = x_voigt(rwv, lines, FWHM=model.fwhm)
  oplot, wv, model, color=clr.cyan, thick=5
  
  ;; LBL for Review/Talks
  ;; H2
  xyouts, 3900., 1.00, 'H!d2!N Lyman-Werner', color=clr.green, charsiz=2.
  oplot, [912., 912., 1115., 1115]*(1+zgrb), 0.55+[0.3,0.4,0.4,0.3], $
         color=clr.green, thick=5
  ;; HI
  wlya = 1215.67 * (1+zgrb)
  xyouts, wlya, 0.95, 'HI Ly!9a!X', color=clr.green, charsiz=2., alignm=0.5
  plotsym, 1, 4., thick=4
  oplot, [wlya], [0.80], psym=8, color=clr.green


  ;; IGM
;  xyouts, mean([xrng[0],4800.]), 1.5, 'Ly!9a!X Forest', $
;          color=clr.blue, alignmen=0.5, charsiz=2.
;  oplot, [912., 912., 1190.0, 1190.0]*(1+zgrb), 1.0+[0.3,0.4,0.4,0.3], $
;         color=clr.blue, thick=2
  
;  if not keyword_set(NOLY) then begin
;      lyman.N = lyalines.N
;      modelLym = x_voigt(rwv, lyman, FWHM=model.fwhm)
;      oplot, wv, modelLym*red_conti, color=clr.blue, thick=2
;  endif


  if keyword_set( PSFILE ) then begin
      x_psclose
      !p.multi=[0,1,1]
  endif

  return
end
      
      
