pro fig_080607_r400, conti, AV=AV, NOSMC=nosmc, NOPS=nops, SUB=sub, REVIEW=review

  ;; Get structure if necessary
  if not keyword_set( CSZ ) then csz = 1.4
  if not keyword_set( LSZ ) then lsz = 1.1
  if not keyword_set( YSTP ) then ystp = 0.
  if not keyword_set( SMTH ) then SMTH = 1
  if not keyword_set(NOPS) then psfile = 'fig_r400.ps'
  if not keyword_set(CONTI) then conti = 0.5
  if not keyword_set(J0_NH2) then J0_NH2 = 21.0
  if not keyword_set(J1_NH2) then J1_NH2 = 21.0
  if not keyword_set(AV) then AV = 0.2
  if not keyword_set(NOSMC) then SMC = 1
  if not keyword_set(C0) then c0 = 1.6 ;; lambda = 1340Ang

  ;; Dust values
  AV = 3.206
  R_V = 4.008
  c3 = 1.33
  c4 = 0.28

  c = x_constants()

  r400 = '~/GRB/data/080607/LRIS/GRB080607_R400X.fits'
  cr400 = '~/GRB/data/080607/LRIS/GRB080607_R400C.fits'

  ;; PLOT
  if keyword_set( PSFILE ) then x_psopen, psfile, /maxs
  !p.multi=[0,1,1]
  clr = getcolor(/load)

  zgrb = 3.03626

  xrng = [5600, 9200.]
  yrng = 10*[0., 0.7]
  ymrg = [4,1]
  if keyword_set(REVIEW) then begin
      xrng = [5600, 7200.]
      yrng = [-2., 7.]
  endif
  ;; Loop
;  xtit = 'Rest Wavelength (Ang)'
  xtit = 'Wavelength (Ang)'
  fx = x_readspec(r400, wav=wv, sig=sig, inflg=2)

  ;; Convert to fnu
  fx = fx * wv * (wv/1d8) / c.c * 1e10 ;; erg/s/cm^2/Hz 1e23
;  yrng = [-0.1, 8.3]
;  yrng = [17.7, 14.5]
;  xrng = [9200, 5600.]

  rwv = wv / (1+zgrb)
  rnu = c.c/rwv ;; Goofy units

  ;; Correct flux for Galactic extinction
  AV_MW = 0.07
  Alambda_MW = AV_MW * x_fitzpalav(wv)
  fx = fx * 10.^(0.4*Alambda_MW)

  ;; Dust
  Al = AV * x_fitzpalav(rwv, R_V=R_V, c3=c3, c4=c4)

  print, 'A_1100 = ', AV * x_fitzpalav(1100., R_V=R_V, c3=c3, c4=c4)
  print, 'A_1400 = ', AV * x_fitzpalav(1400., R_V=R_V, c3=c3, c4=c4)

  ;; Tie to 7140Ang (~1750Ang rest)
  pix = where(wv GT 7124. and wv LT 7152.7)
  medfx = median(fx[pix])
  medwv = median(wv[pix]) / (1+zgrb)
  mednu = c.c/medwv

  Almed = AV * x_fitzpalav(medwv, R_V=R_V, c3=c3, c4=c4)

  trufx0 = medfx * 10.^(0.4*Almed)
  trufx = (rnu/mednu)^(-0.5) * trufx0  ;; Only good for fnu fluxes
;  stop

  ;; Reddened continuum
  red_conti = trufx * 10.^(-0.4*Al)

  ;;

  ;; Continuum
  conti = xmrdfits(cr400, /silen) * 10.^(0.4*Alambda_MW)

  ;; Smooth
;  fx = smooth(fx,SMTH)

  ;; Plot
  if keyword_set(SUB) then begin
      pos = [0.08, 0.10, 0.98, 0.35]
      xmrg = [6,0]
      NOER=1
  endif else begin
      xmrg = [8,2]
  endelse

  plot, [0], [0],  color=clr.black, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, xtitle=xtit, pos=pos, $
        xrange=xrng, yrange=yrng, xtickn=xspaces, NOERASE=NOER, $
        ytickint=2., $
;        ytitle='Magnitude (AB)', $
        ytitle='Relative Flux (erg/s/Ang/cm!u2!N)', $
        /nodata, xstyle=1, ystyle=1


  ;; Convert to magnitude
;  fxmed_AB = -2.5*alog10(medfx*1d-23) - 48.6
;  off = 15. - fxmed_AB
;  fx_AB = -2.5*alog10((fx>1e-3)*1d-23) - 48.6 + off
;  oplot, wv, fx_AB, color=clr.black, psym=10, thick=2
        
  ;; Data
  conv = 1. / wv / (wv/1d8) * c.c / 1e10 ;; erg/s/cm^2/Hz 1e23
  oplot, wv, fx*conv, color=clr.black, psym=10, thick=4
;  oplot, wv, sig, color=clr.green, psym=10, thick=1

  ;; Continuum
;  oplot, wv, conti, color=clr.orange, thick=1, linesty=1
  oplot, wv, red_conti*conv, color=clr.red, thick=3, linesty=1

  ;; Write to disk
  mwrfits, red_conti*conv/(10.^(0.4*Alambda_MW)), 'r400_redconti.fits', /create

  ;; Lines
  lin = [5626.7705d, 5663.8262, 5687.7329,  5731.9604, 5772.6021]
  nlin = n_elements(lin)
  
  asz = 1.7
  lsz = 1.6
  plotsym, 2, asz, thick=4
  offa = 0.20
  offb = 0.40
;  for jj=0L,nlin-1 do begin
;      lbl = strtrim(jj+1,2)
;      yvl = min(fx[where(abs(lin[jj]-wv) LT 5)], imn)
;
;      ;; Plot
;      oplot, [lin[jj]], [yvl-offa], psym=8, color=clr.blue
;      xyouts, lin[jj], yvl-offb, lbl, color=clr.blue, charsiz=lsz, $
;              alignment=0.5
;  endfor

  ;; Atmosphere
  if not keyword_set(REVIEW) then begin
      plotsym, 0, 2., thick=2
      x1 = 7620.
      y1 = 1.4
      if keyword_set(SUB) then dy = 0.34 else dy = 0.09
      dx = 30.
      oplot, [x1], [y1], color=clr.black, psym=8
      oplot, x1+[-1*dx,dx], [y1,y1], color=clr.black, thick=2
      oplot, [x1,x1], y1+[-1*dy,dy], color=clr.black, thick=2
  endif

  ;; Label
;  oplot, [-1e9, 1e9], [0., 0.], color=clr.gray, linesty=2, thick=1

  if keyword_set(REVIEW) then begin
      lsz = 1.6
      ;; Label CO
      xpt = [replicate(1419.044d,2), $
              replicate(1477.565,3), replicate(1509.748,3), replicate(1544.302,2)]
      hrz = -0.7
      ypt = [1.5, hrz, hrz, 1.0, hrz, hrz, 0.5, hrz, hrz, 0.5]

      oplot, xpt*(1+zgrb), ypt, color=clr.darkgreen
      xyouts, 5950., -1.80, 'CO A-X Bandheads', $
              color=clr.darkgreen, charsize=lsz, alignment=0.5

      ;; Label MgII
      z1 = 1.340
      z2 = 1.462
      xpt = [replicate(2796.352d,2)*(1+z1), replicate(2803.531,3)*(1+z1), $
             replicate(2796.352,3)*(1+z2), replicate(2803.531,2)*(1+z2)]
      hrz = -0.6
      ypt = [1.5, hrz, hrz, 1.5, hrz, hrz, 0.5, hrz, hrz, 0.5]

      oplot, xpt, ypt, color=clr.blue
      xyouts, 6750., -1.80, 'MgII:  z=1.340, 1.462', $ 
              color=clr.blue, charsize=lsz, alignment=0.5
  endif

  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]

  return
end
      
      
