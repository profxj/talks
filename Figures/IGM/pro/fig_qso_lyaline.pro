;; This shows a QSO SED + single lya line
pro fig_qso_lyaline

  if not keyword_set(CSZ) then csz = 2.1
  if not keyword_set(lSZ) then lsz = 2.1

  ;; Read in the continuum
  qso_fx = xmrdfits('Data/SDSSJ220758.30+125944.3_F_c.fits')
  npix = n_elements(qso_fx)
  CRVAL1  =        3.59106460703 
  CDELT1  =    1.44762400000E-05 
  wave = 10.d^(CRVAL1 + dindgen(npix)*CDELT1)

  zqso = 3.0909
  wrest = wave / (1+zqso)

  ;; Emission Plot
  x_psopen, 'fig_qso_lyaline.ps', /maxs
  !p.multi = [0,1,1]
  clr = getcolor(/load)


  xmrg = [7,7]
  ymrg = [15,1]
  yrng=[0.0, 2.3]
  xrng=[980., 1600]

  lya = x_setline(1215.6701d)
  N1 = 13.5
  lya.b = 30.
  ytit = 'Relative Intensity'
  xtit = 'Rest Wavelength (Ang)'

  for qq=0,2 do begin

     yplt = qso_fx
     ;; Create the spectrum
     if QQ GT 0 then begin
        lya.N = 14.
        sub = where((wrest-1215.) LT 30)
        tmp = x_voigt(wrest[sub], lya, FWHM=3)
        yplt[sub] = tmp * yplt[sub]
     endif

     if qq GE 2 then xrng = [1205, 1227.]
     

     ;; Plot
     plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
           xmargin=xmrg, ymargin=ymrg, ytitle=ytit, $
           xtitle=xtit, yrange=yrng, thick=5, ytickinter=1.0, $
           xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata;, xtickint=1.

     ;; 
     oplot, wrest, yplt, color=clr.black, thick=7

        ;xyouts, 1213.3, 0.30, 'Hydrogen', color=clr.black, charsi=lsz
     xyouts, 1230., 2.0, 'Ly!9a!X', color=clr.blue, charsi=lsz, align=0.
     xyouts, 1550., 1.0, 'CIV', color=clr.blue, charsi=lsz, align=0.5
     xyouts, 1000., 2.0, 'QSO SED', color=clr.blue, charsi=lsz, align=0.

  endfor

  ;; End
  x_psclose
  !p.multi = [0,1,1]


  return

end
