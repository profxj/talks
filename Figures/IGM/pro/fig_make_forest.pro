;; This builds in Lya at z=2, 3, and 4
pro fig_make_forest

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
  zqso = 4.1
  wobs = wrest * (1+zqso) ;; Make believe

  ;; Emission Plot
  x_psopen, 'fig_make_forest.ps', /maxs
  !p.multi = [0,1,1]
  clr = getcolor(/load)


  xmrg = [7,7]
  ymrg = [15,1]
  yrng=[0.0, 2.3]
  xrng=[4960., 6250]

  lya = x_setline(1215.6701d)
  N1 = 13.5
  lya.b = 30.
  ytit = 'Relative Intensity'
  xtit = 'Rest Wavelength (Ang)'

  for qq=0,4 do begin

     if qq LT 4 then begin
        yplt = qso_fx
        
        ;; Create the spectrum
        if QQ GT 0 then begin
           lya.zabs = 4.
           lya.N = 14.
           sub = where(abs(wobs-1215.*(1+lya.zabs)) LT 30)
           tmp = x_voigt(wobs[sub], lya, FWHM=3)
           yplt[sub] = tmp * yplt[sub]
        endif

        ;; Create the spectrum
        if QQ GT 1 then begin
           lya.zabs = 3.5
           lya.N = 14.
           sub = where(abs(wobs-1215.*(1+lya.zabs)) LT 30)
           tmp = x_voigt(wobs[sub], lya, FWHM=3)
           yplt[sub] = tmp * yplt[sub]
        endif
        
        ;; Create the spectrum
        if QQ GT 2 then begin
           lya.zabs = 3.1
           lya.N = 14.
           sub = where(abs(wobs-1215.*(1+lya.zabs)) LT 30)
           tmp = x_voigt(wobs[sub], lya, FWHM=3)
           yplt[sub] = tmp * yplt[sub]
        endif
     endif

        ;; Plot
        plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
              xmargin=xmrg, ymargin=ymrg, ytitle=ytit, $
              xtitle=xtit, yrange=yrng, thick=5, ytickinter=1.0, $
              xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata ;, xtickint=1.
        
        ;; 
        if qq LT 4 then oplot, wobs, yplt, color=clr.black, thick=7
        
        ;; Actual spectrum
        if qq GE 4 then begin
           fx = x_readspec('~/Keck/ESI/RedData/FJ0747+2739/FJ0747+2739_xF.fits', $
                           wav=wv)
           bad = where(abs(wv-5579.4) LT 2.5)
           fx[bad] = 6.2
           wv = wv/1.025
           oplot, wv, fx/6.5, color=clr.black, psym=10, thick=2
        endif
        

  endfor

  ;; End
  x_psclose
  !p.multi = [0,1,1]


  return

end
