pro fig_lymanwerner

  if not keyword_set(PSFILE) then psfile = 'fig_lymanwerner.ps'

  c = x_constants()
  ;; Read in spectrum
  fx = xmrdfits( getenv('XIDL_DIR')+'/SDSS/LLS/nhi16_19b30.fits',10)
  npnh = n_elements(fx) 
  fx[*] = 1.
  wv_mod = 10^(2.7d + dindgen(npnh)*1e-4)

  ;; H2
  model = { $
        Tex: 0., $
        b: 0., $
        fwhm: 0., $
        NH2: 0. }
  if not keyword_set(DEF_TEX) then def_Tex = 200.  ;; K
  if not keyword_set(DEF_B) then def_b = 3. ;; km/s
  if not keyword_set(DEF_NH2) then def_NH2 = 21.2 ;; log
  if not keyword_set(DEF_FWHM) then def_FWHM = 2. ;; log
  model.NH2 = DEF_NH2
  model.b = DEF_b
  model.fwhm = DEF_FWHM
  model.Tex = DEF_TEX

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
               h2lin.wrest GT min(wv_mod))
  h2lin = h2lin[gdh2]
  
  for ff=0L,nJpp-1 do begin
      idx = where(h2lin.Jpp EQ ff, nidx)
      if nidx EQ 0 then stop
      h2lin[idx].N = alog10( nJ[ff] )
  endfor
;      print, 'NH2 = ', alog10(total(nJ))
  
  h2lin.zabs = 0.  ;; Rest-frame
  h2lin.b = model.b
  fx = x_voigt(wv_mod, h2lin, FWHM=model.fwhm)
  

  ;; Plot
  if keyword_set( PSFILE ) then x_psopen, psfile, /maxs
  !p.multi=[0,1,1]
  clr = getcolor(/load)

  xrng = [900., 1230.]
  yrng = [-0.1, 1.1]
  csz = 2.

  plot, [0],  [0],  color=clr.lightgray, $
    background=clr.black, charsize=csz,$
    xmargin=[7,1.5], ymargin=[5,7], xtitle='Wavelength (Ang)', $
    ytitle='Normalized Flux', /nodata, xthick=7, ythick=7, xstyle=1, ystyle=1, $
    yr=yrng, xr=xrng

  oplot, wv_mod, fx, color=clr.lightgray, thick=7

  ;; Label
  lsz = 2.
;  xyouts, 1200., 0.0, 'Ly!9a!X', color=clr.yellow, charsiz=lsz, alignm=0.5
;  xyouts, 1025., 0.0, 'Ly!9b!X', color=clr.yellow, charsiz=lsz, alignm=0.5
;  xyouts, 972., 0.0, 'Ly!9g!X', color=clr.yellow, charsiz=lsz, alignm=0.5
;  xyouts, 949., 0.0, 'Ly!9d!X', color=clr.yellow, charsiz=lsz, alignm=0.5

  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]

  return
end
