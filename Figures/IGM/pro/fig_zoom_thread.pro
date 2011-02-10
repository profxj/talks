;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Zoom in on a piece of spectrum to visualize a Cosmic Web 'thread'
pro fig_zoom_thread

  compile_opt strictarr

  ;; Get structure if necessary
  if not keyword_set( PSFILE ) then psfile = 'fig_zoom_thread.ps'
  if not keyword_set( CSZ ) then csz = 2.3
  if not keyword_set( LSZ ) then lsz = 1.9
  if not keyword_set( BLSZ ) then blsz = 2.0

  if not keyword_set( MAXCHI ) then MAXCHI = 1.3
  if not keyword_set( WVMIN ) then wvmin = 3900. 
  if not keyword_set(OMEGA_M) then omega_m = 0.3
  if not keyword_set(H0) then H0 = 72.
  if not keyword_set(ZQSO) then zqso = 3.8
  if not keyword_set(FWHM) then fwhm = 2.

  xrng1=[800., 1500.]

  ;; Read in Data
  zem = 3.561
  scale = 1e15
  flux = x_readspec('~/Keck/ESI/RedData/Q2223+20/Q2223+20_xF.fits', wav=wave)
  flux = flux * scale
  conti = xmrdfits('~/Keck/ESI/RedData/Q2223+20/Q2223+20_c.fits')
  conti = conti*scale
  rwave = wave / (1+zem)
  gd = where(rwave GT xrng1[0] and rwave LT xrng1[1])

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; PLOT
  x_psopen, psfile, /maxs
  clr = getcolor(/load)
  !p.multi=[0,1,1]
  
  xmrg = [8,2]
  ymrg = [4.0,4]
  ;xrng2=xrng1*(1+zem)
  xrng=[ [4200., 6600], $
         [4200., 6600], $
         [4530., 5735], $
         [4870., 5490], $ 
         [5000., 5300], $
         [5030., 5180], $
         [5125., 5200], $
         [5154., 5173], $
         [5161., 5167]]

  sz = size(xrng, /dimen)

  thk = 5
;  yrng=[-0.02, 3.3]
  yrng=[-0.0, 2.]
  ymin = 0.9
  ymx = [2., 2., 2., 1.3, replicate(ymin,9)]

  for kk=0,sz[1]-1 do begin
     
     yrng[1] = ymx[kk]

     ;; Start plot
     plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
           xmargin=xmrg, ymargin=ymrg, xtitle='Wavelength (Angstroms)', $
           ytitle='Brightness', yrange=yrng, thick=4, $
           xrange=xrng[*,kk], ystyle=1, xstyle=1, psym=1, /nodata ;, /ylog
     
     oplot, wave[gd], flux[gd], color=clr.black, thick=5, psym=10


     ;; Labels
     xlbl = xrng[0,kk] + 0.2*(xrng[1,kk]-xrng[0,kk])
     ylbl = yrng[1]*0.9
     xyouts, xlbl, ylbl, 'Keck/ESI Spectrum', color=clr.black, charsi=lsz, align=0.

     if kk EQ 1 or kk EQ (sz[1]-1) then $
        oplot, [5162., 5162., 5166.5, 5166.5, 5162.], $
               [0., ymin*0.999, ymin*0.999, 0., 0.], color=clr.darkgreen, thick=13

  endfor

  ;; Close Ps
  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]

  print, 'fig_zoom_thread:  All done!'
       
  return
end
      
      
