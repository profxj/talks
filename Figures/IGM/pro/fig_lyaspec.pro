;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro fig_lyaspec, PSFILE=psfile

  compile_opt strictarr

  ;; Get structure if necessary
  if not keyword_set( PSFILE ) then psfile = 'fig_lyaspec.ps'
  if not keyword_set( CSZ ) then csz = 2.3
  if not keyword_set( LSZ ) then lsz = 1.3
  if not keyword_set( BLSZ ) then blsz = 2.0

  if not keyword_set( MAXCHI ) then MAXCHI = 1.3
  if not keyword_set( WVMIN ) then wvmin = 3900. 
  if not keyword_set(OMEGA_M) then omega_m = 0.3
  if not keyword_set(H0) then H0 = 72.
  if not keyword_set(ZQSO) then zqso = 3.8
  if not keyword_set(FWHM) then fwhm = 2.

  ;; Column for 1 Mpc
  c = x_constants()
  NHI = alog10( 0.04 * c.rhoc / c.mp * (1+zqso)^3 * c.Mpc )  + 0.3 ;; Kludge
  hub = cosm_hubble(zqso, /init) ;/10  ;; Another kludge 
  drdz = 3e5/H0/sqrt(omega_m) / (1+zqso)^(2.5) ;; Mpc

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; PLOT
  x_psopen, psfile, /portrait
  clr = getcolor(/load)
  !p.multi=[0,1,7]
  lclr2 = clr.cyan
  lclr = clr.orange
  
  ;; FIRST PLOT ::  Data and best model
  xmrg = [8,2]
  ymrg = [0,0]
  xrng1=[955., 1220.]
  yrng=[0.0, 1.11]

  thk = 3

  for qq=0L,5 do begin

      !p.multi=[0,1,7]
      ;; Dummy
      plot, [0], [0], color=clr.black, background=clr.black, xsty=4, ysty=4

      ;; Start plot
      plot, [0], [0], color=clr.lightgray, background=clr.lightgray, $
            charsize=csz,$
            xmargin=xmrg, ymargin=ymrg, xtitle='', $
            ytitle='Normalized Flux', yrange=yrng, thick=4, $
            xtickn=replicate(' ', 30), $
            xrange=xrng1, ystyle=1, xstyle=5, psym=1, /nodata ;, /ylog

      axis, xaxis=1, charsiz=csz,  xrang=xrng1, color=clr.lightgray, $
            xtitle='!9l!X!drest!N (Angstroms)', xsty=1


      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Data

      xoff = 7.
      ylbl = 0.90
      ylbl2 = 0.65
      ;;;;;;
      ;; PKS0405 for now
      qso_fil = '~/HST/data/PG1718+481/spectra/PG1718+481_E230M_f.fits'
      zqso = 1.084
      fx = x_readspec(qso_fil, wav=wv)
      fx = smooth(fx,1)
      rwv = wv/(1+zqso)
      cut = where(rwv GT 990)
      oplot, rwv[cut], fx[cut], color=clr.lightgray, psym=10, thick=thk
      xyouts, xrng1[0]+xoff, ylbl, 'HST/STIS',$
              color=lclr2, charsiz=lsz
      xyouts, xrng1[0]+xoff, ylbl2, 'z!dq!N=1.08', $
              color=lclr, charsiz=lsz
      
      ;;;;;;
      ;; J2340-00 for now
      if qq GT 1 then begin
          plot, [0], [0], color=clr.lightgray, background=clr.lightgray, $
                charsize=csz,$
                xmargin=xmrg, ymargin=ymrg, xtitle='', $
                ytitle='Normalized Flux', yrange=yrng, thick=4, $
                xtickn=replicate(' ', 30), $
                xrange=xrng1, ystyle=1, xstyle=1, psym=1, /nodata ;, /ylog
          qso_fil = '~/Keck/HIRES/RedData/J2340-00/J2340-00_f.fits'
          zqso = 2.08
          fx = x_readspec(qso_fil, wav=wv)
          keep = where(wv GT 3070.)
          fx = smooth(fx,5)
          rwv = wv/(1+zqso)
          cut = where(rwv GT 990)
          oplot, rwv[cut], fx[cut], color=clr.lightgray, psym=10, thick=thk
          xyouts, xrng1[0]+xoff, ylbl, 'Keck/HIRES', $
                  color=lclr2, charsiz=lsz
          xyouts, xrng1[0]+xoff, ylbl2, 'z!dq!N=2.09', $
                  color=lclr, charsiz=lsz
      endif

      ;;;;;;
      ;; J1558-0031
      if qq GT 2 then begin
          plot, [0], [0], color=clr.lightgray, background=clr.lightgray, $
                charsize=csz,$
                xmargin=xmrg, ymargin=ymrg, xtitle='', $
                ytitle='Normalized Flux', yrange=yrng, thick=4, $
                xtickn=replicate(' ', 30), $
                xrange=xrng1, ystyle=1, xstyle=1, psym=1, /nodata ;, /ylog
          qso_fil = '~/Keck/HIRES/RedData/J1558-0031/J1558-0031_f.fits'
          zqso = 2.83
          fx = x_readspec(qso_fil, wav=wv)
;          keep = where(wv GT 3070.)
;          wv = wv[keep]
;          fx = fx[keep]
          fx = smooth(fx,5)
          rwv = wv/(1+zqso)
          cut = where(rwv GT 990)
          oplot, rwv[cut], fx[cut], color=clr.lightgray, psym=10, thick=thk
          xyouts, xrng1[0]+xoff, ylbl, 'Keck/HIRES', $
                  color=lclr2, charsiz=lsz
          xyouts, xrng1[0]+xoff, ylbl2, 'z!dq!N=2.83', $
                  color=lclr, charsiz=lsz
      endif

      ;;;;;;
      ;; SDSS0127-00
      if qq GT 3 then begin
          plot, [0], [0], color=clr.lightgray, background=clr.lightgray, $
                charsize=csz,$
                xmargin=xmrg, ymargin=ymrg, xtitle='', $
                ytitle='Normalized Flux', yrange=yrng, thick=4, $
                xtickn=replicate(' ', 30), $
                xrange=xrng1, ystyle=1, xstyle=1, psym=1, /nodata ;, /ylog
          qso_fil = '~/Keck/ESI/RedData/SDSS0127-00/SDSS0127-00_f.fits'
          zqso = 4.06
          fx = x_readspec(qso_fil, wav=wv)
;          keep = where(wv GT 3070.)
;          wv = wv[keep]
;          fx = fx[keep]
          fx = smooth(fx,5)
          rwv = wv/(1+zqso)
          cut = where(rwv GT 990)
          oplot, rwv[cut], fx[cut], color=clr.lightgray, psym=10, thick=thk
          xyouts, xrng1[0]+xoff, ylbl, 'Keck/ESI', $
                  color=lclr2, charsiz=lsz
          xyouts, xrng1[0]+xoff, ylbl2, 'z!dq!N=4.06', $
                  color=lclr, charsiz=lsz
      endif

      ;;;;;;
      ;; J1202+3235
      if qq GT 4 then begin
          plot, [0], [0], color=clr.lightgray, background=clr.lightgray, $
                charsize=csz,$
                xmargin=xmrg, ymargin=ymrg, xtitle='', $
                ytitle='Normalized Flux', yrange=yrng, thick=4, $
                xtickn=replicate(' ', 30), $
                xrange=xrng1, ystyle=1, xstyle=1, psym=1, /nodata ;, /ylog
          qso_fil = '~/Keck/ESI/RedData/J1202+3235/J1202+3235a_f.fits'
          zqso = 5.292
          fx = x_readspec(qso_fil, wav=wv)
          fx = smooth(fx,3)
          rwv = wv/(1+zqso)
          cut = where(rwv GT 990)
          oplot, rwv[cut], fx[cut], color=clr.lightgray, psym=10, thick=thk
          xyouts, xrng1[0]+xoff, ylbl, 'Keck/ESI', $
                  color=lclr2, charsiz=lsz
          xyouts, xrng1[0]+xoff, ylbl2, 'z!dq!N=5.29', $
                  color=lclr, charsiz=lsz
      endif
      

      ;;;;;;;;;;;;;;;;;;;;;;
      ;; Labels
;      cclr = clr.orange
;      xyouts, 1180.*(1+zqso), yrng[1]*0.8, $
;              'Ly!9a!X,NV', color=cclr, charsiz=lsz, align=1.0

      xrng2=xrng1*(1+zqso)/1215.6701 - 1.
      axis, xaxis=0, charsiz=csz,  xrang=xrng2, color=clr.lightgray, $
            xtitle='z!dLy!9a!X!N', xsty=1
  endfor

  ;; Close Ps
  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]


  print, 'fig_sdss:  All done!'
       
  return
end
      
      
