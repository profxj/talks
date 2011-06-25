pro fig_gmos_spectra

  if not keyword_set(TRIM) then trim = 10L
  if not keyword_set(CSZ) then csz = 2.5
  if not keyword_set(CSZ2) then csz2 = 1.3
  if not keyword_set(lSZ) then lsz = 2.0
  if not keyword_set(lSZ2) then lsz2 = 0.7
  if not keyword_set(XNCOLORS) then xncolors=200L
  if not keyword_set(TRIM_PIX) then trim_pix = 45L
  if not keyword_set( NTOT ) then ntot = 14L
  if not keyword_set( XRNG ) then xrng = [2000, 4900.]
  if not keyword_set( DELTA_N ) then delta_N = 0.05 ; Column (dex)

  ;; Data
  dat_dir = '~/Dropbox/gmos/examplespectra/'

  ;; Read in quasars
  files = findfile(dat_dir+'*B600*', count=nfil)

  npx = 1 
  npy = 6

  x_psopen, 'fig_gmos_spectra.ps', /portrait

  ;;starting at p.multi[0], plot p.multi[1] rows and p.multi[2] columns
  ;;with p.multi[3] z dimensions and going from top to bottom (p.multi[4]
  !p.multi=[npx*(npy+1),npx,npy+1,0,1]
  clr = getcolor(/load)
  lclr = clr.white

  for ii=0L,npy-1 do begin
     
     if ii EQ npy then !p.multi=[npy+1,npx,npy+1,0,1]

     ;; QSO
     fx = x_readspec(files[ii], wav=wave, infl=2, sig=sig)
        
     yrng = [-0.05*max(fx), max(fx)*1.2]
     xrng = [5000., 7900]
        ;;;;;;;;;;;;;;;;;;;;
     ;; Plot the data
     if (ii MOD npy) NE (npy-1) and (ii NE npy-1) then begin
        xspaces = replicate(' ',30) 
        xtitle = ''
     endif else begin 
        if keyword_set(XSPACES) then delvarx, xspaces
        xtitle='Wavelength (Ang)'
     endelse
     if (ii GT npy-1) then yspaces = replicate(' ',30) 
     plot, [0], [0], xrange=xrng, $
           yrange=yrng, xtickn=xspaces, xmargin=[6,0], $
           ymargin=[0,0], /NODATA, $ ;ytickn=yspaces, $
           charsize=csz, psym=10, background=clr.white, color=lclr, $
           xstyle=1, ystyle=1, ytickinterval=ytint,$
           xtitle=xtitle, thick=9 ;, xtickint=400.

     ;; Data
     oplot, wave, sig, color=clr.tomato, psym=10, thick=1
     oplot, wave, fx, color=lclr, psym=10, thick=5

     ;; Label
     ipos1 = strpos(files[ii], '/', /reverse_search)
     ipos2 = strpos(files[ii], '_B600')
     xyouts, 5100., yrng[1]*0.65, strmid(files[ii],ipos1+1,ipos2-ipos1-1), $
             color=clr.yellow, charsi=lsz, align=0.
  endfor

  xyouts, 0.03, 0.4, 'Relative Flux (10!u-17!N erg s!u-1!N cm!u-2!N A!u-1!N)', $
          color=lclr, align=0., orient=90., /normal, charsiz=1.6
  x_psclose
  !p.multi = [0,1,1]

  return

end
