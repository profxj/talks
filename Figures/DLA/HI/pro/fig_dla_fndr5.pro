;; fig_fnall, PSFILE='Figures/fig_fnall.ps'
pro fig_dla_fndr5, PSFILE=psfile, GZFIL=gzfil, DLALST=dlalst, ZBIN=zbin, $
               STRCT=strct, GZSTR=gzstr, NOPLT=noplt, PRX=prx, $
               NPLT=nplt, stp=stp, SUBR=subr, NODBL=nodbl, TIFFF=tifff, $
               OUTFIL=outfil, NMAX=nmax, NO_PS=no_ps, N_SNGL=n_sngl, $
               N_GMMA=n_gmma, NOLBL=nolbl, XRNG=xrng, LBLSZ=lblsz, $
               YRNG=yrng, CSZ=csz, N_FIT=n_fit

  ;; Get structure if necessary
  if not keyword_set( PSFILE ) then psfile = 'fig_dla_fndr5.ps'
  if keyword_set(NO_PS) then delvarx, psfile
  if not keyword_set( XZFIL ) then xzfil = getenv('PSDSS')+'/DR5/Analysis/xz_val_L70_M30.fits'
  if not keyword_set( GZFIL ) then $
    gzfil = '~/SDSS/DR5_QSO/dr5_dlagz_s2n4.fits'
  if not keyword_set( NPLT ) then nplt = 18L
  if not keyword_set( STP ) then stp = 0.1
  if not keyword_set( NMIN ) then nmin = 20.3
;  if not keyword_set( NMAX ) then nmax = 1d99
  if not keyword_set( ZBIN ) then zbin = [2.2, 5.5]
  if not keyword_set(LBLSZ) then lblsz = 2.1
  if not keyword_set(CSZ) then csz = 2.1
  if not keyword_set(XRNG) then xrng=[20.2, 22.05]
  if not keyword_set(YRNG) then yrng=[-26.5, -21.]
  
  ;; GZ
  if not keyword_set( GZSTR ) then $
    gzstr = xmrdfits(gzfil, 1, /silent)

  ;; SDSS

  sdss_fndla, GZFIL=GZFIL, STRCT=dr5_fnstr,  STP=stp, XZFIL=xzfil, $
              /SALL, ALLDR=alldr, CL=0.95, NPLT=nplt, NODBL=nodbl

  xval = dr5_fnstr.xval
  yplt = dr5_fnstr.yplt
  bins = dr5_fnstr.bins
  yerr1 = dr5_fnstr.yerr1
  yerr2 = dr5_fnstr.yerr2

  ;; Function
  fplt = 20.2 + (22-20.2)*findgen(1000L)/1000.

  x_psopen, psfile, /portrait
  !p.multi = [0,1,1]
  xmrg = [9,1]
  ymrg = [4.5,1]
  clr = getcolor(/load)
  lclr = clr.white

  plot, [0], [0], color=lclr, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, xtitle='log N!dHI!N', $
        ytitle='log f(N!dHI!N,X)', yrange=yrng, thick=4, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata
  val = where(yplt GT -99., complement=lim, ncomplement=nlim)
  oploterror, [xval[val]], [yplt[val]], [xval[val]-bins[val]],  $
              [yplt[val]-yerr2[val]], psym=3, thick=6, $
              color=lclr, /lobar, errcolor=lclr
  oploterror, [xval[val]], [yplt[val]], [bins[val]+stp-xval[val]],  $
              [yerr1[val]-yplt[val]], psym=3, thick=6, $
              color=lclr, /hibar, errcolor=lclr
  if nlim NE 0 then begin
      plotsym,1, 2., thick=3
      oploterror, [xval[lim]], [yerr1[lim]], [xval[lim]-bins[lim]], $
                  replicate(0.,nlim), $
                  color=lclr, psym=8, /lobar, errcolor=lclr, thick=6
      oploterror, [xval[lim]], [yerr1[lim]], [bins[lim]+stp-xval[lim]], $
                  replicate(0., nlim), $
                  color=lclr, psym=8, /hibar, errcolor=lclr, thick=6
  endif

  ;; Output
  if keyword_set(OUTFIL) then begin
      outx1 = xval-(xval-bins)
      outx2 = xval+(bins+stp-xval)
      outy = yplt
      outy1 = yplt-yerr2
      outy2 = yerr1-yplt
      if nlim NE 0 then begin
          outy[lim] = yerr1[lim]
          outy1[lim] = -99
          outy2[lim] = -99
      endif
      ;; Write
      writecol, outfil, outx1, outx2, outy, outy1, outy2
  endif

  ;; Single power law

  N_SNGL = 1
  if not keyword_set(N_SNGL) then begin
      thy = dr5_fnstr.k1 * (10.^fplt)^(dr5_fnstr.a1)
      oplot, fplt, alog10(thy), linestyle=4, color=clr.green
      oplot, [20.3,20.45], [-24.4,-24.4], linestyle=4, color=clr.green
      xyouts, 20.5, -24.5, 'Single Power-Law', color=clr.green, $
              charsize=lblsz
  endif


  N_GMMA = 1
  ;; Gamma function
  if not keyword_set(N_GMMA) then begin
      thy = dr5_fnstr.k2 * (10.^fplt/dr5_fnstr.Ng)^(dr5_fnstr.a2)* $
            exp(-1.*(10.^fplt/dr5_fnstr.Ng))
      oplot, fplt, alog10(thy), linestyle=2, color=clr.red
      oplot, [20.3,20.45], [-24.8,-24.8], linestyle=2, color=clr.red
      xyouts, 20.5, -24.9, '!9G!X-Function', color=clr.red, $
              charsize=lblsz
  endif
  
  
  ;; Double Power-Law
  if not keyword_set(N_FIT) then begin
      thy = fltarr(n_elements(fplt))
      ns10 = 10.^dr5_fnstr.Nd
      lw = where(fplt LT dr5_fnstr.Nd, complement=hi, ncomplement=numhi)
      thy[lw] = dr5_fnstr.k3 * (10^fplt[lw]/ns10)^dr5_fnstr.a3
      if numhi NE 0 then thy[hi] = dr5_fnstr.k3 * $
        (10^fplt[hi]/ns10)^dr5_fnstr.a4
      oplot, fplt, alog10(thy), linestyle=3, color=clr.cyan
      
      oplot, [20.70,20.88], [-21.45,-21.45], linestyle=3, color=clr.cyan
      xyouts, 20.92, -21.5, 'Double Power-Law', color=clr.cyan, $
              charsize=lblsz
  endif
  
  ;; Close Ps
  if keyword_set( PSFILE ) then x_psclose
  if not keyword_set( SUBR ) then !p.multi = [0,1,1]
  
  print, 'fig_fnall: All done!'

  return
end
      
      
