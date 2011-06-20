; .com ../Analysis/pro/dla_proxdx
; .com ../Analysis/pro/dla_proxstat
; fig_lox 

pro fig_lls_complox, PSFILE=psfile, GZFIL=gzfil, DLALST=dlalst, $
             RSDSS=rsdss, RTOT=rtot, TIMBIN=timbin, STRCT=strct, $
             BINS=bins, OUTFIL=outfil, VPROX=vprox, $
             VMAX=vmax, PLLS=plls, PROX=prox, QSOFIL=qsofil, $
             LLSFIL=llsfil, DLA_STR=dla_str

  compile_opt strictarr

  ;; Initialize
  cd, getenv('LLSPAP')+'/SDSS/DR7/Analysis/pro/', curr=curr
  RESOLVE_ROUTINE, 'max_lozpower'
  lls_dr7_initparm, init
  cd, curr

  if not keyword_set( LLSFIL ) then llsfil = init.llsfil
  if not keyword_set( QSOFIL ) then qsofil = init.qsofil 
  if not keyword_set( VPROX ) then vprox = init.vprox
  if not keyword_set( XZFIL ) then xzfil = init.xzfil
  if not keyword_set( BINS ) then BINS = init.bins
  if not keyword_set( MAXDZ ) then maxdz = init.maxoff 

  if not keyword_set( PSFILE ) then psfile = 'fig_lls_complox.ps'
  if not keyword_set( CSZ ) then csz = 2.1
  if not keyword_set( LSZ ) then lsz = 2.0

  ;; Truncate the last bin
;  szb = size(bins, /dime)
;  bins = bins[*,0:szb[1]-2]

  ;; Stats
  sdss_llslozx, qsofil, llsfil, strct, BINS=bins, XZFIL=xzfil, VPROX=vprox, $
                MAXDZ=maxdz

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; PLOT
  if keyword_set( PSFILE ) then x_psopen, psfile, /maxs
  clr = getcolor(/load)
  !p.multi=[0,1,1]
  lclr = clr.white
  
  xtit = 'z' 
  yrng = [0.0, 1.3]
  plot, [0.], [0.], color=lclr, background=clr.white, charsize=csz,$
        xmargin=[8.0,1], ymargin=[4,1.], xtitle=xtit, thick=9, $
        /nodata, xrange=[3.4, 4.5], ystyle=1, yrange=yrng, xstyle=1, $
        xtickn=spaces

  ;; Label
;  xyouts, 2.2, yrng[1]*0.85, lbl, color=clr.black, chars=1.8
  
  xlbl = 0.05
  !p.font = -1
  xyouts, xlbl, 0.50, '!12l!X', alignment=0.0, color=lclr,$
          charsize=2.6, /normal, orientation=90
  !p.font = 0
  xyouts, xlbl, 0.51, '(X)', alignment=0.0, color=lclr,$
          charsize=2.2, /normal, orientation=90
  

  ;; LLS (tau>2)
  np = n_elements(strct.medz)
  strct[0].medz = bins[0,0]
  strct[np-1].medz = bins[1,np-1]
  x_curvefill, strct.medz, strct.lox+strct.siglx[0], strct.lox-strct.siglx[1], $
               color=clr.cyan
  lowtau2 = strct.lox-strct.siglx[1]
  xyouts, mean(strct[np-2:*].medz), mean(strct[np-2:*].lox), $
          '!9t!X > 2', color=clr.black, charsi=lsz
      
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; DLAS 
  if not keyword_set( GZFIL ) then gzfil = '~/SDSS/DR5_QSO/dr5_dlagz_s2n4.fits'
  if not keyword_set(DLA_STR) then $
    sdss_dlalox, GZFIL=gzfil, STRCT=dla_str, BINS=bins
  x_curvefill, strct.medz, dla_str.dndx, replicate(0., n_elements(strct)), $
               color=clr.tomato

  xyouts, mean(strct[np-2:*].medz), mean([dla_str.dndx[np-2:*], 0,0])-0.02, $
          'DLAs', color=clr.black, charsi=lsz

  ;;;;;;;;;;;;;;;;
  ;; SLLS (19.0-20.3)
  ;; Increment by lX_SLLS on the DLAs (at all z)
  if not keyword_set(lX_SLLS) then lX_SLLS = init.lX_SLLS 
  slls_dndx = dla_str.dndx + lX_SLLS
  x_curvefill, strct.medz, dla_str.dndx, slls_dndx, color=clr.yellow

  xyouts, mean(strct[np-2:*].medz), mean([dla_str.dndx[np-2:*],slls_dndx[np-2:*]]), $
          'SLLS', color=clr.black, charsi=lsz

  ;; LLS
  x_curvefill, strct.medz, slls_dndx, lowtau2, color=clr.cyan
;  xyouts, mean(strct[np-2:*].medz), mean([lowtau2[np-2:*],slls_dndx[np-2:*]]), $
;          'LLS', color=clr.black, charsi=lsz
  oplot, strct.medz, strct.lox+strct.siglx[0], color=lclr, linest=2, thick=7
  oplot, strct.medz, strct.lox-strct.siglx[1], color=lclr, linest=2, thick=7

;  oploterror, strct.medz, strct.lox, $
;              strct.medz-bins[0,*], strct.siglx[0], $
;              psym=1, color=clr.black, errcolor=clr.black, /lobar
;  oploterror, strct.medz, strct.lox, $
;              bins[1,*]-strct.medz, strct.siglx[1], $
;              psym=1, color=clr.black, errcolor=clr.black, /hibar

  ;; Arrows
  xdla = 4.43
  athick = 11.
  hthick = 3.
  hsz = 500.
  arrow, xdla, 0., xdla, dla_str.dndx[np-1], color=clr.tomato, $
         thick=athick, hthick=hthick, $
         /solid, /data, hsize=hsz

  arrow, xdla, dla_str.dndx[np-1], xdla, dla_str.dndx[np-1]+lX_SLLS, $
         color=clr.yellow, thick=athick, hthick=hthick, $
         /solid, /data, hsize=hsz

  arrow, xdla, dla_str.dndx[np-1]+lX_SLLS, xdla, strct[np-1].lox,  $
         color=clr.cyan, thick=athick, hthick=hthick, $
         /solid, /data, hsize=hsz

  ;; Close Ps
  if keyword_set( PSFILE ) then x_psclose
  print, 'fig_complox:  All done!'

  return
end
      
      
