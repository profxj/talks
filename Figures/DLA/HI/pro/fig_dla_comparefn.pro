;; fig_fnall, PSFILE='Figures/fig_fnall.ps'
pro fig_dla_comparefn, PSFILE=psfile, GZFIL=gzfil, DLALST=dlalst, ZBIN=zbin, $
                STRCT=strct, GZSTR=gzstr, NOPLT=noplt, PRX=prx, $
                NPLT=nplt, stp=stp, SUBR=subr, NODBL=nodbl, TIFFF=tifff, $
                OUTFIL=outfil, NMAX=nmax, NO_PS=no_ps, CSZ=csz, LBLSZ=lblsz,$
                NMIN=nmin, NORAO=norao

  ;; Get structure if necessary
  if not keyword_set( PSFILE ) then psfile = 'fig_dla_compare_fn.ps'
  if keyword_set(NMIN) then psfile = 'fig_cut_evolfn.ps'
  if not keyword_set(NMIN) then NMIN = 0.
  if not keyword_set( XZFIL ) then xzfil = getenv('PSDSS')+'/DR5/Analysis/xz_val_L70_M30.fits'
  if not keyword_set( GZFIL ) then $
    gzfil = '~/SDSS/DR5_QSO/dr5_dlagz_s2n4.fits'
  if not keyword_set( NPLT ) then nplt = 17L
  if not keyword_set( STP ) then stp = 0.1
  if not keyword_set(LBLSZ) then lblsz = 2.5
  if not keyword_set(CSZ) then csz = 1.9
  
  seed = -1211L
  ;; GZ
  if not keyword_set( GZSTR ) then $
    gzstr = xmrdfits(gzfil, 1, /silent)

  if not keyword_set( BINS ) then $
    bins = [ [2.2, 2.4], [2.4,2.7], [2.7, 3.0], [3.0, 3.5], [3.5,4.0], [4.0,5.5]]
  if n_elements(bins) EQ 2 then sz = [2,1] else sz = size(bins,/dimensions)
  nbin = sz[1]

  ;; SDSS
  sdss_dlastrct, dla, /all, /stat

  ;; CUT DLA HERE
  dla = dla[where(dla.NHI GT NMIN)]
  ndla = n_elements(dla)
  print, 'fig_evolfn: Ndla = ', ndla

  npt = 100L
  Nval = 20.299 + findgen(npt)*1.7/float(npt-1)
  cumu = fltarr(npt)

  x_psopen, psfile, /portrait
  xmrg = [9,1]
  ymrg = [4.5,1]

  clr = getcolor(/load)
  lclr = clr.white

  plot, [0], [0], color=lclr, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, xtitle='log N!dHI!N', $
        ytitle='Cumulative Fraction', yrange=[0, 1.05], thick=8, $
        xrange=[20.3, 22.05], ystyle=1, xstyle=1, psym=1, /nodata

  ;; Bins
  y1 = 0.90
  stp = 0.065
  off = 0.01
  for qq=0L,nbin-1 do begin
;  for qq=0L,0 do begin
      ;; Cut
      idx = where(dla.zabs GE bins[0,qq] and dla.zabs LT bins[1,qq] AND $
                  dla.NHI GT NMIN, nsub)
      for jj=0L, npt-1 do begin
          a = where(dla[idx].NHI LE Nval[jj], na)
          cumu[jj] = na/float(nsub)
      endfor
      ;; Plot
      psty = (qq MOD 4) + 1
      if qq GT 3 then pclr = clr.tomato else pclr=clr.yellow
      oplot, nval, 1.-cumu, color=pclr, psym=10, linesty=psty, thick=5
      ;; KS Test
      if nsub GT 4 then begin
          kstwo, dla.NHI+0.05*randomu(seed, ndla)-0.025, $
                 dla[idx].NHI+0.05*randomu(seed, nsub)-0.025, $
                 d, prob
          print, bins[0,qq], bins[1,qq], prob
      endif else print, 'fig_evolfn: Not enough DLAs for KS', nsub
      ;; Label
      oplot, [21.25, 21.35], y1-[stp,stp]*qq, color=pclr, linesty=psty, thick=5
      xyouts, 21.4, y1-stp*qq-off, 'z=['+string(bins[0,qq],format='(f3.1)')+$
              ','+string(bins[1,qq],format='(f3.1)')+')', color=pclr, $
              charsize=lblsz
  endfor

  ;; Total SDSS-DR5
  for jj=0L, npt-1 do begin
      a = where(dla.NHI LE Nval[jj], na)
      cumu[jj] = na/float(ndla)
  endfor
  oplot, Nval, 1.-cumu, color=lclr, psym=10, thick=7
  oplot, [21.25, 21.35], y1-[stp,stp]*nbin, color=lclr
  xyouts, 21.4, y1-stp*nbin-off, 'z=[2.2,5.5)', color=lclr,$
            charsize=lblsz

  ;; Rao?
  if keyword_set(RAO) then begin
      cd, '../../Talks/pro/', CURRENT=current
      RESOLVE_ROUTINE, 'parse_rao06'
      cd, current
      parse_rao06, rao_dla
      nrao = n_elements(rao_dla)
      for jj=0L, npt-1 do begin
          a = where(rao_dla.NHI LE Nval[jj], na)
          cumu[jj] = na/float(nrao)
      endfor
      nbin = nbin+1
      oplot, Nval, 1.-cumu, color=clr.gray, psym=10, thick=5, linesty=1
      oplot, [21.25, 21.35], y1-[stp,stp]*nbin, color=clr.gray, linesty=1
      xyouts, 21.4, y1-stp*nbin-off, 'z=[0.2,1.6)', color=clr.gray, charsize=lblsz
      
      kstwo, dla.NHI+0.05*randomu(seed, ndla)-0.025, $
             rao_dla.NHI+0.05*randomu(seed, nsub)-0.025, d, prob
      print, 'Rao: ', prob
  endif
      
  ;; Zwaan
  beta = 1.52
  log_NHIstar = 21.3
  A = 1 - beta

  cumu = ( igamma(A,10.d^(Nval-log_NHIstar), /doubl) $
           - igamma(A,10.d^(20.3 - log_NHIstar), /doubl) ) / $
         ( 1 - igamma(A,10.d^(20.3-log_NHIstar)) )
  nbin = nbin+1
  oplot, Nval, 1.-cumu, color=clr.cyan, psym=10, thick=7
  oplot, [21.25, 21.35], y1-[stp,stp]*nbin, color=clr.cyan
  xyouts, 21.4, y1-stp*nbin-off, 'z=0', color=clr.cyan, charsize=lblsz

  
  ;; Close Ps
  if not keyword_set( NO_PS ) then begin
      x_psclose
      !p.multi = [0,1,1]
  endif

  print, 'fig_fnall: All done!'

  return
end
      
      
