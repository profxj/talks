pro qso_history

  if not keyword_set( PSFILE ) then psfile = 'qso_history.ps'
  if not keyword_set(CSZ) then csz = 3.9
  if not keyword_set(lSZ) then lsz = 1.9
  if not keyword_set(THK) then thk = 11


  ;;; BEGIN PLOTS
;  x_psopen, psfile, /maxs
  clr = getcolor(/load)

  zrng = [0.1, 1e5]
  yrng = [2020, 1960]; 1960, 2020]
  xrng = [0.01, 10] ;; Redshift

  SURFACE, DIST(5), /NODATA, /SAVE, YRANGE=yrng, $
           XRANGE=xrng,  ZRANGE=zrng, XSTYLE=8, $
           YSTYLE=1, ZSTYLE=1, CHARSIZE=csz, /zlog, /xlog, $
           ytitle='Year', ztitle='Number', xtitl='Distance'

  axis, xaxis=1, xrang=xrng, /xlog, charsiz=csz, /t3d, xsty=1, $
        xtitl='Distance', color=clr.white
  
  ;; First discovery
  plots, [0.1, 0.1], $
         [1960, 1960], $
         [0.1, 1], color=clr.yellow, /t3d, psym=-4

  ;; LBQS
  plots, [1.0, 1.0], $
         [1990, 1990], $
         [0.1, 1000], color=clr.blue, /t3d, psym=-4

  ;; SDSS
  plots, [3.0, 3.0], $
         [2005, 2005], $
         [0.1, 10000], color=clr.green, /t3d, psym=-4

  ;; SASIR
  plots, [8.0, 8.0], $
         [2020, 2020], $
         [0.1, 10000], color=clr.red, /t3d, psym=-4

;  if keyword_set( PSFILE ) then x_psclose
  !p.multi = [0,1,1]

  return

end
