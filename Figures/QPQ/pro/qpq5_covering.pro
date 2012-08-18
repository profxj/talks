;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  Plots OThick absorbers
pro qpq5_covering, COVFACT = COVFACT, psfile = psfile, outfile = outfile, NOPS=nops, $
                  CSZ=csz, XMRG=xmrg

  ;; Get structure if necessary
  if not keyword_set( PSFILE ) then psfile = 'qpq5_covering.ps'
  if not keyword_set( CSZ ) then csz = 1.9
  if not keyword_set( LSZ ) then lsz = 1.9

  if not keyword_set(SIG_FX) then sig_fx = 0.05 ;; Error in tau values
  if not keyword_set( MAXCHI ) then MAXCHI = 1.3
  ;  if not keyword_set( WVMIN ) then wvmin = 700. 
  if not keyword_set(OMEGA_M) then omega_m = 0.3

  if not keyword_set( NRNG ) then nrng = 50L

  ;; Read in pair stuff
  qpq_fil = '~/Dropbox/QSOPairs/qpq5_pairs.fits'
  qpq_strct = xmrdfits(qpq_fil, 1)
  npair = n_elements(qpq_strct)

  COVFACT = 1

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; PLOT
  if not keyword_set(NOPS) then begin 
     x_psopen, psfile, /maxs
     !p.multi=[0,1,1]
     xmrg = [7,8]
     ymrg = [4.0,0.5]
     xtit='R!dperp!N (kpc)'
  endif else begin
     if not keyword_set(XMRG) then xmrg = [7,6]
     ymrg = [0., 0]
     xspaces = replicate(' ',30) 
     xtit = ''
  endelse
  clr = getcolor(/load)
  ;lclr = clr.lightgray
  lclr = clr.white

  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Covering fraction

  xrng = [-10., 310]
  yrng = [0., 1.05]

  ytit = 'Covering Fraction'
  plot, [0], [0], color=lclr, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, xtitle=xtit, $
        ytitle=ytit, yrange=yrng, thick=4, xtickn=xspaces, $
        xrange=xrng, ystyle=9, xstyle=1, psym=1, /nodata

  ;; Label
  ;if keyword_set(NOPS) then xyouts, 01., 0.9, '(a)', color=lclr, charsi=lsz

  eclr = clr.yellow
  binsz = 50.

  ;; Loop on Radii
  for kk=0L,round(xrng[1]/binsz)-1 do begin
     x1 = binsz*kk
     x2 = x1 + binsz
     
     gdi = where(qpq_strct.R_phys GE x1 and qpq_strct.R_phys LT x2, npts)
     gdot = where(qpq_strct[gdi].flg_OThick EQ 1, nothk)

     ;; Stats
     cov_frac = float(nothk)/float(npts)
     npoiss = x_poisscl(nothk,sigma=1)
     xbin = mean([x1,x2])

     ;; Plot
     plotsym, 8, 1.7, /fill
     oploterror, [xbin], [cov_frac], [abs(xbin-x1)], [cov_frac - npoiss[1]/float(npts)], $
                 psym=8, errcolor=eclr, errthick=5, /lobar, color=eclr
     oploterror, [xbin], [cov_frac], [abs(xbin-x2)], [(abs(cov_frac - npoiss[0]/float(npts))) < (1.-cov_frac)], $
                 psym = 8, errcolor = eclr, errthick = 5, /hibar, color = eclr

     plotsym, 2, 3.7, thick=5
     oplot, [xbin], [cov_frac], psym=8, color=eclr
     
     ;; ;;;;;;;;;;;;;;;
     ;; Random IGM 

     ;; Ribaudo+11
     lstar = 1.62
     zstar = 3.23
     gamma = 1.83
     
     vsig = 1500.  ;; Error for QSO redshift (km/s)
     z1 = x_relvel(qpq_strct[gdi].z_fg, vsig)
     z2 = x_relvel(qpq_strct[gdi].z_fg, -1*vsig)

     const = lstar / (1+gamma) / (1+zstar)^gamma
     NLLS = const * ( total( (1+z2)^(gamma+1) ) - total( (1+z1)^(gamma+1) ) )

     ;; Assume 20% error
     oploterror, [xbin], [NLLS/float(npts)], 0.2*[NLLS/float(npts)], errcol=clr.cyan
  endfor
  
  xyouts, 200., 0.95, 'Optically Thick Gas', color=clr.yellow, charsi=lsz

  ;; Add Rudie et al. 2012 :: r<R_vir = 30 +/- 14%
  lbgc = clr.red
  oploterror, [75.], [0.3], [25.], [0.14], color=lbgc, errcol=lbgc
  plotsym, 3, 1.7, /fill
  oplot, [75.], [0.3], color=lbgc, psym=8

  ;; Stats to R=200kpc
  gdi = where(qpq_strct.R_phys LT 200., npts)
  gdot = where(qpq_strct[gdi].flg_OThick EQ 1, nothk)
  print, npts, ' systems with R<200 kpc and ', nothk, ' are othick.'
  
  yrng = [1.5, 3.5]
  ytit='z!dfg!N'
  axis, yaxis = 1, color = lclr, charsi = csz, yrang = yrng, $
        ysty = 1, ytit =ytit,  /save

  ;; May need to expand beyond Lya, e.g. to include MgII
  gd = where(qpq_strct.flg_EWLya GT 0, ngd)
  qpq_strct = qpq_strct[gd]
  
  ;; Loop on velocity offset
  vmnx = [ [-10000.,-1000], $
           [-1000, -400], $
           [-400, 400], $
           [400, 1000], $
           [1000, 10000] ]
  v_clrs = [clr.blue, clr.green, clr.black, clr.orange, clr.red]
  ;; JFH Hack for talk
  v_clrs = [clr.black, clr.black, clr.black, clr.black, clr.black]
  vrel = x_relvel(qpq_strct.z_Lya, qpq_strct.z_fg, /rev)

  xlbl = 210.
  ylbl = yrng[1]-1.
  yoff = 1.3
  for ii=0L,n_elements(v_clrs)-1 do begin
     ;; Any?
     idx = where(vrel GT vmnx[0, ii] and vrel LE vmnx[1,ii], na)
     for jj=0L,na-1 do begin
        a = idx[jj]
        if qpq_strct[a].flg_OThick EQ 1 then FILL = 1 else FILL = 0
        if qpq_strct[a].flg_OThick LT 0 then vclr = clr.darkgray else vclr =lclr
        plotsym, 0, 0.7, FILL=fill, thick=3
        oplot, [qpq_strct[a].R_phys], [qpq_strct[a].z_Lya], color=vclr, psym=8
        ;; Outliers
        ;if (abs(vrel[a]) GE 1000) and FILL then print, NHI_qso[a], NHI_zfg[a], vrel[a]
     endfor
     ;; Label
;     xyouts, xlbl, ylbl-yoff*ii, strtrim(round(vmnx[0,ii]),2)+'< !9d!Xv < '+$
;             strtrim(round(vmnx[1,ii]),2), color=v_clrs[ii], charsiz=lsz
  endfor


  ;; Close Ps
  if not keyword_set(NOPS) then begin 
     if keyword_set( PSFILE ) then x_psclose
     !p.multi=[0,1,1]
  endif

  print, 'fig_covering:  All done!'
  return
end
      
      
