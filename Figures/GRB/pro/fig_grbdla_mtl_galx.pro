;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; fig_cheme, psfile='Figures/fig_cheme.ps', DLA=dla
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro boot_z, dla, Z, sigZ, seed

  if not keyword_set( NTRIAL ) then ntrial = 1000L
  if not keyword_set( CL ) then cl = 0.83

  ;; 
  ndla = n_elements(dla)
  sv_mean = fltarr(ntrial)
  sigZ = fltarr(2)

  for q=0,ntrial-1 do begin
      ;; Random DLA
      rdla = long(randomu(seed, ndla)*ndla)
      rsign = randomu(seed, ndla)
      rsigNH = randomn(seed, ndla)
      rsigMH = randomn(seed, ndla)
      a = where(rsigN LT 0.5,na)
      sign = replicate(1., ndla)
      if na NE 0 then sign[a] = -1.

      ;; HI
      HI = dla[rdla].NHI + dla[rdla].sigNHI[1]*sign*rsigNH
      if na ne 0 then $
        HI[a] = dla[rdla[a]].NHI + dla[rdla[a]].sigNHI[0]*sign[a]*rsigNH[a]
      top = total(10^(dla[rdla].mtl + dla[rdla].sigmtl*rsigMH + HI))
      bot = total(10^HI)

      sv_mean[q] = top/bot
  endfor
      
  ;; Central value
  top = total( 10^(dla.mtl + dla.NHI))
  bot = total( 10^(dla.NHI))
  Z = alog10(top/bot)
  print, 'Tot NHI = ', alog10(bot)
            
  ;; STATS
  srt = sort(sv_mean)
  sv_mean = sv_mean[srt]
  sigZ[0] = Z - alog10(sv_mean[long(ntrial*(1.-cl))])
  sigZ[1] = alog10(sv_mean[long(ntrial*cl)]) - Z

  return
end
;;;;;;;
pro boot_sigmean, dla, mn, sigmn, seed

  if not keyword_set( NTRIAL ) then ntrial = 1000L
  if not keyword_set( CL ) then cl = 0.975
  ndla = n_elements(dla)
  sv_mean = fltarr(ntrial)
  sigmn = fltarr(2)

  ;; MONTE
  for q=0L,ntrial-1 do begin
      sum = 0.
      
      rdla = long(randomu(seed, ndla)*ndla)
      rsig = randomn(seed, ndla)
      sum = total( dla[rdla].mtl + dla[rdla].sigmtl*rsig )
      sv_mean[q] = sum / float(ndla)
  endfor
      
  ;; DLA mean
  mn = total( dla.mtl ) / float(ndla)
  ;; Grab c.l.
  srt = sort(sv_mean)
  sv_mean = sv_mean[srt]
  sigmn[0] = mn - sv_mean[long(ntrial*(1.-cl))]
  sigmn[1] = sv_mean[long(ntrial*cl)] - mn

  return
end
;;;
; talk_mtl, psfile='talk_grbmtl.ps', /NODLA
pro fig_grbdla_mtl_galx, PSFILE=psfile, NBIN=nbin, SURVIVE=SURVIVE, $
              ZMAX=zmax, ZMIN=zmin, NODLA=nodla, LBG=lbg, $
              NARROW=narrow, NONO=nono, XRNG=xrng, BW=bw, TWOPG=twopg, $
              CSZ=csz, LSZ=lsz

  if not keyword_set(PSFILE) then psfile = 'fig_grbdla_mtl_galx.ps'
  if not keyword_set( ZMAX ) then zmax = 999.
  if not keyword_set( ZMIN ) then zmin = 1.65
  if not keyword_set( CSZ ) then csz = 3.1
  if not keyword_set( LSZ ) then lsz = 1.8
  if not keyword_set( NBIN ) then nbin = 5
  ;; Parse DLA files
  if not keyword_set(DLA) then parse_dlalst, dla, $
    '/u/xavier/DLA/Lists/metal704_dla.lst', /noel, ROOT='/u/xavier/DLA/'
  ;; Parse GRB files
  if not keyword_set(GRB) then parse_dlalst, grb, $
    '/u/xavier/DLA/Lists/grb_dla.lst', /noel, ROOT='/u/xavier/DLA/'

  ;; Plot
  gdd = where(dla.NHI GE 20.299 AND dla.flgmtl NE 0 AND $
              dla.zabs GE ZMIN and dla.zabs LT ZMAX, ndla)
  srt = sort(dla[gdd].zabs)
  dla = dla[gdd[srt]]

  ;; GRB
  gdd = where(grb.NHI GE 20.3 AND grb.flgmtl NE 0 and grb.zabs LT ZMAX, ngrb)
  grb = grb[gdd]

  tot_age = cosm_time(100.,/W06MAP, /INIT)/1e9
  age_d = tot_age-cosm_time(dla.zabs)/1e9
  age_g = tot_age-cosm_time(grb.zabs)/1e9

  close, /all
  openw, 44, 'dla.dat'

  ;; Plot
  if keyword_set( PSFILE ) then x_psopen, psfile, /maxs
  clr = getcolor(/load)
  if keyword_set(TWOPG) then !p.multi=[0,1,2] else !p.multi=[0,1,1]


  ymrg = [15,4]
  csz = 1.7

  if not keyword_set(XRNG) then xrng = [5, 0.]
  yrng = [-3., 0.5]
  ;; Overall
  plot, dla.zabs, dla.mtl, color=clr.lightgray, $
    background=clr.black, charsize=csz,$
    xmargin=[9,1], ymargin=ymrg, xtitle='Age of Universe (Gyr)', $
    ytitle='[M/H]', /nodata, xthick=7, ythick=7, xstyle=9, ystyle=1, $
    yr=yrng, xr=xrng

  xspaces = replicate(' ',30) 
  axis, xaxis=1, charsize=csz, xstyle=1, xtickn=xspaces, xminor=1, xthick=5, $
        color=clr.lightgray

  ;; DLA
  if keyword_set(BW) then qclr = clr.gray else qclr = clr.green

  for ii=0L,ndla-1 do begin
      ssize = 0.3 + (10^(dla[ii].NHI-20.3))/10.
;      plotsym, 8, 0.5;, ssize
;      oplot, [dla[ii].zabs], [dla[ii].mtl], color=clr.darkgreen, psym=1 ;green, psym=8
;      if dla[ii].nhi GE 21 then symp = 4 else symp = 1
      symp = 1 
      if not keyword_set(NODLA) then $
        oplot, [age_d[ii]], [dla[ii].mtl], color=qclr,  $
               psym=symp, symsiz=0.7        ;green, psym=8
   endfor
;  writecol, blah, dla.zabs, dla.mtl, FILNUM=44 

  ;; Cosmic average
  nbin = 3
  for jj=0L,nbin-1 do begin
      ;; Grab em
      gdi = where(age_d GT (jj+1) and age_d LT (jj+2), ngdi)
      if ngdi EQ 0 then continue
      ;; Stats
      boot_z, dla[gdi], Z, sigZ, seed
      printf, 44, jj+1, jj+2, Z, sigZ[0], sigZ[1]
      ;; Plot
      plotsym, 8, 2, /fill
      if not keyword_set(NODLA) then begin
          oploterror, [jj+1.5], [Z], [0.5], [sigZ[0]], $
                      color=qclr,  psym=1, /lobar, thick=8, $
                      errcolor=qclr 
          oploterror, [jj+1.5], [Z], [0.5], [sigZ[1]], $
                      color=qclr, psym=1, /hibar, thick=8, $
                      errcolor=qclr 
      endif
  endfor



  ;; GRB
  if keyword_set(BW) then gclr = clr.black else gclr = clr.cyan
  flg_lim = intarr(ngrb)
  for ii=0L,ngrb-1 do begin
      ssize = 0.3 + (10^(grb[ii].NHI-20.3))/10.
      if grb[ii].sigNHI[0] GE 0.6 then begin
          if keyword_set(NONO) then continue
          symsz=1. 
      endif else symsz=3.5

      if grb[ii].flgmtl LT (-10) or grb[ii].sigNHI[0] GE 0.6 then begin
          plotsym, 2, symsz,thick=6       ;, ssize
;          oplot, [grb[ii].zabs], [grb[ii].mtl], color=gclr, psym=8
          oplot, [age_g[ii]], [grb[ii].mtl], color=gclr, psym=8
          flg_lim[ii] = 1
      endif else begin
;          oploterror, [grb[ii].zabs], [grb[ii].mtl], [grb[ii].sigmtl], $
          oploterror, [age_g[ii]], [grb[ii].mtl], [grb[ii].sigmtl], $
            color=gclr, errcolor=gclr, psym=2
      endelse
  endfor

  LBG = 1
  ;; LBG
  if keyword_set(LBG) then begin
      age1 = tot_age-cosm_time(2)/1e9
      oploterror, [age1], [-0.3], [0.2], color=clr.yellow, psym=4, errcol=clr.yellow
      age2 = tot_age-cosm_time(3)/1e9
      oploterror, [age2], [-0.5], [0.2], color=clr.yellow, psym=4, errcol=clr.yellow
  endif

  ;; Solar
  oplot, [-1e9, 1e9], [0, 0.], color=clr.lightgray, linest=2, thick=5

  ;; Label
;  oplot, [4.7], [-2.65], color=clr.green, psym=1
  lsz = 1.8
  xyouts, xrng[0]-0.2, -2.7, 'QSO-DLA', color=qclr,  charsiz=lsz, align=0.
  oplot, [4.7], [-2.40], color=clr.cyan, psym=1
  xyouts, xrng[0]-0.2, -2.35, 'GRB-DLA', color=gclr, charsiz=lsz, align=0.
  if keyword_set(LBG) then $
    xyouts, xrng[0]-0.2, -2.0, '(Bright) LBG', color=clr.yellow, charsiz=lsz, align=0.

  ;; Label the top axis
  xgrd = 1. + 5.5* findgen(1000L)/1e3
  tgrd = tot_age - cosm_time(xgrd)/1e9
  lblt = lindgen(round(xrng[0]+1))
  lblt = lblt[1:*]
  nlbl = n_elements(lblt)
  for jj=0L,nlbl-1 do begin
      mn = min(abs(tgrd-lblt[jj]),imn)
      xyouts, lblt[jj], yrng[1]+0.1, string(xgrd[imn],format='(f3.1)'), $
              color=clr.lightgray, charsize=csz, alignment=0.5
  endfor
  xyouts, xrng[0]/2,  yrng[1]+0.4, 'z!dabs!N', color=clr.lightgray, charsiz=csz,$
          alignment=0.5

  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]

  close, /all

  return
end
