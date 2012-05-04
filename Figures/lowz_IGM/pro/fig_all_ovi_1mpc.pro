;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro fig_all_ovi_1mpc, summ_fil, ALL_GAL=all_gal, NOPS=nops, NONHI=nonhi, $
                  YMRG=ymrg, CSZ=csz, LSZ=lsz, POS=pos, NOCOS=nocos

  if not keyword_set(binsz) then binsz = 0.01
  if not keyword_set(lsz) then lsz = 1.8 
  if not keyword_set(csz) then csz = 2.3 
  if not keyword_set(dlim) then dlim = 5. ;arcmin

  ;; Input LCO survey
  if not keyword_set(SUMM_FIL) then $
     summ_fil = '/u/xavier/paper/OVI/Galaxies/Analysis/all_galabs_1Mpc_strct.fits'
  struct = xmrdfits(summ_fil,1)
  ngal = n_elements(struct)

  ;; Luminosity cuts
  lum_cuts = [ [0., 0.1], $
               [0.1, 1.], $
               [1., 1e9] ]
  sz = size(lum_cuts, /dimen)

  ;; Plot
  psfile = 'fig_all_ovi_1mpc.ps'
  x_psopen, psfile, /maxs
  !p.multi=[0,1,1]
  clr = getcolor(/load)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; EW of Lya
  xrng = [10., 1300.] ; kpc
  yrng = [3., 1000] ; mA
  xtitle='!9r!X (kpc)'
  bclr = clr.black
  fclr = clr.lightgray
  xmrg = [8,1]
  ymrg = [4,1]

  for ss=0,1 do begin
  
  plot,[0.],[0.],xrange=xrng,yrange=yrng,color=fclr,$
       background=bclr, charsize=csz,  xtickn=spcs, $
       POS=pos, xtitl=xtitle, $
       xmargin=xmrg, ymargin=ymrg, ytitle='W!u1031!N (mA)',$
       /nodata,xstyle=1,ystyle=1, /xlog, /ylog

  all_dra = [0.]
  all_EW = [0.]
  val_dra = [0.]
  val_EW = [0.]
  flg = [0L]

  ;xyouts, 1000., 1300., '(a)', color=clr.black, align=0.5, charsi=lsz

  epsym = 2
  lpsym = 4
  upsym = 1

  ;; Loop on Luminosity
  for qq=0L,sz[1]-1 do begin
     case qq of
        0: begin
           ppsym = 1
           gclr = clr.yellow
           lbl = '< 0.1L*'
        end
        1: begin
           gclr = clr.green
           ppsym = 2
           lbl = '(0.1-1) L*'
        end
        2: begin
           gclr = clr.tomato
           ppsym = 4
           lbl = '> L*'
        end
        else: stop
     endcase

     ;; Limit?
     lim = where(struct.magerr[2] GT 0 AND struct.magerr[5] GT 99. AND $
                 (struct.ddec LT lum_cuts[1,qq] and struct.ddec GT lum_cuts[0,qq]), nlim)
     plotsym, 1, 1.9, thick=5
     if nlim GT 0 then begin
        ;; Plot
        early = where(strmatch(strtrim(struct[lim].gal_type,2), 'Early'), nearly) 
        late = where(strmatch(strtrim(struct[lim].gal_type,2), 'Late'), nlate) 
        unkn = where(strmatch(strtrim(struct[lim].gal_type,2), 'Unkn'), nunkn) 

        if nearly GT 0 then oplot, [struct[lim[early]].dra], [struct[lim[early]].magerr[4]], color=gclr, psym=8, symsiz=sysz
        if nlate GT 0 then oplot, [struct[lim[late]].dra], [struct[lim[late]].magerr[4]], color=gclr, psym=8, symsiz=sysz
        if nunkn GT 0 then oplot, [struct[lim[unkn]].dra], [struct[lim[unkn]].magerr[4]], color=gclr, psym=8, symsiz=sysz

        all_dra = [all_dra, struct[lim].dra]
        all_EW = [all_EW, struct[lim].magerr[4]]
        flg = [flg, replicate(2L,nlim)]
     endif

     ;; Value
     val = where(struct.magerr[2] GT 0 AND struct.magerr[5] LT 99. AND $
                 (struct.ddec LT lum_cuts[1,qq] and struct.ddec GT lum_cuts[0,qq]), nval)
     if nval GT 0 then begin
        val_dra = [val_dra, struct[val].dra]
        val_EW = [val_EW, struct[val].magerr[4]]
        all_dra = [all_dra, struct[val].dra]
        all_EW = [all_EW, struct[val].magerr[4]]
        flg = [flg, replicate(1L,nval)]

        early = where(strmatch(strtrim(struct[val].gal_type,2), 'Early'), nearly) 
        late = where(strmatch(strtrim(struct[val].gal_type,2), 'Late'), nlate) 
        unkn = where(strmatch(strtrim(struct[val].gal_type,2), 'Unkn'), nunkn) 

        if nearly GT 0 then oplot, [struct[val[early]].dra], [struct[val[early]].magerr[4]], color=gclr, psym=epsym, symsiz=sysz
        if nlate GT 0 then oplot, [struct[val[late]].dra], [struct[val[late]].magerr[4]], color=gclr, psym=lpsym, symsiz=sysz
        if nunkn GT 0 then oplot, [struct[val[unkn]].dra], [struct[val[unkn]].magerr[4]], color=gclr, psym=upsym, symsiz=sysz
     endif

     ;; Label
     yval = 80/(1.3^qq)
;     oplot, [15.], [yval], psym=ppsym, color=pclr
;     xyouts, 18., yval/1.05, lbl, color=pclr, charsiz=lsz 
  endfor
  spearman = r_correlate(val_dra[1:*], val_EW[1:*])  ;; Values only
  print, 'All values Spearman: ', spearman

  r200 = where(val_dra GT 200.)
  spear200 = r_correlate(val_dra[r200], val_EW[r200])  ;; Values only
  print, 'r>200 values Spearman: ', spear200

  r200a = where(all_dra GT 200.)
  spear200a = r_correlate(all_dra[r200a], all_EW[r200a])  ;; All
  print, 'All r>200 Spearman: ', spear200a
  writecol, 'rho_vs_EW.dat', all_dra[r200a], all_EW[r200a], flg[r200a]

  ;; Stats
  val_300 = where(val_dra[1:*] LT 300., nval_300)
  all_300 = where(all_dra[1:*] LT 300., nall_300)
  print, 'nval_300: ', nval_300, 'nall_300: ', nall_300

  ;; Median EW
  med_bins = [ [0, 50.], $
               [50., 100], $
               [100., 300], $
               [300., 500], $
               [500., 1000] ]

  szm = size(med_bins, /dimen)
  xmed = [25., 75., 200., 400., 750]
  ymed = fltarr(szm)
  for qq=0L,szm[1]-1 do begin
     gal = where(struct.magerr[2] GT 0 AND struct.dra GT med_bins[0,qq] AND $
                 struct.dra LE med_bins[1,qq], ngal)
     ymed[qq] = median(struct[gal].magerr[4])
  endfor
  plotsym, 0, 1.3, thick=4
  ;oplot, xmed, ymed, color=clr.gray, psym=-8

  ;; Label
  xlbl = 15.
  ylbl = 700.
  xyouts, xlbl, ylbl, 'Dwarf', color=clr.yellow, charsi=lsz
  xyouts, xlbl, ylbl/1.5, 'Sub-L*', color=clr.green, charsi=lsz
  xyouts, xlbl, ylbl/(1.5)^2, 'L*', color=clr.tomato, charsi=lsz

  if ss GE 1 then begin
     oplot, replicate(300., 2), yrng, color=clr.gray, linesty=2, thick=4
  endif

endfor

  x_psclose
  !p.multi=[0,1,1]
  return

end
