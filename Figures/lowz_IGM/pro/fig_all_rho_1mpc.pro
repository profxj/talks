;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro fig_all_rho_1mpc, summ_fil, ALL_GAL=all_gal, NOPS=nops, NONHI=nonhi, $
                  YMRG=ymrg, CSZ=csz, LSZ=lsz, POS=pos, NOCOS=nocos

  if not keyword_set(binsz) then binsz = 0.01
  if not keyword_set(lsz) then lsz = 1.5 
  if not keyword_set(csz) then csz = 2.3 
  if not keyword_set(dlim) then dlim = 5. ;arcmin

  ;; Input LCO survey
  if not keyword_set(SUMM_FIL) then $
     summ_fil = '/u/xavier/paper/OVI/Galaxies/Analysis/all_galabs_1Mpc_strct.fits'
  struct = xmrdfits(summ_fil,1)
  ngal = n_elements(struct)

  ;; Input COS-Halos Survey
  ldir = getenv('DROPBOX_DIR')+'/COS-Halos/lowions/'
  restore, ldir+'/cosmetals_megastructure.sav'

  ;; The Lyalpha equivalent width is then:
  tname = strtrim(megastruct.ion[1,1].trans.name, 2)
  tlam =  strtrim(string(megastruct.ion.trans.lambda, format = '(i4)'), 2)
  lya = where (tname EQ 'HI' AND tlam EQ '1215')
  lyaind = (lya / 27) MOD 16 
  COS_EW = megastruct.ion[1,1].trans[lyaind].wrest  ;; mA
  COS_sigEW = megastruct.ion[1,1].trans[lyaind].sigwrest  ;; mA

  gdcos = where(cos_sigEW GT 0., ncos)
  cos_str = replicate(struct[0], ncos)

  cos_str.ddec = 1.5 ;; L*
  cos_str.mag[2] = 1 ;; Flag (I think)
  cos_str.mag[4] = COS_EW[gdcos] ;; mA
  cos_str.mag[5] = COS_sigEW[gdcos] ;; mA
  cos_str.dra = megastruct[gdcos].rhofinal ;; kpc
  lim = where(cos_ew[gdcos] LT 3.*COS_sigEW[gdcos], nlim)
  if nlim GT 0 then begin
     cos_str[lim].mag[4] = 3. * COS_sigEW[gdcos[lim]]
     cos_str[lim].mag[5] = 100.
  endif
  red = where(megastruct[gdcos].galaxy.sfr_uplim eq 'yes', nred, complement = blue, ncomplement = nblue)
  cos_str[red].gal_type = 'Early'
  cos_str[blue].gal_type = 'Late'

  ;; Append
  if not keyword_set(NOCOS) then struct = [struct, cos_str]


  ;; Luminosity cuts
  lum_cuts = [ [0., 0.1], $
               [0.1, 1.], $
               [1., 1e9] ]
  sz = size(lum_cuts, /dimen)

  ;; Plot
  clr = getcolor(/load)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; EW of Lya
  xrng = [10., 1300.]
  yrng = [10., 4200] ; mA
  xtitle='!9r!X (kpc)'
  bclr = clr.black
  fclr = clr.lightgray
  
  plot,[0.],[0.],xrange=xrng,yrange=yrng,color=fclr,$
       background=bclr, charsize=csz,  xtickn=spcs, $
       POS=pos, xtitl=xtitle, $
       xmargin=xmrg, ymargin=ymrg, ytitle='W!uLy!9a!X!N (mA)',$
       /nodata,xstyle=1,ystyle=1, /xlog, /ylog

  all_dra = [0.]
  all_EW = [0.]
  val_dra = [0.]
  val_EW = [0.]
  flg = [0L]

  ;xyouts, 1000., 1300., '(a)', color=clr.black, align=0.5, charsi=lsz

  ;; Loop on Luminosity
  for qq=0L,sz[1]-1 do begin
     case qq of
        0: begin
           ppsym = 1
           lbl = '< 0.1L*'
        end
        1: begin
           ppsym = 2
           lbl = '(0.1-1) L*'
        end
        2: begin
           ppsym = 4
           lbl = '> L*'
        end
        else: stop
     endcase

     ;; Limit?
     lim = where(struct.mag[2] GT 0 AND struct.mag[5] GT 99. AND $
                 (struct.ddec LT lum_cuts[1,qq] and struct.ddec GT lum_cuts[0,qq]), nlim)
     plotsym, 1, 1.9, thick=5
     if nlim GT 0 then begin
        ;; Plot
        early = where(strmatch(strtrim(struct[lim].gal_type,2), 'Early'), nearly) 
        late = where(strmatch(strtrim(struct[lim].gal_type,2), 'Late'), nlate) 
        unkn = where(strmatch(strtrim(struct[lim].gal_type,2), 'Unkn'), nunkn) 

        if nearly GT 0 then oplot, [struct[lim[early]].dra], [struct[lim[early]].mag[4]], color=clr.tomato, psym=8
        if nlate GT 0 then oplot, [struct[lim[late]].dra], [struct[lim[late]].mag[4]], color=clr.cyan, psym=8
        if nunkn GT 0 then oplot, [struct[lim[unkn]].dra], [struct[lim[unkn]].mag[4]], color=clr.gray, psym=8

        all_dra = [all_dra, struct[lim].dra]
        all_EW = [all_EW, struct[lim].mag[4]]
        flg = [flg, replicate(2L,nlim)]
     endif

     ;; Value
     val = where(struct.mag[2] GT 0 AND struct.mag[5] LT 99. AND $
                 (struct.ddec LT lum_cuts[1,qq] and struct.ddec GT lum_cuts[0,qq]), nval)
     if nval GT 0 then begin
        val_dra = [val_dra, struct[val].dra]
        val_EW = [val_EW, struct[val].mag[4]]
        all_dra = [all_dra, struct[val].dra]
        all_EW = [all_EW, struct[val].mag[4]]
        flg = [flg, replicate(1L,nval)]

        early = where(strmatch(strtrim(struct[val].gal_type,2), 'Early'), nearly) 
        late = where(strmatch(strtrim(struct[val].gal_type,2), 'Late'), nlate) 
        unkn = where(strmatch(strtrim(struct[val].gal_type,2), 'Unkn'), nunkn) 

        if nearly GT 0 then oplot, [struct[val[early]].dra], [struct[val[early]].mag[4]], color=clr.tomato, psym=ppsym
        if nlate GT 0 then oplot, [struct[val[late]].dra], [struct[val[late]].mag[4]], color=clr.cyan, psym=ppsym
        if nunkn GT 0 then oplot, [struct[val[unkn]].dra], [struct[val[unkn]].mag[4]], color=clr.gray, psym=ppsym
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
     gal = where(struct.mag[2] GT 0 AND struct.dra GT med_bins[0,qq] AND $
                 struct.dra LE med_bins[1,qq], ngal)
     ymed[qq] = median(struct[gal].mag[4])
  endfor
  plotsym, 0, 1.3, thick=4
  oplot, xmed, ymed, color=clr.gray, psym=-8

  return

end
