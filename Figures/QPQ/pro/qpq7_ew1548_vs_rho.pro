;+ 
; NAME:
; fuse_velplt
;  V1.1
;
;   lbg_trans_ew, 1548.195, WINDEW=1., VMNX=[-300, 400]
;------------------------------------------------------------------------------
pro qpq7_ew1548_vs_rho, wrest, FOREST=forest

  if not keyword_set( LSZ ) then lsz = 1.8
  if not keyword_set( LTHICK ) then lthick = 4.
  if not keyword_set( PTHICK ) then pthick = 2.
  if not keyword_set( BLSZ ) then blsz = 1.7
  if not keyword_set( CSIZE ) then csize = 2.1

  if not keyword_set( CONTI ) then conti = 0.8

  ;; Grab transition info
  if not keyword_set(WREST) then wrest = 1548.195d ;
  getion, wrest, ioni, elm, Z=ionz, NM=strioni

  ;; QPQ7
  qpq_fil = '~/Dropbox/QSOPairs/qpq7_pairs.fits'
  qpq_strct = xmrdfits(qpq_fil, 1)

  if not keyword_set(FOREST) then begin ;; Cut out Lya forest
     good = where( (qpq_strct.z_fg+1)*wrest GT ((qpq_strct.z_bg+1)*1215.6701 + 20.), ngd)
     qpq_strct = qpq_strct[good]
  endif

  fwr = fix(wrest)
  cfwr = strtrim(fix(wrest),2)
  if not keyword_set(PSFILE) then psfile = 'qpq7_ew1548_vs_rho.ps'

  ;; Read in LBG
  if keyword_set(LBG) then begin
     lls_struct, lbg, '/u/xavier/LLS/Lists/lbg_lls.lst', ROOT='/u/xavier/LLS/', /ew
     ion = where(abs(lbg.systems[0].ion[ionZ].state[ioni,*].lambda-wrest) LT 1e-3, nlbg)
     lbg_idx = ion / 30
     nobj = nlbg  + 2 ;; Stacked spectra
  endif
  

  ;; Plot
  if keyword_set( PSFILE ) then x_psopen, psfile, /maxs
  !p.multi=[0,1,1]
  clr = getcolor(/load)
  lbgc = clr.red
  fclr = clr.white
  ;fclr = clr.lightgray


  ;; ;;;;;;;;;;;;;;;;;;
  ;; Histograms
  yrng = [0.01, 5.] ; Ang
  xrng = [30., 1000] ; Number
  xmrg = [8, 13]
  ymrg = [4, 1]
  thisletter = byte(94)
  perpletter = '!9' + string(thisletter) + '!X'
  xtit = 'R!d'+perpletter+'!N (kpc)'
  plot, [0], [0], xrange=xrng, $
        yrange=yrng, /NODATA, xmarg=xmrg, ymarg=ymrg, $
        charsize=csize, psym=10, background=clr.white, color=fclr, $
        xstyle=1, ystyle=1, thick=pthick, ytickinterval=5.,$
        ytitle='W!d'+cfwr+'!N (Ang)', xtit=xtit, /noeras, /ylog, /xlog

  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; QPQ7
  gd = where( abs(qpq_strct.metal_wrest-wrest) LT 1e-3 AND $
              qpq_strct.flg_metal_EW EQ 1, ngd)
  print, 'Number of measurements = ', ngd
  sz = size(qpq_strct.metal_wrest,/dimen)
  gdi = gd / sz[0]
  QPQ7_EW = (qpq_strct.metal_EW)[gd] > 0.013
  plotsym, 8, 1.4, /fill
  oplot, qpq_strct[gdi].R_phys, QPQ7_EW, color=fclr, psym=8

  ;; ;;;;;;;
  ;; Limits
  ulim = where( abs(qpq_strct.metal_wrest-wrest) LT 1e-3 AND $
                qpq_strct.flg_metal_EW EQ 3, nulim)
  print, 'Number of upper limits = ', nulim
  ui = ulim / sz[0]
  QPQ7_EW = (qpq_strct.metal_EW)[ulim] 
  QPQ7_sigEW = (qpq_strct.metal_sigEW)[ulim] 
  plotsym, 8, 1.4
  oplot, qpq_strct[ui].R_phys, QPQ7_EW > 0.013, color=fclr, psym=8
  plotsym, 1, 1.4
  oplot, qpq_strct[ui].R_phys, QPQ7_EW, color=fclr, psym=8

  ;; 2 sigma??
  tsig = where( QPQ7_EW GT 2.*QPQ7_sigEW, ntsig) 
  if ntsig GT 0 then oplot, [qpq_strct[ui[tsig]].R_phys], [QPQ7_EW[tsig]], color=clr.gray, psym=2

  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; LBG
  if keyword_set(LBG) then begin
     lbgc = clr.darkgreen

     ;; EW
     LBG_EW = (lbg.systems[0].ion[ionZ].state[ioni,*].clm)[ion]  
     LBG_sigEW = (lbg.systems[0].ion[ionZ].state[ioni,*].sigclm)[ion]  
     lim = where(LBG_EW LT 3*LBG_sigEW OR LBG_EW LT 0.02, complement=gdval)
     
     LBG_R = lbg[lbg_idx]. srvy_mag
     if not  keyword_set(NOLBG) then begin
        oplot, LBG_R[gdval], LBG_EW[gdval], color=lbgc, psym=1
        plotsym, 1, 1.4, thick=5
        oplot, LBG_R[lim], LBG_EW[lim] > 0.02, color=lbgc, psym=8
     endif
  endif

  if keyword_set( PSFILE ) then x_psclose

  !p.multi=[0,1,1]
  
  return
end
