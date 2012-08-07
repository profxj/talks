;+ 
; NAME:
; fuse_velplt
;  V1.1
;------------------------------------------------------------------------------
pro qpq5_ew1334, wrest, FOREST=forest, NOPS=nops, CSZ=csz, NOLBG=nolbg, XMRG=xmrg

  if not keyword_set( LSZ ) then lsz = 1.8
  if not keyword_set( LTHICK ) then lthick = 4.
  if not keyword_set( PTHICK ) then pthick = 2.
  if not keyword_set( BLSZ ) then blsz = 1.7
  if not keyword_set( CSZ ) then CSZ = 2.5

  if not keyword_set( CONTI ) then conti = 0.8

  ;; Grab transition info
  if not keyword_set(WREST) then wrest = 1334.5323  ;; CII 1334
  getion, wrest, ioni, elm, Z=ionz, NM=strioni

  ;; QPQ7
  qpq_fil = '~/Dropbox/QSOPairs/qpq5_pairs.fits'
  qpq_strct = xmrdfits(qpq_fil, 1)

  if not keyword_set(FOREST) then begin ;; Cut out Lya forest
     good = where( (qpq_strct.z_fg+1)*wrest GT ((qpq_strct.z_bg+1)*1215.6701 + 20.), ngd)
     qpq_strct = qpq_strct[good]
  endif

  fwr = fix(wrest)
  cfwr = strtrim(fix(wrest),2)
  if not keyword_set(PSFILE) then psfile = 'fig_qpq5_ew1334.ps'


  ;; Plot
  if not keyword_set(NOPS) then begin 
     if keyword_set( PSFILE ) then x_psopen, psfile, /maxs
     !p.multi=[0,1,1]
     xmrg = [8, 1]
     ymrg = [4, 1]
  endif else begin
     ymrg = [0., 0]
     ;xspaces = replicate(' ',30) 
     ;xtit = ''
  endelse
  xtit='R!dphys!N (kpc)'
  clr = getcolor(/load)
  lbgc = clr.red
  lclr = clr.white
  ;lclr = clr.lightgray
  qclr = clr.yellow


  ;; ;;;;;;;;;;;;;;;;;;
  ;; Histograms
  yrng = [0.01, 5.] ; Ang
  xrng = [-10., 310]
  plot, [0], [0], xrange=xrng, $
        yrange=yrng, /NODATA, xmarg=xmrg, ymarg=ymrg, $
        charsize=CSZ, psym=10, background=clr.white, color=lclr, $
        xstyle=1, ystyle=1, thick=pthick, ytickinterval=5., xtickn=xspaces, $
        ytitle='W!d'+cfwr+'!N (Ang)', xtit=xtit, /noeras, /ylog

  ;if keyword_set(NOPS) then xyouts, 01., 2.2, '(b)', color=lclr, charsi=lsz, align=0.

  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; QPQ7
  gd = where( abs(qpq_strct.metal_wrest-wrest) LT 1e-3 AND $
              qpq_strct.flg_metal_EW EQ 1)
  sz = size(qpq_strct.metal_wrest,/dimen)
  gdi = gd / sz[0]
  QPQ7_EW = (qpq_strct.metal_EW)[gd] > 0.013
  plotsym, 8, 1.4, /fill
  oplot, qpq_strct[gdi].R_phys, QPQ7_EW, color=qclr, psym=8

  ;; ;;;;;;;
  ;; Limits
  ulim = where( abs(qpq_strct.metal_wrest-wrest) LT 1e-3 AND $
                qpq_strct.flg_metal_EW EQ 3)
  ui = ulim / sz[0]
  QPQ7_EW = (qpq_strct.metal_EW)[ulim] 
  QPQ7_sigEW = (qpq_strct.metal_sigEW)[ulim] 
  plotsym, 8, 1.4, thick=4
  oplot, qpq_strct[ui].R_phys, QPQ7_EW > 0.013, color=qclr, psym=8
  plotsym, 1, 1.4, thick=4
  oplot, qpq_strct[ui].R_phys, QPQ7_EW, color=lclr, psym=8

  ;; 2 sigma??
  tsig = where( QPQ7_EW GT 2.*QPQ7_sigEW, ntsig) 
  if ntsig GT 0 then oplot, [qpq_strct[ui[tsig]].R_phys], [QPQ7_EW[tsig]], color=clr.lightgray, psym=2

  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; LBG

  ;; Read in LBG
  lls_struct, lbg, '/u/xavier/LLS/Lists/lbg_lls.lst', ROOT='/u/xavier/LLS/', /ew
  ion = where(abs(lbg.systems[0].ion[ionZ].state[ioni,*].lambda-wrest) LT 1e-3, nlbg)
  lbg_idx = ion / 30

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

  if not keyword_set(NOPS) then begin 
     if keyword_set( PSFILE ) then x_psclose
     !p.multi=[0,1,1]
  endif
  
  return
end
