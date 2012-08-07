;+ 
; NAME:
; fuse_velplt
;  V1.1
;
;   lbg_trans_ew, 1548.195, WINDEW=1., VMNX=[-300, 400]
;------------------------------------------------------------------------------
pro lbg_cii1334_ew, wrest, WINDEW=windEW, VMNX=vmnx

; lowzovi_prsdat -- Reads in DLA data to a structure

;if (N_params() LT 0) then begin 
;    print,'Syntax - ' + $
;      'fuse_velplt, strct, instr_list, vel_fil, NTOT=, CSIZE=, ' + $
;      'LSIZE=, PSFILE=, XTINT= (v1.1)' 
;    return
;endif 

;  
  if not keyword_set( NTOT ) then ntot = 12L
  if not keyword_set( LSZ ) then LSZ = 2.3
  if not keyword_set( LTHICK ) then lthick = 4.
  if not keyword_set( PTHICK ) then pthick = 2.
  if not keyword_set( BLSZ ) then blsz = 2.3
  if not keyword_set( CSIZE ) then csize = 2.5

  if not keyword_set( CONTI ) then conti = 0.8

  ;; Grab transition info
  wrest = 1334.5323
  getion, wrest, ioni, elm, Z=ionz, NM=strioni
  ;ionz = 6
  ;ioni = 2
  ;wrest = 1334.5323

  clambda = strtrim(fix(wrest),2)
  if not keyword_set(PSFILE) then $
     psfile = 'lbg_'+elm+strioni+clambda+'_ew.ps'

  x_psopen, psfile, /maxs
  clr = getcolor(/load)
  lclr = clr.white

  xmrg = [9,6]
  ymrg = [4.0,3.5]

  ;; Read in LBG
  lls_struct, lbg, '/u/xavier/LLS/Lists/lbg_lls.lst', ROOT='/u/xavier/LLS/', /ew

  ;; ;;;;;;;;;;;;;;;;;;
  ;; Histograms
  xrng = alog10([0.01, 5.]) ; Ang
  yrng = [0., 5] ; Number
  plot, [0], [0], xrange=xrng, xmar=xmrg, ymar=ymrg, $
        yrange=yrng, /NODATA, $
        charsize=csize, psym=10, background=clr.white, color=lclr, $
        xstyle=1, ystyle=1, thick=pthick, ytickinterval=5.,$
        xtitle='W!d'+clambda+'!N (Ang)', ytit='N', /noeras

  ;; All
  ion = where(abs(lbg.systems[0].ion[ionZ].state[ioni,*].lambda-wrest) LT 1e-3, nlbg)
  lbg_idx = ion / 30

  ;; R<100kpc
  lowR = where(lbg[lbg_idx].srvy_mag LT 100., complement=highR)
  LBG_EW = (lbg.systems[0].ion[ionZ].state[ioni,*].clm)[ion[lowR]]  > 0.01
  plothist, alog10(LBG_EW), bin=0.1, fcolor=clr.yellow, /fill, /overpl

  ;; R>100kpc
  LBG_EW = (lbg.systems[0].ion[ionZ].state[ioni,*].clm)[ion[highR]]  > 0.01
  plothist, alog10(LBG_EW), bin=0.1, color=clr.red, /overpl

  ;; Label
  xyouts, -1.7, 4., 'LBG/QSO Analysis', color=lclr, charsi=lsz, align=0.
  xyouts, 0., 4., 'R<100 kpc', color=clr.yellow, charsi=lsz, align=0.
  xyouts, 0., 3.5, 'R>100 kpc', color=clr.red, charsi=lsz, align=0.

  ;; Grid again
  plot, [0], [0], xrange=xrng, xmar=xmrg, ymar=ymrg, $
        yrange=yrng, /NODATA, $
        charsize=csize, psym=10, background=clr.white, color=lclr, $
        xstyle=1, ystyle=1, thick=pthick, ytickinterval=5.,$
        xtitle='W!d'+clambda+'!N (Ang)', ytit='N', /noeras

  if keyword_set( PSFILE ) then x_psclose

  !p.multi=[0,1,1]
  
  return
end
