;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; fig_cheme, psfile='Figures/fig_cheme.ps', DLA=dla
pro fig_hist1526, PSFILE=psfile, NBIN=nbin, DLA=dla, GRB=grb, LBL=lbl, $
                SZLBL=lblsz, CSIZE=csize, CUMUL=cumul, TALK=talk

  psfile = 'fig_hist1526.ps'
  if not keyword_set( NBIN ) then nbin = 5
  if not keyword_set(CSIZE) then csize = 2.3
  if not keyword_set(LBLSZ) then lblsz = 1.9

  ;; Parse DLA files
  if not keyword_set(DLA) then parse_dlalst, dla, $
    '/u/xavier/DLA/Lists/metal704_dla.lst', /ew, ROOT='/u/xavier/DLA/'
  ;; Parse GRB files
  if not keyword_set(GRB) then parse_dlalst, grb, $
    '/u/xavier/DLA/Lists/grb_dla.lst', /ew, ROOT=getenv('DLA')

  ;; Plot
  gdd = where(dla.NHI GE 20.3, ndla)
  srt = sort(dla[gdd].zabs)
  dla = dla[gdd[srt]]
;  ngrb = n_elements(grb)


  ;; Plot
  if keyword_set( PSFILE ) then begin
      x_psopen, psfile, /maxs
      !p.multi=[0,1,1]
  endif
  clr = getcolor(/load)
  grb_clr=clr.green

  ;; Overall
  yrng = [0., 10.]
  xrng = [0, 2.7]
  xmrg = [7,1]
  plot, [0], [0], color=clr.lightgray, $
    background=clr.black, charsize=csize,$
    xmargin=xmrg, ymargin=[3.5,0.5], ytitle='Number', $
        xtitle='W!d1526!N (Ang)', /nodata, xthick=7, $
        ythick=7, xstyle=1, ystyle=1,  yr=yrng, xr=xrng

  ;; DLA
  gd = where(abs(dla.ion[14].state[2,*].lambda-1526.) LT 1.)
  plothist, (dla.ion[14].state[2,*].clm)[gd], bin=0.03, /overplot, $
            color=grb_clr, fcolor=grb_clr, /halfbin, thick=7
  qval = (dla.ion[14].state[2,*].clm)[gd]
  writecol, 'dla_zew.dat', (dla.zabs)[gd], $
            (dla.ion[14].state[2,*].clm)[gd] 

  ;; GRB
  gdg = where(abs(grb.ion[14].state[2,*].lambda-1526.) LT 1.)
  plothist, (grb.ion[14].state[2,*].clm)[gdg], bin=0.2, /overplot, $
            color=clr.cyan, fcolor=clr.cyan, /halfbin, thick=7;, /fill
  gval = (grb.ion[14].state[2,*].clm)[gdg]


  ;;;;;;;;;;;;;;;;
  ;; Label
  xyouts, xrng[0]+0.1, yrng[1]-2, 'GRB-DLA', color=clr.cyan, charsiz=lblsz 
  xyouts, xrng[0]+0.1, yrng[1]-2.7, 'QSO-DLA', color=grb_clr, charsiz=lblsz 

  kstwo, (dla.ion[14].state[2,*].clm)[gd] ,(grb.ion[14].state[2,*].clm)[gdg], $
         d, prob
  print, 'ks = ', prob

;  xyouts, xrng[1]-0.7, yrng[1]-5, 'P!dKS!N = '+string(prob,format='(f6.4)'), $
;          color=clr.black, charsiz=lblsz 

  if keyword_set(LBL) then $
    xyouts, xrng[0]+0.05, yrng[1]-1, LBL, color=clr.black, charsiz=lblsz

  ;; Cumulative
  if keyword_set(CUMUL) then begin
      if not keyword_set(POS) then pos = [0.5, 0.27, 0.95, 0.47]
      plot, [0.], [0.], color=clr.black, $
            background=clr.white, charsize=1.2,$
            xmargin=[7,1], ymargin=[3.5,0.5], xtitle='W!d1526!N (Ang)', $
            ytitle='Cumul Fraction', /nodata, xthick=5, $
            ythick=5, xstyle=1, ystyle=1, pos=pos, $
            yr=[0., 1.], xr=xrng, /noerase

      ;; GRB
      srt = sort(gval)
      ngrb = n_elements(gval)
      oplot, [0.,gval[srt]], (lindgen(ngrb+1))/float(ngrb), color=clr.cyan

      ;; QSO
      srt = sort(qval)
      nqdla = n_elements(qval)
      oplot, [0.,qval[srt]], (lindgen(nqdla+1))/float(nqdla), color=grb_clr
  endif
      
  if keyword_set( PSFILE ) then begin
      x_psclose
      !p.multi=[0,1,1]
  endif
  

  return
end
