;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; fig_cheme, psfile='Figures/fig_cheme.ps', DLA=dla
pro fig_histnhi, PSFILE=psfile, NBIN=nbin, DLA=dla, GRB=grb, LBL=lbl, $
                SZLBL=lblsz, CSIZE=csize, CUMUL=cumul, TALK=talk

  psfile = 'fig_histnhi.ps'
  if not keyword_set( NBIN ) then nbin = 5
  if not keyword_set(CSIZE) then csize = 2.3
  if not keyword_set(LSZ) then lsz = 2.5

  ;; Parse GRB-DLA files
  if not keyword_set(GRB) then parse_dlalst, grb_dla, $
    '/u/xavier/DLA/Lists/grb_dla.lst', ROOT=getenv('DLA')
  if not keyword_set(GRB) then lls_struct, grb_lls, $
    '/u/xavier/LLS/Lists/grb_lls.lst', ROOT=getenv('LLSTREE')

  x_psopen, psfile, /maxs
  !p.multi=[0,1,1]
  clr = getcolor(/load)
  grb_clr=clr.green

  ;; Overall
  yrng = [0., 20.]
  xrng = [16.5, 23.]
  xmrg = [7,1]
  ymrg = [6,7]

  for qq=0,4 do begin

     plot, [0], [0], color=clr.lightgray, $
           background=clr.black, charsize=csize,$
           xmargin=xmrg, ymargin=ymrg, ytitle='Number', $
           xtitle='log N!dHI!N', /nodata, xthick=7, $
           ythick=7, xstyle=1, ystyle=1,  yr=yrng, xr=xrng

     ;; Histogram
     plothist, [grb_dla.nhi,grb_lls.nhi], bin=0.5, color=clr.cyan, $
               /fill, fcolor=clr.cyan, /overpl

     ;; GMC
     if qq GT 0 then begin
        xyouts, 22.5, 17., '!9S!X!dGMC!N', color=clr.yellow, charsiz=lsz, align=0.5
        oplot, [21.8, 21.8], yrng, color=clr.lightgray, linesty=2, thick=6
     endif

     ;; HI disk
     if qq GT 1 then begin
        xyouts, 20.9, 17., 'HI disk', color=clr.green, charsiz=lsz, align=0.5
        oplot, [20.0, 20.0], yrng, color=clr.lightgray, linesty=2, thick=6
     endif

     ;; Halo disk
     if qq GT 2 then begin
        xyouts, 18.7, 17., 'Halo?', color=clr.orange, charsiz=lsz, align=0.5
        oplot, [17.5, 17.5], yrng, color=clr.lightgray, linesty=2, thick=6
     endif

     ;; Optically thin
     if qq GT 3 then begin
        xyouts, 17.3, 5., 'Optically thin!!', color=clr.pink, charsiz=lsz, align=0.0, $
                orient=90.
     endif
  endfor

  x_psclose
  !p.multi=[0,1,1]
  

  return
end
