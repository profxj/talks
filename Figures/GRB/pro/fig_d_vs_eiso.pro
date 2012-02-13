;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; fig_cheme, psfile='Figures/fig_cheme.ps', DLA=dla
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
pro fig_d_vs_eiso, PSFILE=psfile, NBIN=nbin, SURVIVE=SURVIVE, $
              ZMAX=zmax, ZMIN=zmin, NODLA=nodla, LBG=lbg, $
              NARROW=narrow, NONO=nono, XRNG=xrng, BW=bw, TWOPG=twopg, $
              CSZ=csz, LSZ=lsz

  if not keyword_set(PSFILE) then psfile = 'fig_d_vs_esio.ps'
  if not keyword_set( CSZ ) then csz = 2.0
  if not keyword_set( LSZ ) then lsz = 1.8

  ;; Data
  GRBnm = ['020813', '050730', '051111', '060418', '080310', '080319B', $
           '080330', '090426', '090926']
  dgas = [75., 125., 200., 480, 160, 1000, 80, 80, 680]
  sig_dgas = [25., 20, 100., 60, 20, 200, 20, 20, 40]
  ngrb = n_elements(grbnm)
  Eiso = [76., 9., 6.2, 10., 3.4, 141., 0.41, 2.6, 189]

  ;; Plot
  if keyword_set( PSFILE ) then x_psopen, psfile, /maxs
  clr = getcolor(/load)
  !p.multi=[0,1,1] 


  ymrg = [4,1]
  xrng = [0.1, 300]
  yrng = [1., 1e4]
  ;; Overall
  plot, [0], [0],  color=clr.black, $
    background=clr.black, charsize=csz,$
    xmargin=[9,1], ymargin=ymrg, xtitle='E!diso!N (10!u52!N)', $
    ytitle='d!dgas!N (pc)', /nodata, xthick=5, ythick=5, xstyle=1, ystyle=1, $
    yr=yrng, xr=xrng, /ylog, /xlog

  plotsym, 0, 1.5, /fill
  oploterror, [Eiso], [dgas], [sig_dgas], color=clr.blue, psym=8, errcolo=clr.blue


  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]

  close, /all

  return
end
