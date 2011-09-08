;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro fig_hist_dv_lyaovi, summ_fil, STRUCT=struct, PSFILE=psfile

  if not keyword_set(psfile) then psfile = 'fig_hist_dv_lyaovi.ps'
  if not keyword_set(binsz) then binsz = 0.01
  if not keyword_set(lsz) then lsz = 2.3 
  if not keyword_set(csz) then csz = 2.3 
  if not keyword_set(dlim) then dlim = 5. ;arcmin
  if not keyword_set(RHO_GAL) then rho_gal = 200. ; kpc

  ;; Input
  if not keyword_set(SUMM_FIL) then summ_fil = '~/paper/OVI/Galaxies/Analysis/all_galabs_1Mpc_strct.fits'
  struct = xmrdfits(summ_fil,1)
  ngal = n_elements(struct)

  ;; Plot
  if keyword_set(psfile) then x_psopen,psfile,/maxs
  !p.multi=[0,1,1]
  clr = getcolor(/load)
  lclr = clr.white

  ;;;;;;;;;;
  yrng = [0., 10]
  xrng = [-400., 400]  ; km/s
  xmrg = [8,1.5]
  ymrg = [4,0.5]
  
  ;;
  plot,[0.],[0.],xrange=xrng,yrange=yrng,color=lclr,$
       background=clr.white, charsize=csz, $
       xmargin=xmrg, ymargin=ymrg,xtitle='!9d!Xv (km s!u-1!N)', $
       ytitle='Number', /nodata,xstyle=1,ystyle=1;, /xlog

  ;; HI
  val = where(struct.mag[2] GT 0 AND $
              struct.dra LT 300. AND $
              (struct.mag[5] LT 99. OR $
               (struct.mag[9] LT 9. and $
                struct.mag[9] GT 0)), nval) 
  plothist, x_relvel(struct[val].magerr[3],struct[val].z,/rev), bin=40., $
                           color=lclr, /overplo, fcolor=lclr

  ;; OVI
  val = where(struct.magerr[2] GT 0 AND $
              struct.dra LT 300. AND $
              (struct.magerr[5] LT 99. OR $
               (struct.magerr[9] LT 9. and $
                struct.magerr[9] GT 0)), nval) 
  plothist, x_relvel(struct[val].magerr[3],struct[val].z,/rev), bin=40., $
                           color=clr.cyan, /overplo, fcolor=clr.cyan, /fill

  xyouts, -320., 8.7, 'Ly!9a!X', color=lclr, charsi=lsz
  xyouts, -320., 8., 'OVI', color=clr.cyan, charsi=lsz

  plot,[0.],[0.],xrange=xrng,yrange=yrng,color=lclr,$
       background=clr.white, charsize=csz, $
       xmargin=xmrg, ymargin=ymrg,xtitle='!9d!Xv (km s!u-1!N)', $
       ytitle='Number', /nodata,xstyle=1,ystyle=1, /noeras ;, /xlog

  x_psclose
  !p.multi=[0,1,1]
  return

end
