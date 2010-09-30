pro hist_stars, RREAL=rreal

  if not keyword_set( PSFILE ) then psfile = 'hist_stars.ps'
  if not keyword_set(CSZ) then csz = 1.9
  if not keyword_set(lSZ) then lsz = 1.9
  if not keyword_set(THK) then thk = 11

  ;; Mass interval
  npt = 50L
  m0 = 0.01 * 10.^(4*findgen(npt)/(npt-1))
  m1 = m0[1:*]
  m0 = m0[0:npt-2]

  beta = -1.35
  nstar = m0^beta - m1^beta

  ;;; BEGIN PLOTS
  x_psopen, psfile, /maxs
  clr = getcolor(/load)

  yrng=[1., 3e3]
  xrng=[1e-2, 1.5e2]

  for qq=0,4 do begin

     if qq GE 4 then yrng=[1., 3e5]
  ;; Plot
  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        xtitle='Mass of Star', ymarg=[4,1], xmarg=[7,1], $
        ytitle='Number', yrange=yrng, thick=thk, xthic=thk, ythic=thk, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata, /xlog, /ylog
     

  ;; Sun
  mn = min( abs([m0+m1]/2.- 1.), isun)
  nsun = nstar[isun]
  nrm_sun = 1000.
  if qq LT 3 then $
     x_curvefill, [m0[isun], m1[isun]], $
                  [1e-10, 1e-10], [nrm_sun,nrm_sun], color=clr.black
  xyouts, 1., nrm_sun*1.5, 'Sun', color=clr.blue, align=0.5, charsi=lsz  

  ;; 10x Sun
  if qq GE 1 then begin
     mn = min( abs([m0+m1]/2.- 10.), i10)
     n10 = nstar[i10]
     nrm_sun = 1000.
     yval = nrm_sun*n10/nsun
     if qq LT 3 then begin 
        x_curvefill, [m0[i10], m1[i10]], $
                     [1e-10, 1e-10], [yval, yval], color=clr.black
        xyouts, 10., yval*1.5, '10x Sun', color=clr.red, align=0.5, charsi=lsz  
     endif
  endif

  ;; 100x Sun
  if qq GE 2 then begin
     mn = min( abs([m0+m1]/2.- 100.), i100)
     n100 = nstar[i100]
     nrm_sun = 1000.
     yval = nrm_sun*n100/nsun
     if qq LT 3 then begin 
        x_curvefill, [m0[i100], m1[i100]], $
                     [1e-10, 1e-10], [yval, yval], $
                     color=clr.black
        xyouts, 70., yval*1.5, '100x Sun', color=clr.red, align=0.5, charsi=lsz  
     endif
  endif

  if qq GE 3 then begin
     a = where(m1 GE 1., na)
     xplt = [m0[a], m1[a]]
     yplt = [nstar[a],nstar[a]]
     srt = sort(xplt)
     xplt = xplt[srt]
     yplt = yplt[srt]
     x_curvefill, m0[a], nstar[a]*nrm_sun/nsun, replicate(1e-10,na), $
                  color=clr.black
     xyouts, 20., 350.,  'Massive Stars', $
             color=clr.red, align=0.5, charsi=lsz  
     xyouts, 20., 200.,  'Rare, but dominate the light', $
             color=clr.red, align=0.5, charsi=lsz  
  endif

  if qq GE 4 then begin
     a = where(m0 LE 1., na)
     oplot, m0[a], nstar[a]*nrm_sun/nsun, color=clr.black, linesty=1, thick=thk
     xyouts, 0.08, 350.,  'Low Mass Stars', $
             color=clr.cyan, align=0.5, charsi=lsz  
     xyouts, 0.08, 200.,  'These have most of the mass!', $
             color=clr.cyan, align=0.5, charsi=lsz  
  endif
endfor

  if keyword_set( PSFILE ) then x_psclose
  !p.multi = [0,1,1]

  return

end
