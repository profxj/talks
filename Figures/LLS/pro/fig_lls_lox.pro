; .com ../Analysis/pro/dla_proxdx
; .com ../Analysis/pro/dla_proxstat
; fig_lox 

pro fig_lls_lox, PSFILE=psfile, GZFIL=gzfil, DLALST=dlalst, $
             RSDSS=rsdss, RTOT=rtot, TIMBIN=timbin, STRCT=strct, $
             BINS=bins, OUTFIL=outfil, VPROX=vprox, $
             VMAX=vmax, PLLS=plls, PROX=prox, QSOFIL=qsofil, $
             LLSFIL=llsfil

  compile_opt strictarr

  ;; Initialize
  cd, getenv('LLSPAP')+'/SDSS/DR7/Analysis/pro/', curr=curr
  RESOLVE_ROUTINE, 'max_lozpower'
  lls_dr7_initparm, init
  cd, curr

  if not keyword_set( LLSFIL ) then llsfil = init.llsfil
  if not keyword_set( QSOFIL ) then qsofil = init.qsofil  else init.qsofil = qsofil
  if not keyword_set( MAXDZ ) then maxdz = init.maxoff 
  if not keyword_set( VPROX ) then vprox = init.vprox
  if not keyword_set( XZFIL ) then xzfil = init.xzfil
  if not keyword_set( ZEM_MIN ) then zem_min = init.zem_min
  if not keyword_set( BINS ) then BINS = init.all_bins

  ;; Calculate the FIT
  max_lozpower, STRCT=fit_strct, Z0=z0, INIT=init, MAXDZ=maxdz, LLS=lls
  nlls = n_elements(LLS)

  ;; Get structure if necessary
  if not keyword_set( PSFILE ) then psfile = 'fig_lls_lox.ps'
  if not keyword_set( CSZ ) then csz = 2.1
  if not keyword_set( LSZ ) then lsz = 2.1

  ;; Calculate from the data
  sdss_llslozx, qsofil, llsfil, strct, BINS=bins, XZFIL=xzfil, VPROX=vprox, $
                ZEM_MIN=zem_min, MAXDZ=maxdz

  sdss_llslozx, qsofil, llsfil, lowz_strct, BINS=bins[*,0:1], $
                XZFIL=xzfil, VPROX=vprox, MAXDZ=maxdz ;; No zem cut



  ;; PLOT
  if keyword_set( PSFILE ) then x_psopen, psfile, /maxs
  clr = getcolor(/load)
  lclr = clr.white
  !p.multi=[0,1,1]
  
  xtit = 'z' 
  yrng = [0.0, 2.5]
  xrng = [3.2, 6.1]
  for qq=0,1 do begin
          
     Sscl = 1.1
     ;;
     plot, [0.], [0.], color=lclr, background=clr.white, charsize=csz,$
           xmargin=[8.0,1], ymargin=[4,1.], xtitle=xtit, $
           /nodata, xrange=xrng, ystyle=1, $
           yrange=yrng, xstyle=1, xtickn=spaces

  ;; Label
;  xyouts, 2.2, yrng[1]*0.85, lbl, color=clr.black, chars=1.8
  
     xlbl = 0.05
;  xlbl = 0.08
     !p.font = -1
     xyouts, xlbl, 0.50, '!12l!X', alignment=0.0, color=lclr,$
;  xyouts, xlbl, 0.75, '!12l!X', alignment=0.0, color=clr.black,$
             charsize=2.6, /normal, orientation=90
     !p.font = 0
     xyouts, xlbl, 0.51, '!d!9t!X>=1!N(X)', alignment=0.0, color=lclr,$
;  xyouts, xlbl, 0.76, '!dLLS!N(z)', alignment=0.0, color=clr.black,$
             charsize=2.2, /normal, orientation=90
     
     ;; Plot 'good' points
     good = where(bins[0,*] GE 3.5, complement=bad)
     oploterror, strct[good].medz, strct[good].lox*Sscl, $
                 strct[good].medz-bins[0,good], strct[good].siglx[0], $
                 psym=1, color=lclr, errcolor=lclr, /lobar
     oploterror, strct[good].medz, strct[good].lox*Sscl, $
                 bins[1,good]-strct[good].medz, strct[good].siglx[1], $
                 psym=1, color=lclr, errcolor=lclr, /hibar
     printcol, strct[good].medz, strct[good].lox

     ;; Plot low z points
     oploterror, strct[bad].medz, strct[bad].lox*Sscl, $
                 strct[bad].medz-bins[0,bad], strct[bad].siglx[0], $
                 psym=1, color=clr.gray, errcolor=clr.gray, /lobar
     oploterror, strct[bad].medz, strct[bad].lox*Sscl, $
                 bins[1,bad]-strct[bad].medz, strct[bad].siglx[1], $
                 psym=1, color=clr.gray, errcolor=clr.gray, /hibar

     ;; Label
     xyouts, 5.0, 0.1, 'SDSS-DR7: POW10', color=lclr, charsi=lsz, align=0.
     if qq EQ 1 then begin ;; Songaila
        print, 'Show SC06!!'
        z_dndz = [ [4.2, 5.], $
                   [5., 6.] ]
        mz_dndz = [4.4, 5.7]
        dndz = [4.04, 8.91]
        sig_dndz = [0.7, 3.5]
        for jj=0,1 do begin
           dxdz = (cosm_xz(mz_dndz[jj]+0.05, /w05map) - $
                  cosm_xz(mz_dndz[jj]-0.05, /w05map))  / 0.1
           oploterror, [mz_dndz[jj]],  [dndz[jj]]/dxdz, $
                       [mz_dndz[jj]]-z_dndz[0,jj], sig_dndz[jj]/dxdz, $
                       psym=1, color=clr.cyan, errcolor=clr.cyan, /lobar
           oploterror, [mz_dndz[jj]],  [dndz[jj]]/dxdz, $
                       z_dndz[1,jj]-[mz_dndz[jj]], sig_dndz[jj]/dxdz, $
                       psym=1, color=clr.cyan, errcolor=clr.cyan, /hibar
;           if jj EQ 1 then stop
        endfor
        xyouts, 5.0, 0.3, 'Keck/ESI: SC11', color=clr.cyan, charsi=lsz, align=0.
     endif

  endfor
     
  ;; Close Ps
  if keyword_set( PSFILE ) then x_psclose
  print, 'fig_lox:  All done!'

  return
end
      
      
