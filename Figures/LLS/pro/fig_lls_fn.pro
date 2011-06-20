; .com ../Analysis/pro/dla_proxdx
; .com ../Analysis/pro/dla_proxstat
; fig_lox 

pro fig_lls_fn, PSFILE=psfile, GZFIL=gzfil, DLALST=dlalst, $
            RSDSS=rsdss, RTOT=rtot, TIMBIN=timbin, STRCT=strct, $
            BINS=bins, OUTFIL=outfil, VPROX=vprox, $
            VMAX=vmax, PLLS=plls, PROX=prox, QSOFIL=qsofil, $
            LLSFIL=llsfil, DR5_FNSTR=dr5_fnstr, DLA_STR=dla_str, $
            LX_SLLS=lx_slls, NOPLOT=noplot, BOOST=boost, thk=thk

  compile_opt strictarr

  ;; Get structure if necessary
  if not keyword_set( PSFILE ) then psfile = 'fig_lls_fn.ps'
  ;; Initialize
  cd, getenv('LLSPAP')+'/SDSS/DR7/Analysis/pro/', curr=curr
  RESOLVE_ROUTINE, 'max_lozpower'
  lls_dr7_initparm, init
  cd, curr

  if not keyword_set( LLSFIL ) then llsfil = init.llsfil
  if not keyword_set( QSOFIL ) then qsofil = init.qsofil 
  if not keyword_set( VPROX ) then vprox = init.vprox
  if not keyword_set( XZFIL ) then xzfil = init.xzfil
  if not keyword_set( BINS ) then BINS = init.bins
  if not keyword_set( MAXDZ ) then maxdz = init.maxoff 
  if not keyword_set( ZEM_MIN ) then zem_min = init.zem_min
  if not keyword_set( CSZ ) then csz = 2.1
  if not keyword_set( LSZ ) then lsz = 2.1
  if not keyword_set( STP ) then stp = 0.2
  if not keyword_set(NINT) then nint = 200L 
  if not keyword_set(BOOST) then boost = 1.

  if not keyword_set(THK) then thk = 7.

  if BOOST GT 1. then psfile = 'fig_boostfn.ps'

  ;; QSOs
  qsos = xmrdfits(qsofil, 1, /silent)
  nqso = n_elements(qsos)

  ;; LLS
  all_lls = xmrdfits(llsfil, 1)
  idx = sdss_llsstat(all_lls, qsos, VPROX=vprox, PARTIAL=PLLS, PROX=prox, $
                    MOCK=mock, MAXDZ=maxdz, ZEM_MIN=zem_min)
  lls = all_lls[idx]


  ;; DLA
  if not keyword_set( GZFIL ) then gzfil = '~/SDSS/DR5_QSO/dr5_dlagz_s2n4.fits'
  bins = [ [2.4,2.6], [3.4,4.0]]
  if not keyword_set(DLA_STR) then $
    sdss_dlalox, GZFIL=gzfil, STRCT=dla_str, BINS=bins
  if not keyword_set(DR5_FNSTR) then $
    sdss_fndla, GZFIL=GZFIL, STRCT=dr5_fnstr,  STP=0.2, XZFIL=xzfil, $
                /SALL, ALLDR=alldr, CL=0.95, NPLT=nplt, ZBIN=[3.4, 4.0]
  lX_DLA = dla_str.dndx[1]
  

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; PLOT
  if keyword_set( PSFILE ) and not keyword_set(NOPLOT) then begin
      x_psopen, psfile, /maxs
  endif
  clr = getcolor(/load)
  !p.multi=[0,1,1]
  lclr = clr.white
      
  xmrg = [8,2]
  ymrg = [4.0,1]
  if not keyword_set(YRNG) then yrng=[-23, -18.]
  if not keyword_set(XRNG) then xrng=[17.4, 21.]
  plot, [0], [0], color=lclr, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, xtitle='log N!dHI!N', $
        ytitle='log f(N!dHI!N,X)', yrange=yrng, thick=9, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata, $
        xthi=thk, ythi=thk

  ;;;;;;;;;;;;;;;;;;;;
  ;; DLA
  xval = dr5_fnstr.xval
;  xval[0] = 20.3
  yplt = dr5_fnstr.yplt
  bins = dr5_fnstr.bins
  yerr1 = dr5_fnstr.yerr1
  yerr2 = dr5_fnstr.yerr2
;  x_curvefill, xval, yerr1, yerr2, color=clr.red;, /line_fill

  val = where(yplt GT -99., complement=lim, ncomplement=nlim)
  oploterror, [xval[val]], [yplt[val]], [xval[val]-bins[val]],  $
              [yplt[val]-yerr2[val]], psym=3, thick=7, $
              color=clr.tomato, /lobar, errcolor=clr.tomato, errthi=6
  oploterror, [xval[val]], [yplt[val]], [bins[val]+stp-xval[val]],  $
              [yerr1[val]-yplt[val]], psym=3, thick=7, $
              color=clr.tomato, /hibar, errcolor=clr.tomato, errthi=6
  if nlim NE 0 then begin
      plotsym,1, 2., thick=7
      oploterror, [xval[lim]], [yerr1[lim]], [xval[lim]-bins[lim]], $
                  replicate(0.,nlim), $
                  color=clr.tomato, psym=8, /lobar, errcolor=clr.tomato, thick=7, errthi=6
      oploterror, [xval[lim]], [yerr1[lim]], [bins[lim]+stp-xval[lim]], $
                  replicate(0., nlim), $
                  color=clr.tomato, psym=8, /hibar, errcolor=clr.tomato, thick=7, errthi=6
  endif



  fplt = 20.3 + (22-20.3)*findgen(1000L)/1000.
  thy = fltarr(n_elements(fplt))
  ns10 = 10.^dr5_fnstr.Nd
  lw = where(fplt LT dr5_fnstr.Nd, complement=hi, ncomplement=numhi)
  thy[lw] = dr5_fnstr.k3 * (10^fplt[lw]/ns10)^dr5_fnstr.a3
  if numhi NE 0 then thy[hi] = dr5_fnstr.k3 * $
    (10^fplt[hi]/ns10)^dr5_fnstr.a4
  oplot, fplt, alog10(thy), color=lclr, thick=thk

  print, 'DLA double: ', dr5_fnstr.a3, dr5_fnstr.k3, dr5_fnstr.Nd

  xyouts, 20.40, -21.5, '!9b!X!dDLA!N ~ -2', color=clr.tomato, charsiz=lsz
  
  ;;;;;;;;;;;;;;;;;;;
  ;; SLLS
  if not keyword_set(lX_SLLS) then lX_SLLS = init.lX_SLLS 
  a_SLLS = -1.2    ;; Full sample
  N1 = 1d19
  N2 = 2d20
  k_SLLS = alog10(lX_SLLS * (1+a_SLLS) / (N2^(a_SLLS+1) - N1^(1+a_SLLS)))
  fplt = 19.0 + 1.3*findgen(1000)/1000.
  thy0 = 10.^k_SLLS * (10.^fplt)^a_SLLS

;  lX_SLLS = 10.^k_SLLS * (N2^(a_SLLS+1) - N1^(1+a_SLLS)) / $
;            (1+a_SLLS)
;  print, 'SLLS lX: ', lX_SLLS
  a1_SLLS = -1.4
  k1_SLLS = alog10(lX_SLLS * (1+a1_SLLS) / (N2^(a1_SLLS+1) - N1^(1+a1_SLLS)))
  thy1 = 10.^k1_SLLS * (10.^fplt)^a1_SLLS
;  oplot, fplt, alog10(thy1), color=clr.green, linesty=1
  a2_SLLS = -1.01
  k2_SLLS = alog10(lX_SLLS * (1+a2_SLLS) / (N2^(a2_SLLS+1) - N1^(1+a2_SLLS)))
  thy2 = 10.^k2_SLLS * (10.^fplt)^a2_SLLS
;  oplot, fplt, alog10(thy2), color=clr.green, linesty=2
  x_curvefill, [fplt[0],fplt[999]], $
               alog10([thy1[0],thy2[999]]), alog10([thy2[0],thy1[999]]), $
               color=clr.yellow
  oplot, fplt, alog10(thy0), color=lclr, thick=thk

  logfn_19 = alog10(10.^k_SLLS * (10.^19.)^a_SLLS)

  xyouts, 19.5, -20.3, '!9b!X!dSLLS!N = '+string(a_SLLS,format='(f4.1)')+$
          '!S!d-0.2!R!N!u+0.2!N', $
          color=clr.yellow, charsiz=lsz

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; LLS
  lX_LLS = 0.5 * BOOST    ;; Good for z=3.7; Error is 0.1
;  avec = -0.3 - 0.65*findgen(1000)/1000.
  Nt2 = 10.d^(17.5)
  N19 = 1d19
  avec = -0.10005 - 1.5*dindgen(100)/100.

  alls = fltarr(3,3)
  alls = fltarr(3,3)
  klls = fltarr(3,3)
  klls = fltarr(3,3)

  ;; MFP stuff
  nmfp = 15
  NHI_MFP = 15. + (2.5)*findgen(nmfp)/(nmfp-1)
  Nmin = 1d12
  Nmax = 1d22
  dlnN = alog(Nmax/Nmin) / (nint-1)
  Nval = Nmin * exp( dlnN * dindgen(nint) )

  print, 'logfn_19 = ', logfn_19
  ;;;;;;;
  ;; LLS constraint only
  for ii=1,2 do begin
      ;; Set SLLS and kvec
      if ii EQ 1 then f_SLLS = logfn_19-0.2 else f_SLLS = logfn_19+0.2
      kvec = f_SLLS - avec*19
      ;; Vary LX_LLS
      for jj=1,2 do begin
          if jj EQ 1 then LX_tLLS = 0.15 else LX_tLLS = 0.35 ;; 17.5-19 only
          ;; BOOST
          LX_tLLS = LX_tLLS + (lX_LLS - 0.5)
          ;;
          lvec = 10.^kvec / (1+avec) * ( N19^(1+avec) - Nt2^(1+avec) )
          mn = min(abs(lvec-LX_tLLS), imn)
          ;; Save 
          alls[ii,jj] = avec[imn]
          klls[ii,jj] = kvec[imn]

      endfor
  endfor

  ;;;;;;;;;
  ;; Plot
  npt = 101L
  fplt = 17.5 + 1.5*findgen(npt)/(npt-1)
  up_lls = replicate(-9999, npt)
  low_lls = replicate(1., npt)
  for ii=1,2 do begin
      for jj=1,2 do begin
          up_lls = up_lls > (kLLS[ii,jj] + fplt*aLLS[ii,jj]) 
          low_lls = low_lls < (kLLS[ii,jj] + fplt*aLLS[ii,jj]) 
      endfor
  endfor
  x_curvefill, fplt, up_lls, low_lls, color=clr.cyan
  for ii=1,2 do begin
      for jj=1,2 do begin
          oplot, fplt, kLLS[ii,jj] + fplt*aLLS[ii,jj], $ 
                 color=lclr,  linesty=ii, thick=thk
      endfor
  endfor

  ;; Central
  lX_LLS = 0.5 * BOOST   ;; Good for z=3.7; Error is 0.1
  kvec = k_SLLS + 19*(a_SLLS-avec)
  lvec = 10.^kvec / (1+avec) * ( N19^(1+avec) - Nt2^(1+avec) )
  lX_tLLS = LX_LLS - lX_SLLS - lX_DLA
  mn = min(abs(lvec-LX_tLLS), imn)
  ;; Best
  a_LLS = avec[imn]
  k_LLS = kvec[imn]
  fplt = 17.5 + 1.5*findgen(101)/100.
  oplot, fplt, k_LLS + fplt*a_LLS, color=lclr, thick=thk

  print, 'best: ', a_LLS, k_LLS 

  print, 'beta: ', a_LLS, alls

  xyouts, 18.1, -18.8, '!9b!X!dLLS!N = '+string(a_LLS,format='(f4.1)')+$
          '!S!d-'+string(abs(a_lls-min(alls[1:2,1:2],max=mxa))/2., $
                         format='(f4.1)')+'!R!N!u+'+$
          string(abs(a_lls-mxa)/2.,format='(f4.1)')+'!N', $
          color=clr.cyan, charsiz=lsz

  strct = { $
          a_LLS: a_lls, $
          k_LLS: k_lls, $
          a_SLLS: a_slls, $
          k_SLLS: k_slls, $
          dla_str: dr5_fnstr $
          }
          
  ;; Close Ps
  if keyword_set( PSFILE ) then x_psclose
  print, 'fig_fn:  All done!'
       
  return
end
      
