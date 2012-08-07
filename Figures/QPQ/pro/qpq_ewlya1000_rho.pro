;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  Generates histograms of Lya EW vs. Rho
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro qpq_ewlya1000_rho, STRCT=strct,  NOSTOP=nostop

  compile_opt strictarr

  c=x_constants()
  cd, getenv('QSO_DIR')+'tex/QPQ6/Analysis/AbsLines/pro/', curr=curr
  resolve_routine, 'qpq6_mkdata', /compile_full_file, /is_function
  resolve_routine, 'qpq6_grab_idlines', /compile_full_file, /is_function
  cd, curr
  cd, getenv('QSO_DIR')+'tex/QPQ6/Analysis/pro/', curr=curr
  resolve_routine, 'qpq6_anly_cut', /compile_full_file, /is_function
  resolve_routine, 'qpq6_get_instrument', /compile_full_file, /is_function
  cd, curr

  ;; Get structure if necessary
  if not keyword_set( PSFILE ) then psfile = 'qpq_ewlya1000_rho.ps'
  if not keyword_set( CSZ ) then csz = 2.3
  if not keyword_set( LSZ ) then lsz = 1.5

  if not keyword_set(SIG_FX) then sig_fx = 0.05 ;; Error in tau values
  if not keyword_set( MAXCHI ) then MAXCHI = 1.3
  if not keyword_set( NRNG ) then nrng = 50L
  if not keyword_set(VINT) then VINT = 1000.

  qsodir = getenv('QSO_DIR')

  ;; Read in pair structure
  qpq_fil = '~/Dropbox/QSOPairs/qpq6_pairs.fits'
  qpq_strct = xmrdfits(qpq_fil, 1)

  gdp = qpq6_anly_cut(qpq_strct, NPAIRS=npair)
  gd_strct = qpq_strct[gdp]

  ;; ;;;;;;;;;;;;;;;;;;;;;;;
  ;; Data
  ewi = round( median( where(abs(gd_strct.vEW - VINT) LT 1.)  $
                       MOD n_elements(gd_strct[0].vEW) ) )
  all_EWLya = gd_strct.vEWLya[ewi]
  all_Rphys = gd_strct.R_phys

  Rcuts = [ $
          [0, 300], $
          [0, 100], $
          [100, 200], $
          [200, 300] ]
  sz_R = size(Rcuts, /dimen)

  npx = 1
  npy = sz_R[1]
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; PLOT
  x_psopen, psfile, /portrait
  !p.multi=[npx*(npy+1),npx,npy+1,0,1]
  clr = getcolor(/load)
  ;lclr = clr.lightgray
  lclr = clr.white

  ;; LOOP
  for tt=0L,npy-1 do begin
     if (tt mod npy) EQ 0 then begin
        if tt GT 0 then offy = (nVINT MOD 2) else offy = 0
        !p.multi=[npx*(npy+1) - (npy+1)*(tt/npy) - offy, $
                  npx,npy+1,0,1]
     endif

     xmrg = [8,11]
     ymrg = [0,0.]

     xrng = [1e-1, 20]
     yrng = [0., 24]

     ytit='N'

     ;; Plot
     ysty=1
     if (tt MOD npy) NE (npy-1) and (tt NE npy-1) then begin
        xspaces = replicate(' ',30) 
        xtitle = ''
     endif else begin 
        if keyword_set(XSPACES) then delvarx, xspaces
        xtit='W!dLy!9a!X!N (Ang)'
     endelse
     if (tt GT npy-1) then yspaces = replicate(' ',30) 
          

     plot, [0], [0], color=lclr, background=clr.white, charsize=csz,$
           xmargin=xmrg, ymargin=ymrg, xtitle=xtit, $
           ytitle=ytit, yrange=yrng, thick=4, xtickn=xspaces, $
           xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata,  /xlog

     EWMIN = 0.1
     BSZ = 0.1         

     ;; Label
     ypos = 0.8
     dumc = strtrim(round(VINT), 2)
     ;xyouts, 0.12, yrng[1]*ypos, '!9d!Xv = [-'+dumc+','+dumc+'] km s!u-1!N', $
     ;        color=clr.black, charsiz=lsz, align=0.

     ypos = 0.7
     xyouts, 0.12, yrng[1]*ypos, 'R!dperp!N = ['+strtrim(round(Rcuts[0,tt]),2)+','+$
             strtrim(round(Rcuts[1,tt]),2)+'] (kpc)', color=lclr, charsiz=lsz, align=0.

     ;; Cut on R
     idxR = where( all_Rphys GE Rcuts[0,tt] and all_Rphys LT Rcuts[1,tt])
     EWLya = all_EWlya[idxR]

     ;; Histogram
     plothist, alog10(EWLya>EWMIN), all_xhist, all_yhist, bin=BSZ, /fill, fcolor=lclr, /overpl, /noplo
     nbin = n_elements(all_xhist)

     ypt = fltarr(2*nbin)
     ypt[lindgen(nbin)*2] = all_yhist
     ypt[lindgen(nbin)*2+1] = all_yhist
     
     xpt = fltarr(2*nbin)
     xpt[lindgen(nbin)*2] = 10.^(all_xhist-BSZ/2)
     xpt[lindgen(nbin)*2+1] = 10.^(all_xhist+BSZ/2)
     
     ;; Plot
     oplot, xpt, ypt, color=lclr
     flg_fill = 1
     if flg_fill then x_curvefill, xpt, ypt, fltarr(2*nbin), color=lclr

     ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     ;; D_A
     ;; Using Kirkman et al. 2005 (best for our redshifts)
     ;; But augmenting it to include LLS, DLA and metals
     ;; K05_DA = A*(1+z)^gamma with 
     A = 0.0062 
     gamma = 2.75
     ;; DA = DA8s + DA6s + DMs = DA8s + 0.01 + 0.023  [Tytler et al. 2004]
     all_DA = A * (1+gd_strct.z_fg)^2.75  + 0.03
     dlam = (2*VINT)/3e5 * 1215.6701  ;; Ang
     avg_DA = mean(all_DA)
     print, 'Average D_A = ', avg_DA, mean(gd_strct.z_fg)
     ;; Mark
     plotsym, 1, 3.0, thick=4
     ;oplot, dlam*[avg_DA], [yrng[1]*0.93], color=pclr, psym=8
     oplot, replicate(dlam*avg_DA,2), yrng, color=clr.green, linesty=2, thick=2
     ;; Stats
     high = where(gd_strct.EWLya GT (all_DA * dlam), nhigh)
     print, 'Excess = ', nhigh, ' of ', n_elements(all_DA), $
            '     = '+strtrim(float(nhigh)/n_elements(all_DA),2)+'%'

     ;; Empirical Random
     ;ewi = round( median( where(abs(gd_strct.vEW - VINT[tt]) LT 1.)  $
     ;                     MOD n_elements(gd_strct[0].vEW) ) )
     gd_ran = where(abs(gd_strct[idxR].bg_vEWLya[ewi,*]) GT 1e-3, nran)
     
     plothist, alog10((gd_strct[idxR].bg_vEWLya[ewi,*])[gd_ran] > EWMIN), all_xhist, all_yhist, $
               bin=BSZ, /fill, fcolor=clr.blue, /overpl, /noplo
     nbin = n_elements(all_xhist)
     ypt = fltarr(2*nbin)
     ypt[lindgen(nbin)*2] = all_yhist*(float(npair)/nran)
     ypt[lindgen(nbin)*2+1] = all_yhist*(float(npair)/nran)
     
     xpt = fltarr(2*nbin)
     xpt[lindgen(nbin)*2] = 10.^(all_xhist-BSZ/2)
     xpt[lindgen(nbin)*2+1] = 10.^(all_xhist+BSZ/2)
     oplot, xpt, ypt, color=clr.gray;, linest=1

     ;; ;;;;;;;;;;;;;;
     ;; P_KS
     kstwo,  (gd_strct[idxR].bg_vEWLya[ewi,*])[gd_ran], EWLya, d, PKS
     print, 'PKS = ', PKS
     if PKS LT 1e-5 then plbl = 'P!dKS!N < 0.0001'  $
     else  plbl = 'P!dKS!N = 10!u'+string(alog10(PKS), format='(f5.1)')+'!N'

     xyouts, 10.0, yrng[1]*ypos, plbl, color=lclr, charsiz=lsz, align=0.5

  endfor
  ;; Close Ps
  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]

  print, 'fig_hist_ewlya:  All done!'
       
  return
end
      
      
