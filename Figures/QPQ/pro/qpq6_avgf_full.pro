;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  Avg F figure for the full QPQ6 sample
;;;    Includes cumulative comparison
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro qpq6_avgf_full

  compile_opt strictarr

  c=x_constants()
  ;cd, getenv('QSO_DIR')+'tex/QPQ6/Analysis/AbsLines/pro/', curr=curr
  ;resolve_routine, 'qpq6_mkdata', /compile_full_file, /is_function
  ;resolve_routine, 'qpq6_grab_idlines', /compile_full_file, /is_function
  ;cd, '../Analysis/pro/'
  ;resolve_routine, 'qpq6_anly_cut', /compile_full_file, /is_function
  ;resolve_routine, 'qpq6_get_instrument', /compile_full_file, /is_function
  ;cd, curr

  ;; Get structure if necessary
  if not keyword_set( PSFILE ) then psfile = 'qpq6_avgf_full.ps'
  if not keyword_set( CSZ ) then csz = 1.8
  if not keyword_set( LSZ ) then lsz = 1.5

  if not keyword_set(SIG_FX) then sig_fx = 0.05 ;; Error in tau values
  if not keyword_set( MAXCHI ) then MAXCHI = 1.3
  if not keyword_set( NRNG ) then nrng = 50L
  if not keyword_set( SHIFT ) then shift = 0.
  if not keyword_set( SCALE ) then scale = 1.
  ;if n_elements(XLOG) EQ 0 then xlog = 1

  qsodir = getenv('QSO_DIR')

  ;; Read in pair structure
  qpq_fil = '~/Dropbox/QSOPairs/qpq6_final_cut.fits'
  qpq_strct = xmrdfits(qpq_fil, 1)

  ;gdp = qpq6_anly_cut(qpq_strct, NPAIRS=npair)
  ;gd_strct = qpq_strct[gdp]
  gd_strct = qpq_strct
  npair = n_elements(gd_strct)

  ;; Velocity Intervals
  aVINT = [500., 1000., 1600.]
  nVINT = n_elements(aVINT)


  npx = 1
  npy = 1
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; PLOT
  x_psopen, psfile, /maxs
  !p.multi=[0,1,1]
  clr = getcolor(/load)
  fclr = clr.white
  pclr = clr.cyan

  ;; Histogram first
  xmrg = [18,21]
  ymrg = [14,4.]

  xrng = [0., 1]
  yrng = [0., 150]

  ytit='N'
  xtit='<F>!d1000!N'
          
  plot, [0], [0], color=fclr, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, xtitle=xtit, $
        ytitle=ytit, yrange=yrng, thick=4, xtickn=xspaces, $
        xrange=xrng, ystyle=1, xstyle=9, psym=1, /nodata,  xlog=xlog

  FMIN = 0.0
  FMAX = 0.999
  BSZ = 0.05         

  ;; Label
  ;ypos = 0.8
  ;dumc = strtrim(round(VINT[tt]), 2)

  ypos = 0.9
  xyouts, 0.1, yrng[1]*ypos, 'R!dphys!N < 1Mpc', color=fclr, charsiz=lsz, align=0.0

  ypos = 0.75
  xyouts, 0.1, yrng[1]*ypos, 'CGM', color=clr.cyan, charsiz=lsz, align=0.0
  ypos = 0.65
  xyouts, 0.1, yrng[1]*ypos, 'Control', color=clr.gray, charsiz=lsz, align=0.0

  ;; ;;;;;;;;;;;;;;;;;;;;;;;
  ;; Data
  VINT = 1000.
  dLam = 1215.6701 * (2*VINT) / 3e5 ; Ang
  ewi = round( median( where(abs(gd_strct.vEW - VINT) LT 1.)  $
                       MOD n_elements(gd_strct[0].vEW) ) )
  EWLya = gd_strct.vEWLya[ewi]
  avgF = 1. - EWLya/dLam
  ;stop

  ;; Histogram
  plothist, (FMAX < avgF>FMIN), all_xhist, all_yhist, bin=BSZ, /fill, fcolor=pclr, /overpl, noplo=XLOG, $
            color=pclr
  
  ;; Plot
  ;oplot, xpt, ypt, color=clr.black
  ;flg_fill = 1
  ;if flg_fill then x_curvefill, xpt, ypt, fltarr(2*nbin), color=clr.black

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
  oplot, [1-avg_DA], [yrng[1]*0.93], color=clr.gray, psym=8

  ;print, tt, 'DA = ', dlam*avg_DA
  ;; Stats
  high = where(gd_strct.EWLya GT (all_DA * dlam), nhigh)
  print, 'Excess = ', nhigh, ' of ', n_elements(all_DA), $
         '     = '+strtrim(float(nhigh)/n_elements(all_DA),2)+'%'

  ;; Empirical Random
  gd_ran = where(abs(gd_strct.bg_vEWLya[ewi,*]) GT 1e-3, nran)
  rEWLya = (gd_strct.bg_vEWLya[ewi,*])[gd_ran] 
  ravgF = 1. - rEWLya/dLam
  
  plothist, (FMIN > ravgF < FMAX), all_xhist, all_yhist, $
            bin=BSZ, /fill, fcolor=clr.blue, /overpl, /noplo
  nbin = n_elements(all_xhist)
  ypt = fltarr(2*nbin)
  ypt[lindgen(nbin)*2] = all_yhist*(float(npair)/nran)
  ypt[lindgen(nbin)*2+1] = all_yhist*(float(npair)/nran)
  
  xpt = fltarr(2*nbin)
  xpt[lindgen(nbin)*2] = (all_xhist-BSZ/2)
  xpt[lindgen(nbin)*2+1] = (all_xhist+BSZ/2)
  oplot, (xpt+SHIFT)*SCALE, ypt, color=clr.gray, linest=((SCALE GT 1.) + (SHIFT GT 0.))
  
  
  ;; EW axis
  xrng = dLam * (1.-[0., 1])
  axis, xaxis = 1, color = fclr, charsi = csz, xrang = xrng, $
        xsty = 1, xtit = 'W!d1000!N (Ang)', /save


  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; CUMULATIVE PLOT
  xrng = [0., 1]
  yrng = [0., 1.]

  ytit='Cumulative'
          
  plot, [0], [0], color=fclr, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, xtitle=xtit, $
        ytitle=ytit, yrange=yrng, thick=4, xtickn=xspaces, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata,  xlog=xlog

  ;; Label
  ypos = 0.75
  xyouts, 0.1, yrng[1]*ypos, 'CGM', color=clr.cyan, charsiz=lsz, align=0.0
  ypos = 0.65
  xyouts, 0.1, yrng[1]*ypos, 'Control', color=clr.gray, charsiz=lsz, align=0.0

  ;; Data
  x_cumulative, avgF, xval, yval
  oplot, xval, yval, color=pclr

  ;; Other velocity intervals
  if keyword_set(SHOW_ALL) then begin
     for tt=0L,nVINT-1 do begin
        VINT = aVINT[tt]
        dLam = 1215.6701 * (2*VINT) / 3e5 ; Ang
        ewi = round( median( where(abs(gd_strct.vEW - VINT) LT 1.)  $
                             MOD n_elements(gd_strct[0].vEW) ) )
        EWLya = gd_strct.vEWLya[ewi]
        avgF = 1. - EWLya/dLam
        
        x_cumulative, avgF, xval, yval
        oplot, xval, yval, color=fclr, linesty=tt+1
     endfor
  endif

  ;; Random
  x_cumulative, ravgF, xval, yval
  oplot, xval, yval, color=clr.gray

  ;; ;;;;;;;;;;;;;;
  ;; P_KS
  kstwo,  avgF, ravgF, d, PKS
  print, 'PKS = ', PKS
  if PKS LT 1e-5 then plbl = 'P!dKS!N < 0.0001'  $
  else  plbl = 'P!dKS!N = '+string(PKS, format='(f6.4)')
;    else  plbl = 'P!dKS!N = 10!u'+string(alog10(PKS), format='(f5.1)')+'!N'
  ypos = 0.9
  xyouts, 0.1, yrng[1]*ypos, plbl, color=fclr, charsiz=lsz, align=0.0
  ;                              ;stop

  ;; Close Ps
  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]

  print, 'fig_hist_ewlya:  All done!'
       
  return
end
      
      
