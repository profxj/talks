;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  Plots sightline distribution
pro qpq5_survey, COVFACT = COVFACT, psfile = psfile, outfile = outfile, NOPS=nops, $
                  CSZ=csz, XMRG=xmrg

  compile_opt strictarr

  ;; Get structure if necessary
  if not keyword_set( PSFILE ) then psfile = 'qpq5_survey.ps'
  if not keyword_set( CSZ ) then csz = 1.9
  if not keyword_set( LSZ ) then lsz = 1.9

  if not keyword_set(SIG_FX) then sig_fx = 0.05 ;; Error in tau values
  if not keyword_set( MAXCHI ) then MAXCHI = 1.3
  ;  if not keyword_set( WVMIN ) then wvmin = 700. 
  if not keyword_set(OMEGA_M) then omega_m = 0.3

  if not keyword_set( NRNG ) then nrng = 50L

  ;; Compile
  cd, getenv('QSO_DIR')+'tex/QPQ5/Analysis/pro/', curr=curr
  resolve_routine, 'qpq5_get_instrument', /compile_full_file, /is_function
  cd, curr
  qso_dir =  getenv('QSO_DIR')

  ;; Read in pair stuff
  qpq_fil = '~/Dropbox/QSOPairs/qpq5_pairs.fits'
  qpq5 = xmrdfits(qpq_fil, 1)
  npair = n_elements(qpq5)

  ;; Get angles
  posang, 1, qpq5.rad/15., qpq5.decd, qpq5.rad_bg/15., qpq5.decd_bg, pa

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; PLOT
  x_psopen, psfile, /maxs
  !p.multi=[0,1,1]
  xmrg = [7,17]
  ymrg = [4.0,0.5]
  clr = getcolor(/load)
  lclr = clr.black
  ;lclr = clr.white

  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Covering fraction

  xrng = [-303, 303]
  yrng = xrng

  plot, [0], [0], color=lclr, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, xtitle=xtit, $
        ytitle=ytit, yrange=yrng, thick=4, xtickn=xspaces, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata


  ;; Loop on Radii
  for qq=0L,npair-1 do begin
     
     instr = qpq5_get_instrument( strtrim(qso_dir+qpq5[qq].lya_fil,2) )
     case instr of
        'SDSS': clbl = 'S'
        'BOSS': clbl = 'B'
        'GMOS': clbl = 'G'
        'LRIS': clbl = 'L'
        'LRISb': clbl = 'L'
        'MagE': clbl = 'M'
        'MIKE': clbl = 'M'
        'ESI': clbl = 'E'
        else: stop
     endcase

     if qpq5[qq].z_fg LT 2.5 then cclr = clr.blue else cclr = clr.red

     xpos = qpq5[qq].R_phys * cos(pa[qq]*!pi/180.)
     ypos = qpq5[qq].R_phys * sin(pa[qq]*!pi/180.)

     xyouts, xpos, ypos, clbl, color=cclr, charsi=lsz, align=0.5

  endfor

  ;; Background
  radius = 300.
  points = (2 * !PI / 1000.0) * FINDGEN(1000)
  x = radius * COS(points )
  y = radius * SIN(points )
  x = [x,x[0],x[1]]
  y = [y,y[0],y[1]]
  oplot, x,y,color=clr.gray, thick=4

  ;; Close Ps
  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]

  print, 'fig_covering:  All done!'
  return
end
      
      
