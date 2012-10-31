;; qpq6_lyaspec_image, IMG=img, DAT_STRCT=dat_strct, vmnx=vmnx
pro qpq6_lyaspec_image, wrest, dv, grid_file, psfile, YMNX=ymnx, XRNG=xrng, YRNG=yrng, $
                   IMG=img, DAT_STRCT=dat_strct, VMNX=vmnx

  if not keyword_set(PAD_FRAC) then pad_frac = 0.1
  if not keyword_set(CSZ) then csz = 1.1
  if not keyword_set(CSZ3) then csz3 = 0.9
  if not keyword_set(lSZ) then lsz = 1.5
  if not keyword_set( PSFILE ) then psfile = 'qpq6_lyaspec_image.ps'

  ;; Grab the stack image
  compile_opt strictarr

  Rcuts = [ [0., 100], $  ;; First one is necessary to get things to run
            [0., 1000] ]

  c=x_constants()
  cd, getenv('QSO_DIR')+'tex/QPQ6/Analysis/AbsLines/pro/', curr=curr
  resolve_routine, 'qpq6_mkdata', /compile_full_file, /is_function
  resolve_routine, 'qpq6_grab_idlines', /compile_full_file, /is_function
  cd, getenv('QSO_DIR')+'tex/QPQ6/Analysis/pro/'
  resolve_routine, 'qpq6_anly_cut', /compile_full_file, /is_function
  resolve_routine, 'qpq6_get_instrument', /compile_full_file, /is_function
  cd, getenv('QSO_DIR')+'tex/QPQ6/Analysis/Stacks/pro/'
  resolve_routine, 'qpq6_stack_lya', /compile_full_file
  cd, curr

  if not keyword_set(IMG) then begin
     qpq6_stack_lya, /nocut, S2N=8., /NEW, STACK_IMG=img, /NOWRITE, RCUTS=rcuts, $
                     VMNX=vmnx, DAT_STRCT=dat_strct
     img = transpose(img)
  endif

  if not keyword_set(XRNG) then xrng=[0., n_elements(dat_strct)-1]
  if not keyword_set(YRNG) then yrng=vmnx

  ;; IFU shots next
  x_psopen, psfile, /maxs
  clr = getcolor(/load)
  !p.multi=[0,1,1]
  fclr = clr.white

  devicefactor=2540.  ;; cm to inch
  x0 = 0.8
  y0 = 0.7
;  stretch_lo = 0.01


  yrange = [0.3, 1.0]

  common xcommon_color, r_vector, g_vector, b_vector
  xicmm_colors
  state = { $
          ncolors: 0L, $
          brightness: 0.5, $
          contrast: 0.35 $
          }
  ;loadct, 0, /silent
  ncolors = 160L
  ctload, 0, ncolors=ncolors
  ;ncolors = !d.table_size - 9
  state.ncolors=ncolors
  r_vector = bytarr(ncolors)
  g_vector = bytarr(ncolors)
  b_vector = bytarr(ncolors)
  ximgd_getct, state, 0;, /CLR
  ;ximgd_stretchct, state


  for ss=0,1 do begin
     case ss of
        0: begin ;; "Raw"
           xpos1 = [0.8]
           ypos1 = [5.]
           uimg = img
           yimsize = 1.9 ;; inches
        end
        1: begin ;; "Smoothed
           ypos1 = [1.5]
           ismth = 10L
           sz = size(img, /dimen)
           uimg = (smooth(img,ismth))[ismth/2:sz[0]-ismth/2-1, ismth/2:sz[1]-ismth/2-1]
           yimsize = yimsize * 0.92
        end
     endcase

     ;; 
     scaled = bytscl(uimg, min=yrange[0], max=yrange[1], $
                     /nan, $
                     top=(ncolors - 1))
     
     ;; Plot
     qq = 0
                                ;ctload, 0, ncolors=xncolor, /rever
     tv, scaled, xpos1[qq], ypos1[qq], ysize=yimsize, /inches
  
     ;; Axes
     dims = size(scaled,/dim)
     xlabel = findgen(dims[0]+1)
     ylabel = findgen(dims[1]+1)
     thisPosition = devicefactor*[xpos1[qq], $
                                  ypos1[qq], $
                                  xpos1[qq]+(yimsize*dims[0]/dims[1]), $
                                  ypos1[qq]+yimsize]
     ;ctload, 0, ncolors=xncolor
     plot, xlabel, ylabel, /nodata, /device, /noerase, $
           position=thisPosition, color=fclr, $
           xrange=[min(xlabel),max(xlabel)], $
           yrange=[min(ylabel),max(ylabel)], $
           xstyle=5, ystyle=5

     clr = getcolor(/load)
     yspaces = replicate(' ', 30)
     xspaces = replicate(' ', 30)
     ytit = 'Relative Velocity (km s!u-1!N)'
     xtit = 'R!dphys!N (kpc)'
     plot, [0], [0], /device, /noerase, xrange=xrng, color=fclr, $
           yrange=yrng, xtitle=xtit, charsiz=csz2, $ ;ytickn=yspaces, $
           xtickn=yspaces, $
           ytitle=ytit, /nodata, xsty=13, ysty=1, $
           position=thisPosition  
     
     
     ;; X-labels
     xyouts, mean(xrng), yrng[0]-0.2*(yrng[1]-yrng[0]), xtit, color=fclr, charsi=csz, align=0.0
     xlbl = [round(min(dat_strct.R)), 100, 200, 500, 750, 1000]
     nxl = n_elements(xlbl)
     for ii=0L,nxl-1 do begin
        mn = min(abs(dat_strct.R-xlbl[ii]), imn)
        xyouts, imn, yrng[0]-0.1*(yrng[1]-yrng[0]), strtrim(round(xlbl[ii]),2), $
                color=fclr, charsi=csz, align=0.5
     endfor
  endfor

  ;; Color bar
  coyote_colorbar, pos=[0.48, 0.10, 0.52, 0.52], range=yrange, $
                   ncolor=ncolors, FORMAT='(f4.1)', color=fclr, $
                   title='Normalized Flux',charsiz=csz3 

  ;; Y-labels
     
  ;if round(dv[qq]) LT 0. then pclr = clr.blue else pclr=clr.red
  ;if abs(round(dv[qq])) LT 10. then pclr = clr.black
  ;xyouts, ymnx[0]+0.05*(ymnx[1]-ymnx[0]), ymnx[0]+0.85*(ymnx[1]-ymnx[0]),  $
  ;        'v='+strtrim(round(dv[qq]),2),color=pclr, charsiz=lsz
  
  if keyword_set( PSFILE ) then x_psclose
  !p.multi = [0,1,1]

  return

end
