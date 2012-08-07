pro lbg_cog     

  common x_calccog_cmm, cog_bval, cog_tau, cog_strct

  if not keyword_set(PSFILE) then psfile = 'lbg_cog.ps'
  c = x_constants()
  Nnval = 50L
  ew = fltarr(Nnval)
  cog_Nval = 11. + findgen(Nnval)*11./(NNval-1)
  nplt = 10000L
  lsz = 2.2

  ;; Calculate

  ;; Plot
  if keyword_set( PSFILE ) then x_psopen, psfile, /maxs
  clr = getcolor(/load)
  plot, [0.], [0.], color=clr.black, background=clr.white, charsize=1.8,$
        xmargin=[9.5,2], ymargin=[5,5], xtitle='Column Density (cm!u-2!N)', $
        ytitle='Equivalent Width  (Ang)', /ylog, $
        /nodata, xrange=[1e11, 1e20], ystyle=1, yrange=[0.001, 100.],  xstyle=1, /xlog, $
        ytickformat='x_logticks'

  ;; v=10km/s

  cogmax_fil = getenv('XIDL_DIR')+'/Spec/Analysis/cogmax_tab.fits'
  cog_strct = xmrdfits(cogmax_fil, 1, /silent)
  compile_opt strictarr
  resolve_routine, 'x_calccog', /COMPILE_FULL_FILE, /EITHER
  cog_bval = 30.
  
  for ss=0,2 do begin
      ;; Init 
      case ss of
          0: begin
             wrest = 1215.6701d
             cc = clr.red
             gamma = 6.265d8
             lsty = 0
          end
          1: begin
             wrest = 1334.5323
             cc = clr.blue
             gamma = 2.870e8
             lsty = 0
          end
          2: begin
             wrest = 1215.6701d
             cc = clr.red
             lsty = 2
             gamma = 6.265d8
             cog_bval = 55.
          end
          else: stop
       endcase

      getfnam, wrest, fval, nam
      flambda = wrest*fval
      xrng = wrest+[-100., 100]


      wplt = xrng[0] + findgen(nplt)*(xrng[1]-xrng[0])/(nplt-1)
      nu0 = c.c / (wrest * 1e-8)
      dnu = (wplt-wrest) / wrest * nu0 

      ;; Brute force
      for jj=0L,Nnval-1 do begin
          tau0 = !pi * c.e^2 * fval / c.me / c.c * 10^cog_nval[jj]
          ;; make the Voigt
          nuD = nu0 * cog_bval*1e5 / c.c
          a = gamma / (4*!dpi*nuD)
          u = dnu / nuD
          tau_nu = tau0 * voigt(a,u) / nuD / sqrt(!dpi)
          ;; EW
          fx = exp(-1.*tau_nu)
          dwv = wplt - shift(wplt,1)
          ew[jj] = total((1.-fx)*dwv)
      endfor
          
      oplot, 10^cog_nval, EW, color=cc, linesty=lsty
      ;; Label
      if ss LE 1 then $
         xyouts, 10^(17.), 0.007*10^(0.35*ss), nam, $
                 color=cc, charsiz=lsz
  endfor

  if keyword_set( PSFILE ) then x_psclose

  return
end
