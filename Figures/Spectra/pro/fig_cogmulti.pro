pro fig_cogmulti

  common x_calccog_cmm, cog_bval, cog_tau, cog_strct

  if not keyword_set(PSFILE) then psfile = 'fig_cogmulti.ps'
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
  plot, [0.], [0.], color=clr.lightgray, background=clr.black, charsize=1.8,$
    xmargin=[9.5,2], ymargin=[5,5], xtitle='N', $
    ytitle='W!d!9l!X!N (Ang)', /ylog, $
    /nodata, xrange=[1e11, 1e22], ystyle=1, yrange=[0.001, 100.],  xstyle=1, /xlog

  ;; v=10km/s

  cogmax_fil = getenv('XIDL_DIR')+'/Spec/Analysis/cogmax_tab.fits'
  cog_strct = xmrdfits(cogmax_fil, 1, /silent)
  compile_opt strictarr
  resolve_routine, 'x_calccog', /COMPILE_FULL_FILE, /EITHER
  cog_bval = 10.
  
  for ss=0,3 do begin
      ;; Init 
      case ss of
          0: begin
             wrest = 1215.6701d
             cc = clr.cyan
             gamma = 6.265d8
          end
          1: begin
             wrest = 1808.0130
             cc = clr.yellow
             gamma = 6.749e6
          end
          2: begin
             wrest = 1250.584
             cc = clr.green
             gamma = 5.260e7
          end
          3: begin
             wrest = 1302.1685
             cc = clr.orange
             gamma = 5.750e8
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
          
      oplot, 10^cog_nval, EW, color=cc
      ;; Label
      xyouts, 10^(15.), 0.7*10^(0.35*ss), nam, $
              color=cc, charsiz=lsz
  endfor

  if keyword_set( PSFILE ) then x_psclose

  return
end
