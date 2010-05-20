pro fig_cog

  common x_calccog_cmm, cog_bval, cog_tau, cog_strct

  if not keyword_set(PSFILE) then psfile = 'fig_cog.ps'
  c = x_constants()
  Nnval = 50L
  ew = fltarr(Nnval)
  cog_Nval = 12. + findgen(Nnval)*10./(NNval-1)
  nplt = 10000L
  lsz = 2.2

  ;; Calculate

  ;; Plot
  if keyword_set( PSFILE ) then x_psopen, psfile, /portrait
  clr = getcolor(/load)
  plot, [0.], [0.], color=clr.black, background=clr.white, charsize=1.8,$
    xmargin=[9.5,2], ymargin=[4,1], xtitle='N!dHI!N', $
    ytitle='Rest EW (Ang)', /ylog, $
    /nodata, xrange=[1e12, 1e22], ystyle=1, yrange=[0.01, 100.],  xstyle=1, /xlog

  ;; v=10km/s

  cogmax_fil = getenv('XIDL_DIR')+'/Spec/Analysis/cogmax_tab.fits'
  cog_strct = xmrdfits(cogmax_fil, 1, /silent)
  compile_opt strictarr
  resolve_routine, 'x_calccog', /COMPILE_FULL_FILE, /EITHER
  
  wrest = 1215.6701d
  getfnam, wrest, fval, nam
  flambda = wrest*fval

  for ss=0,2 do begin
      ;; Init 
      case ss of
          0: begin
              cc = clr.blue
              cog_bval = 10.
              xrng = 1215.67+[-100., 100]
          end
          1: begin
              cc = clr.red
              cog_bval = 30.
              xrng = 1215.67+[-100., 100.]
          end
          2: begin
              cc = clr.darkgreen
              cog_bval = 100.
              xrng = 1215.67 + [-100,100]
          end
          else: stop
      endcase

      wplt = xrng[0] + findgen(nplt)*(xrng[1]-xrng[0])/(nplt-1)
      nu0 = c.c / (1215.67 * 1e-8)
      dnu = (wplt-1215.67) / 1215.67 * nu0 
      gamma = 6.265d8

      ;; Brute force
      for jj=0L,Nnval-1 do begin
          tau0 = !pi * c.e^2 * 0.4164 / c.me / c.c * 10^cog_nval[jj]
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
      xyouts, 10^(13.1), 10*10^(0.3*ss), 'b = '+ $
              string(round(cog_bval),format='(i3)')+'km/s', $
              color=cc, charsiz=lsz
  endfor

  if keyword_set( PSFILE ) then x_psclose

  return
end
