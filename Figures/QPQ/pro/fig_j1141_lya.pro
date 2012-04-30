;+ 
; NAME:
; fig_lya   
;   Version 1.1
;
; PURPOSE:
;    Plots any array interactively
;
; CALLING SEQUENCE:
;   
;   fig_lya, ydat, [head], XSIZE=, YSIZE=, TITLE=, WAVE=
;
; INPUTS:
;   ydat       - Values 
;   [head]     - Header
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   xsize      - Draw window xsize (pixels)
;   ysize      - Draw window ysize (pixels)
;   wave       - wavelength array
;   ERR        - Error array (fits or image)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_fitdla, 'spec.fits'
;
;
; PROCEDURES/FUNCTIONS CALLED:
;  XGETX_PLT
;  XGETY_PLT
;  XGETXPIX_PLT
;  XGETYPIX_PLT
;
; REVISION HISTORY:
;   17-Oct-2002 Written by JXP
;-
;------------------------------------------------------------------------------


pro fig_lya_UpdatePlot, state, NOWAV=nowav, LTHK=lthk
  
  clr = getcolor(/load)
  if not keyword_set( LTHK ) then lthk = 3

  if keyword_set( NOWAV ) then begin
      ymrg = [0,0] 
      xtn = replicate(' ',30)
  endif else begin
      ymrg = [3.5,0.5]
  endelse

  plot, state.wave[state.xpmnx[0]:state.xpmnx[1]], $
    state.fx[state.xpmnx[0]:state.xpmnx[1]]*state.fscale, psym=state.psym, $
    xrange=[state.xymnx[0],state.xymnx[2]], $
    yrange=[state.xymnx[1],state.xymnx[3]], xstyle=1, ystyle=1, $
    title=state.title, charsize=2.2, xtickinterval=50., $
    background=clr.white, xtickn=xtn, $
    color=clr.black, thick=4, xtitle='Observed Wavelength ('+string("305B)+')', $
    /nodata, ytitle='Relative Flux', $
;    ytitle='Normalized Flux', $
    xmargin=[7,1], ymargin=ymrg

; Crude ERROR
  if state.crude NE 0 and state.nlin NE 0 then begin
      ;; LOW
;      tmp = state.lines[0:state.nlin-1]
;      dla = where(tmp.N GT 20.0 and round(tmp.wrest) EQ 1216L, ndla)
;      tmp[dla].N = tmp[dla].N + state.crude_val[1]
      state.lines.N = state.lines.N + state.crude_val[0]
;      crude_hi = x_allvoigt(state.wave, tmp, SIGMA=state.FWHM)
      x_fitline_updfit, state, /exact, wvoff=200.
      y1 = state.conti[state.xpmnx[0]:state.xpmnx[1]]* $
        state.fit[state.xpmnx[0]:state.xpmnx[1]]*state.fscale 
      state.lines.N = state.lines.N - state.crude_val[0]
;      oplot, state.wave[state.xpmnx[0]:state.xpmnx[1]], $
;        state.conti[state.xpmnx[0]:state.xpmnx[1]]* $
;        state.fit[state.xpmnx[0]:state.xpmnx[1]]*state.fscale, color=clr.red,$
;        thick=lthk
      ;; HIGH
      state.lines.N = state.lines.N + state.crude_val[1]
;      crude_hi = x_allvoigt(state.wave, tmp, SIGMA=state.FWHM)
      x_fitline_updfit, state, /exact, wvoff=200.
      y2 = state.conti[state.xpmnx[0]:state.xpmnx[1]]* $
        state.fit[state.xpmnx[0]:state.xpmnx[1]]*state.fscale 
      state.lines.N = state.lines.N - state.crude_val[1]

      ;; Plot
      x_curvefill, state.wave[state.xpmnx[0]:state.xpmnx[1]], $
                   y1, y2, color=clr.tan
  endif

  ;; Data
  oplot, state.wave[state.xpmnx[0]:state.xpmnx[1]], $
    state.fx[state.xpmnx[0]:state.xpmnx[1]]*state.fscale, psym=state.psym, $
    thick=4, color=clr.black

  ; Plot Error array
  if state.flg_sig EQ 1 then $
    oplot, state.wave[state.xpmnx[0]:state.xpmnx[1]], $
    state.sig[state.xpmnx[0]:state.xpmnx[1]]*state.fscale, $
    psym=state.psym, color=clr.orange, thick=lthk, linestyle=1

  ; Zero line
    oplot, [0., 50000.], [0., 0.], color=clr.gray, linestyle=3, thick=5


  ;; FIT
  if state.nlin NE 0 then begin
      state.lines.N = state.lines.N 
      x_fitline_updfit, state, /exact, wvoff=200.
      oplot, state.wave[state.xpmnx[0]:state.xpmnx[1]], $
        state.conti[state.xpmnx[0]:state.xpmnx[1]]* $
        state.fit[state.xpmnx[0]:state.xpmnx[1]]*state.fscale, color=clr.blue,$
        thick=3

      ;; Mark all DLA
      for i=0L,state.nlin-1 do begin
          if state.lines[i].N GT 20.0 and round(state.lines[i].wrest) EQ 1216L then $
            oplot, replicate( (state.lines[i].zabs+1.)*$
                              state.lines[i].wrest, 2), $
            [state.xymnx[1], state.xymnx[3]], color=clr.green, linestyle=2, $
            thick=3 $
          else $
            oplot, replicate( (state.lines[i].zabs+1.)*$
                              state.lines[i].wrest, 2), $
            [state.xymnx[1], state.xymnx[3]], color=clr.gray, linestyle=4, $
            thick=3
      endfor
  endif

; CONTINUUM
  ;; Line
  oplot, state.wave[state.xpmnx[0]:state.xpmnx[1]], $
    state.conti[state.xpmnx[0]:state.xpmnx[1]]*state.fscale, color=clr.red, $
    linestyle=1, thick=4
  ;; Points
;  if state.cstr.npts NE 0 then begin
;      gdc = where(state.cstr.msk EQ 1)
;      oplot, [state.cstr.xval[gdc]], [state.cstr.yval[gdc]*state.fscale], psym=1, $
;        color=clr.orange, symsize=2
;  endif


  ;; Label
  if tag_exist(state,'csize') then cs = state.csize else cs = 0.7
  if strlen(strtrim(state.label,2)) NE 0 then begin
      xyouts, state.xymnx[0]+(state.xymnx[2]-state.xymnx[0])*0.03, $
        state.xymnx[3]*0.88, $
        state.label, color=clr.black, charsize=cs
  endif

end


;;;;;;;;;;;;;;;;;;;;
;  Reset
;;;;;;;;;;;;;;;;;;;;

pro fig_lya_Reset, state


; Plotting
  state.xymnx = state.svxymnx
  state.old_xymnx = state.svxymnx

end

;;;;;;;;;;;;;;;;;;;;
;  INIT FIT
;;;;;;;;;;;;;;;;;;;;

pro fig_lya_inifit, state, inifit


  ;; File?
  a = findfile(inifit, count=na)
  if na EQ 0 then begin
      print, 'x_fitdla: FILE ', inifit, ' does not exist!'
      stop
  endif

  ;; Restore
  restore, inifit

  ;; Continuum
  if keyword_set(conti) then state.conti[*] = conti $
  else state.conti[*] = 1.

  ;; CSTR
  if size(cstr,/type) EQ 8 then state.cstr = temporary(cstr)

  ;; Lines
  state.nlin = n_elements(lines)
;  state.lines[0:state.nlin-1] = lines
  tmp = state.lines[0:state.nlin-1]
  copy_struct, lines[0:state.nlin-1], tmp
  state.lines[0:state.nlin-1] = tmp
  state.curlin = 0L

  ;; nset
  state.nset = state.lines[state.nlin-1].set

  ;; 
  x_fitline_updfit, state, /exact, wvoff=200.

end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro fig_j1141_lya, yin, ysin, FWHM=fwhm, INIFIT=inifit, INFLG=inflg, XYMNX=xymnx, $
             PSFILE=psfile, FSCALE=fscale, CRUDE=crude


;
  ;if  N_params() LT 1  then begin 
  ;  print,'Syntax - ' + $
  ;           'x_fitdla, fx, ys_in INIFIT=, INFLG=) [v1.1]'
  ;  return
  ;endif 

;  Optional Keywords

  if not keyword_set( XSIZE ) then    xsize = 1200
  if not keyword_set( YSIZE ) then    ysize = 800
  if not keyword_set( FWHM ) then    fwhm = 2.
  if not keyword_set( FSCALE ) then    fscale = 1.

; Read in the Data
  label = 'J1141+0724: z=3.544,  N!dHI!N = 10!u19.9!N cm!u-2!N'
  psfile = 'fig_j1141_lya.ps'
  bg_qso = '~/Dropbox/QSOPairs/data/lya_spec_all/LRIS_GMOS_PRELIM/J1141+0724A.fits'
  ydat = x_readspec(bg_qso, INFLG=2, head=head, NPIX=npix, WAV=xdat,  SIG=ysig)

  if not keyword_set(YSIG) then ysig = replicate(1., npix)

  tmp1 = { newabslinstrct }

  tmp2 = { conti_str, $
           npts: 0L, $
           xval: fltarr(100), $
           yval: fltarr(100), $
           msk: lonarr(100) }


;    STATE
  state = { fx: ydat, $
            wave: xdat, $
            sig: ysig, $
            npix: npix, $
            flg_sig: 0, $
            flg_zoom: 0, $
            pos: [0.1,0.1,0.95,0.95], $ ; Plotting
            FWHM: fwhm, $   ; FWHM of instrument (pix)
            nlin: 0, $ ; DLA flag
            lines: replicate(tmp1,20), $ ; DLA flag
            fit: fltarr(n_elements(ydat)) + 1., $
            conti: replicate(1.,npix), $
            label: label, $
            csize: 2., $
            cstr: tmp2, $
            curlin: 0L, $
            nset: -1, $
            crude: 1, $  ; Crude Error
            crude_val: [-0.1, 0.1], $
            xpos: 0.d, $
            ypos: 0.d, $
            flg_dum: 1, $
            psfile: 0, $ ; Postscript
            svxymnx: [0., $     ; xmin
                      min(ydat)-0.01*abs(max(ydat)-min(ydat)), $ ; ymin
                      float(n_elements(ydat)-1), $ ; xmax
                      max(ydat)+0.01*abs(max(ydat)-min(ydat))], $ ; ymax
            xymnx: fltarr(4), $
            old_xymnx:fltarr(4), $
            tmpxy: fltarr(4), $
            xpmnx: lonarr(2), $
            psym: 10, $
            title: '', $
            help: strarr(30), $
            size: lonarr(2), $
            xcurs: 0.0, $
            ycurs: 0.0, $
            fscale: fscale $
          }

; WAVE
  if keyword_set(WAVE) then state.wave = wave
  if keyword_set( YSIG ) then state.flg_sig = 1

  resolve_routine, 'x_specplot', /NO_RECOMPILE
  resolve_routine, 'x_fitline', /NO_RECOMPILE

; Update
  fig_lya_Reset, state

  ;; Init Fit
  inifit = 'j1141_lya.idl'
  if keyword_set(INIFIT) then fig_lya_inifit, state, inifit

; CRUDE
  if keyword_set(CRUDE) then begin
      state.crude = 1
      state.crude_val = crude
  endif

; XYMNX
  xymnx = [5480., -25., 5576, 224]
  if not keyword_set(XYMNX) then begin
      state.svxymnx[0] = min(state.wave)
      state.svxymnx[2] = max(state.wave)
  endif else begin
      state.svxymnx[0] = xymnx[0]
      state.svxymnx[2] = xymnx[2]
      state.svxymnx[1] = xymnx[1]
      state.svxymnx[3] = xymnx[3]
      state.xymnx = state.svxymnx
  endelse

  ;; Set pmnx
  state.xpmnx = x_getxpmnx(state)

; PSFILE
  x_psopen, psfile, /maxs
  fig_lya_UpdatePlot, state  ;; PLOT
  x_psclose

  return
end
