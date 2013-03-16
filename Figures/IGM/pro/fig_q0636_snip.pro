;; Snippets of the Q0636 spectrum then the whol thing
pro fig_q0636_snip

  if not keyword_set(CSZ) then csz = 1.7
  if not keyword_set(lSZ) then lsz = 2.0
  if not keyword_set(lSZ2) then lsz2 = 1.5

  ;; Read
  q0636 = x_readspec('~/Keck/HIRES/RedData/Q0636+680/Q0636+680.fits',infl=2, /struct)

  ;; Sky
  sky = where( q0636.wv GT 5578.5 and q0636.wv LT 5579.8)
  q0636.fx[sky] = 395.
  sky = where( q0636.wv GT 5890.6 and q0636.wv LT 5892.4)
  q0636.fx[sky] = 459
  sky = where( q0636.wv GT 5896.6 and q0636.wv LT 5898.6)
  q0636.fx[sky] = 465
  sky = where( q0636.wv GT 6301.5 and q0636.wv LT 6303.2)
  q0636.fx[sky] = 466

  ;; Emission Plot
  x_psopen, 'fig_q0636_snip.ps', /maxs
  !p.multi = [0,1,1]
  clr = getcolor(/load)

  ;; Eta formulation from Steidgman 2007
  xmrg = [9,1]
  ymrg = [3,0.1]
  pos = [0.15, 0.58, 0.9, 0.9]

  xspaces = replicate(' ',30) 


  axrng=[[4120, 4170], $
         [4500, 4580], $
         [4250, 4320], $
         [3150, 7000]]
  ayrng = [[-5., 100], $
           [-5, 200], $
           [-5, 150], $
           [-5, 650]]
  sz = size(axrng,/dimen)

  ytit = 'Relative Flux'
  xtit = 'Wavelength (Ang)'

  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Q1946
  for qq=0L,sz[1]-1 do begin
     xrng = axrng[*,qq]
     yrng = ayrng[*,qq]
     if qq EQ (sz[1]-1) then begin
        delvarx, xspaces
        thick=2
     endif else begin
        thick = 4
     endelse
     plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
           xmargin=xmrg, ymargin=ymrg, ytitle=ytit, $
           xtitle=xtit, yrange=yrng, thick=5, xtickn=xspaces, $
           pos = pos, $
           xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata
     
     ;;
     oplot, xrng, [0., 0], color=clr.gray, linesty=2, thick=3

     if qq EQ (sz[1]-1) then begin
        for ii=0L,sz[1]-2 do begin
           xrng = axrng[*,ii]
           yrng = ayrng[*,ii]
           x_curvefill, xrng, replicate(yrng[1],2), [0., 0], color=clr.cyan
        endfor
     endif
     ;; Data
     oplot, q0636.wv, q0636.fx, psym=10, color=clr.black, thick=thick
     
  endfor

     
  x_psclose
  !p.multi = [0,1,1]

  return

end
