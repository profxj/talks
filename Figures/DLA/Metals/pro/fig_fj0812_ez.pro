;; Shows Best estimates of D/H and thereby Omega b
pro fig_fj0812_ez, NARROW=narrow

  if not keyword_set(CSZ) then csz = 2.1
  if not keyword_set(lSZ) then lsz = 1.6
  if not keyword_set(lSZ2) then lsz2 = 1.5

  ;; Parse
  if not keyword_set(DLA) then parse_dlalst, dla, $
    '/u/xavier/DLA/Lists/mstrong.lst', ROOT=getenv('DLA')
  
  mt = where(strtrim(dla.qso,2) EQ 'J0812+3208' AND $
             abs(dla.zabs-2.6263) LT 1e-3, nmt)
  if nmt NE 1 then stop
  fj0812 = dla[mt]

  ;; Emission Plot
  psfile= 'fig_fj0812_ez.ps'
  if keyword_set(NARROW) then psfile = 'fig_fj0812_ez_narrow.ps'

  x_psopen, psfile, /maxs
  !p.multi = [0,1,1]
  clr = getcolor(/load)

  ;; Eta formulation from Steidgman 2007
  xmrg = [8,17]
  ymrg = [7,3]

  yrng=[0., 10]
  xrng=[4, 38]

  ytit = '!9e!X(X)'
  xtit = 'Atomic Number'

  asz = 1.9

  flgelm = intarr(100)

  deplete = ['Cl', 'Ti', 'Fe', 'Ni', 'Cr', 'Co','Al']
  deplete2 = ['Mn', 'Cu']

  ;; Plot
  pos1=[0.15, 0.2, 0.7, 0.7]
  if keyword_set(NARROW) then begin
     pos1[3] = 0.55
     csz = 1.6
     lsz = 1.4
  endif

  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, ytitle=ytit, $
        xtitle=xtit, yrange=yrng, thick=5, $
        pos=pos1, $
        xrange=xrng, ystyle=9, xstyle=1, psym=1, /nodata
  
  pts = where(fj0812.XH.flgclm GT 0, npts)
  pts = pts[1:*] ;; Remove Hydrogen
  npts = npts-1
  sv_abnd = fltarr(npts)

  for kk=0L,npts-1 do begin
     getabnd, elm, pts[kk], abnd, flag=1
     sv_abnd[kk] = abnd
  endfor
  oplot, pts, sv_abnd+fj0812.XH[8].clm, color=clr.gray, linesty=2, thick=5

  for kk=0L,npts-1 do begin
     if pts[kk] GT xrng[1] then continue
     case fj0812.XH[pts[kk]].flgclm of
        1: begin
           plotsym, 0, 1.5, /fill
           pclr = clr.blue
           yoff = 0.
        end
        2: begin ;; Lower limit
           plotsym, 2, asz, /fill, thick=5
           pclr = clr.darkgreen
           yoff = 0.2
        end
        3: begin
           plotsym, 1, asz, /fill, thick=5
           pclr = clr.red
           yoff = -0.2
        end
        else: stop
     endcase


           
     getabnd, elm, pts[kk], abnd, flag=1

     ;; Depleted?
     doff = 0.
     if total( strmatch(deplete, elm) ) GT 0 then begin
        plotsym, 2, asz, /fill, thick=5
        pclr = clr.darkgreen
        doff = 0.7
        yoff += 0.2
     endif
     if total( strmatch(deplete2, elm) ) GT 0 then doff = 0.7

     ;; Plot
     oplot, [pts[kk]], [fj0812.XH[pts[kk]].clm+abnd+doff], color=pclr, psym=8
     ;; Label
     if (pts[kk] MOD 2) EQ 0 then yoff += 0.4 else yoff +=-0.9

     xyouts, pts[kk], fj0812.XH[pts[kk]].clm+abnd+yoff+doff, elm, $
             color=clr.black, charsi=lsz, align=0.5
  endfor


  oplot, [20., 23], [9., 9], color=clr.gray, linesty=2, thick=5
  xyouts, 24., 8.8, 'Scaled Solar', color=clr.darkgray, charsi=lsz, align=0.
  
  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Plot high Z elements
  pos2=pos1
  pos2[0] = pos1[2] + 0.02
  pos2[2] = pos1[2] + 0.10
  xrng=[48, 90]

  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, $
        yrange=yrng, thick=5, $
        pos = pos2, xtickint=20, $
        xrange=xrng, ystyle=5, xstyle=1, psym=1, /nodata, /noerase
  
  pts = where(fj0812.XH.flgclm GT 0, npts)
  pts = pts[1:*] ;; Remove Hydrogen
  npts = npts-1
  sv_abnd = fltarr(npts)

  for kk=0L,npts-1 do begin
     if pts[kk] LT xrng[0] then continue
     case fj0812.XH[pts[kk]].flgclm of
        1: begin
           plotsym, 0, 1.5, /fill
           pclr = clr.blue
        end
        2: begin ;; Lower limit
           plotsym, 2, asz, /fill, thick=5
           pclr = clr.darkgreen
        end
        3: begin
           plotsym, 1, asz, /fill, thick=5
           pclr = clr.red
        end
        else: stop
     endcase


     if (pts[kk] MOD 2) EQ 0 then yoff = 0.4 else yoff=-0.8
           
     getabnd, elm, pts[kk], abnd, flag=1

     ;; Depleted?
     doff = 0.
     if total( strmatch(deplete, elm) ) GT 0 then begin
        plotsym, 2, asz, /fill, thick=5
        pclr = clr.darkgreen
        doff = 0.7
     endif
     sv_abnd[kk] = abnd
     ;; Plot
     oplot, [pts[kk]], [fj0812.XH[pts[kk]].clm+abnd+doff], color=pclr, psym=8
     ;; Label
     xyouts, pts[kk], fj0812.XH[pts[kk]].clm+abnd+yoff+doff, elm, $
             color=clr.black, charsi=lsz, align=0.5
  endfor

  oplot, pts, sv_abnd+fj0812.XH[8].clm, color=clr.gray, linesty=2, thick=5
  axis, yaxis = 1, color = clr.black,  yrang = yrng, ysty = 1, ytickn=replicate(' ', 30)

  x_psclose
  !p.multi = [0,1,1]

  return

end
