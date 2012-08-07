pro fig_define_igm

 ;; 3D plot

  if not keyword_set( PSFILE ) then psfile = 'fig_define_igm.ps'
  if not keyword_set( CSZ ) then csz = 1.9
  if not keyword_set( LSZ ) then lsz = 1.9
  if not keyword_set( BLSZ ) then blsz = 2.0

  ;x_psopen, psfile, /maxs
  ;clr = getcolor(/load)
  ;!p.multi=[0,1,1]

  ;; ;;;;;;;;;;;;;;;;;
  ;; ISM

  aISMp = [3,6]
  aISMR = [0,1]
  aISMT = [1,4]

  ISMp = [aISMp[0], $
          replicate(aISMp[1],2), $
          replicate(aISMp[0],3), $
          replicate(aISMp[1],2), $
          replicate(aISMp[0],2), $
          replicate(aISMp[1],3), $
          replicate(aISMp[0],2)]

  ISMR = [replicate(aISMR[0], 2), $
          replicate(aISMR[1], 2), $
          replicate(aISMR[0], 3), $
          replicate(aISMR[1], 2), $
          replicate(aISMR[0], 3), $
          replicate(aISMR[1], 3)]

  ISMT = [ replicate(aISMT[0],5), $
           replicate(aISMT[1],6), $
           replicate(aISMT[0],3), $
           aISMT[1]]
  
  
  xrange = [-5,6]
  zrange = [1,8]
  p = PLOT3D(ISMp, ISMR, ISMT, '-r3o', /SYM_FILLED, $
             XRANGE=xrange, ZRANGE=zrange, $
             YRANGE=[0,5], charsiz=0.5, $
             AXIS_STYLE=2, MARGIN=[0.2,0.3,0.1,0],$
             XMINOR=0, YMINOR=0, ZMINOR=0, $
             DEPTH_CUE=[0,2], /PERSPECTIVE, $
             COLOR="red", $
             XTHICK=4, YTHICK=4, ZTHICK=4, $
             SHADOW_COLOR="deep sky blue", $
             ;XY_SHADOW=1, YZ_SHADOW=1, XZ_SHADOW=1,$
             XTITLE='log !9r/<r>!X', YTITLE='log R (kpc)', ZTITLE='log T (K)')

  ;; ;;;;;;;;;;;;;;;;;
  ;; CGM

  aCGMp = [1,3]
  aCGMR = [1,2.5]
  aCGMT = [4,6]

  CGMp = [aCGMp[0], $
          replicate(aCGMp[1],2), $
          replicate(aCGMp[0],3), $
          replicate(aCGMp[1],2), $
          replicate(aCGMp[0],2), $
          replicate(aCGMp[1],3), $
          replicate(aCGMp[0],2)]

  CGMR = [replicate(aCGMR[0], 2), $
          replicate(aCGMR[1], 2), $
          replicate(aCGMR[0], 3), $
          replicate(aCGMR[1], 2), $
          replicate(aCGMR[0], 3), $
          replicate(aCGMR[1], 3)]

  CGMT = [ replicate(aCGMT[0],5), $
           replicate(aCGMT[1],6), $
           replicate(aCGMT[0],3), $
           aCGMT[1]]
  
  
  p = PLOT3D(CGMp, CGMR, CGMT, '-b3o', /SYM_FILLED, $
             /overplot, /current, $
             XRANGE=xrange, ZRANGE=zrange, $
             YRANGE=[0,5], charsiz=1.0, $
             MARGIN=[0.2,0.3,0.1,0],$
             XMINOR=0, YMINOR=0, ZMINOR=0, $
             DEPTH_CUE=[0,2], /PERSPECTIVE, $
             COLOR="blue", $
             SHADOW_COLOR="deep sky blue")
             ;XY_SHADOW=1, YZ_SHADOW=1, XZ_SHADOW=1,$

  ;; ;;;;;;;;;;;;;;;;;
  ;; ICM

  aICMp = [1,2]
  aICMR = [2,3]
  aICMT = [6,8]

  ICMp = [aICMp[0], $
          replicate(aICMp[1],2), $
          replicate(aICMp[0],3), $
          replicate(aICMp[1],2), $
          replicate(aICMp[0],2), $
          replicate(aICMp[1],3), $
          replicate(aICMp[0],2)]

  ICMR = [replicate(aICMR[0], 2), $
          replicate(aICMR[1], 2), $
          replicate(aICMR[0], 3), $
          replicate(aICMR[1], 2), $
          replicate(aICMR[0], 3), $
          replicate(aICMR[1], 3)]

  ICMT = [ replicate(aICMT[0],5), $
           replicate(aICMT[1],6), $
           replicate(aICMT[0],3), $
           aICMT[1]]
  
  
  p = PLOT3D(ICMp, ICMR, ICMT, '-r3o', /SYM_FILLED, $
             /overplot, /current, $
             XRANGE=xrange, ZRANGE=zrange, $
             YRANGE=[0,5], charsiz=1.0, $
             MARGIN=[0.2,0.3,0.1,0],$
             XMINOR=0, YMINOR=0, ZMINOR=0, $
             DEPTH_CUE=[0,2], /PERSPECTIVE, $
             COLOR="green", $
             SHADOW_COLOR="deep sky blue")
             ;XY_SHADOW=1, YZ_SHADOW=1, XZ_SHADOW=1,$

  ;; ;;;;;;;;;;;;;;;;;
  ;; IGM

  aIGMp = [-5,1]
  aIGMR = [2.5,4]
  aIGMT = [4,5]

  IGMp = [aIGMp[0], $
          replicate(aIGMp[1],2), $
          replicate(aIGMp[0],3), $
          replicate(aIGMp[1],2), $
          replicate(aIGMp[0],2), $
          replicate(aIGMp[1],3), $
          replicate(aIGMp[0],2)]

  IGMR = [replicate(aIGMR[0], 2), $
          replicate(aIGMR[1], 2), $
          replicate(aIGMR[0], 3), $
          replicate(aIGMR[1], 2), $
          replicate(aIGMR[0], 3), $
          replicate(aIGMR[1], 3)]

  IGMT = [ replicate(aIGMT[0],5), $
           replicate(aIGMT[1],6), $
           replicate(aIGMT[0],3), $
           aIGMT[1]]
  
  
  p = PLOT3D(IGMp, IGMR, IGMT, '-r3o', /SYM_FILLED, $
             /overplot, /current, $
             XRANGE=xrange, ZRANGE=zrange, $
             YRANGE=[0,5], charsiz=1.0, $
             MARGIN=[0.2,0.3,0.1,0],$
             XMINOR=0, YMINOR=0, ZMINOR=0, $
             DEPTH_CUE=[0,2], /PERSPECTIVE, $
             COLOR="cyan", $
             SHADOW_COLOR="deep sky blue")
             ;XY_SHADOW=1, YZ_SHADOW=1, XZ_SHADOW=1,$

  ;; ;;;;;;;;;;;;;;;;;
  ;; WHIM

  aWHIMp = [-5,1]
  aWHIMR = [3,5]
  aWHIMT = [5,7]

  WHIMp = [aWHIMp[0], $
          replicate(aWHIMp[1],2), $
          replicate(aWHIMp[0],3), $
          replicate(aWHIMp[1],2), $
          replicate(aWHIMp[0],2), $
          replicate(aWHIMp[1],3), $
          replicate(aWHIMp[0],2)]

  WHIMR = [replicate(aWHIMR[0], 2), $
          replicate(aWHIMR[1], 2), $
          replicate(aWHIMR[0], 3), $
          replicate(aWHIMR[1], 2), $
          replicate(aWHIMR[0], 3), $
          replicate(aWHIMR[1], 3)]

  WHIMT = [ replicate(aWHIMT[0],5), $
           replicate(aWHIMT[1],6), $
           replicate(aWHIMT[0],3), $
           aWHIMT[1]]
  
  
  p = PLOT3D(WHIMp, WHIMR, WHIMT, '-r3o', /SYM_FILLED, $
             /overplot, /current, $
             XRANGE=xrange, ZRANGE=zrange, $
             YRANGE=[0,5], charsiz=1.0, $
             MARGIN=[0.2,0.3,0.1,0],$
             XMINOR=0, YMINOR=0, ZMINOR=0, $
             DEPTH_CUE=[0,2], /PERSPECTIVE, $
             COLOR="orange", $
             SHADOW_COLOR="deep sky blue")
             ;XY_SHADOW=1, YZ_SHADOW=1, XZ_SHADOW=1,$

  return
end
