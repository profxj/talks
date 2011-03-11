pro fig_080607

  psfile = 'fig_080607.ps'

  if keyword_set( PSFILE ) then x_psopen, psfile, /maxs
  !p.multi=[0,1,1]
  clr = getcolor(/load)

  ;; H2 first
  fig_080607_h2, /sub, /NOPS, /NOLY, /REVIEW;, /NOFIT, /NOERR

  ;; R400
  fig_080607_r400, /sub, /nops, /REVIEW

  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]

  return
end
