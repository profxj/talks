pro fig_080607, NOLBL=nolbl

  psfile = 'fig_080607.ps'
  if not keyword_set(NOLBL) then begin
     REVIEW=1 
  endif else begin
     REVIEW=0
     NOERR=1
     NOFIT=1
  endelse

  if keyword_set( PSFILE ) then x_psopen, psfile, /maxs
  !p.multi=[0,1,1]
  clr = getcolor(/load)

  ;; H2 first
  fig_080607_h2, /sub, /NOPS, /NOLY, NOFIT=nofit, NOERR=noerr, REVIEW=REVIEW;, /NOFIT, /NOERR

  ;; R400
  fig_080607_r400, /sub, /nops, REVIEW=review

  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]

  return
end
