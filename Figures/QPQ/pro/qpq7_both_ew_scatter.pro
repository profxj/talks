;; Show the CII and CIV EW scatter plots side-by-side
pro qpq7_both_ew_scatter

  if not keyword_set( XMRG ) then xmrg = [7,1]
  if not keyword_set( XLBL ) then xlbl = 215.
  if not keyword_set( LSZ ) then lsz = 1.3
  if not keyword_set( LTHICK ) then lthick = 4.
  if not keyword_set( PTHICK ) then pthick = 2.
  if not keyword_set( BLSZ ) then blsz = 1.7
  if not keyword_set( CSIZE ) then csize = 1.7

  if not keyword_set( CONTI ) then conti = 0.8

  if not keyword_set( RVIR ) then rvir = 160. ; kpc

  ;; QPQ7
  if not keyword_set(PSFILE) then psfile = 'qpq7_both_ew_scatter.ps'

  ;; Plot
  if keyword_set( PSFILE ) then x_psopen, psfile, /maxs
  !p.multi=[0,2,1]

  ;; 1334
  qpq7_ew_vs_rho, 1334.5323d, /nops, csize=csize, LBL='CII 1334'

  ;; 1548
  qpq7_ew_vs_rho, 1548.195d, /nops, csize=csize, LBL='CIV 1548'

  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]
  x_psclose

  return
end
