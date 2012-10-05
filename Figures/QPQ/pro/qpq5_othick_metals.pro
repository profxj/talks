;; Stacks several plots together
;; qpq5_othick_metals, psfile='qpq5_othick_metals_joe.ps', /nozfg, /neil_lbg
pro qpq5_othick_metals, LBG=lbg, NOZFG=nozfg, PSFILE=psfile, NEIL_LBG=neil_lbg

  if not keyword_set( LSZ ) then lsz = 1.8
  if not keyword_set( LTHICK ) then lthick = 4.
  if not keyword_set( PTHICK ) then pthick = 2.
  if not keyword_set( BLSZ ) then blsz = 1.7
  if not keyword_set( CSIZE ) then csize = 2.8
  if not keyword_set(XMRG) then xmrg = [7,6]

  if not keyword_set( CONTI ) then conti = 0.8

  ;; QPQ5
  qpq_fil = '~/Dropbox/QSOPairs/qpq5_pairs.fits'
  qpq5_strct = xmrdfits(qpq_fil, 1)

  if not keyword_set(PSFILE) then psfile = 'qpq5_othick_metals.ps'
  if keyword_set(LBG) then psfile = 'qpq5_othick_metals_lbg.ps'

  ;; Plot
  if keyword_set( PSFILE ) then x_psopen, psfile, /maxs
  !p.multi=[0,1,3]

  ;; Othick
  qpq5_covering, /nops, csz=csize, XMRG=xmrg, NOZFG=nozfg, neil_lbg=neil_lbg

  ;; CII EW
  qpq5_ew1334, /nops, csz=csize, XMRG=xmrg, LBG=lbg

  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]

  return
end
