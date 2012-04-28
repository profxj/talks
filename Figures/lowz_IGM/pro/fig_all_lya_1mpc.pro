;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro fig_all_lya_1mpc, summ_fil, ALL_GAL=all_gal, CUT_GAL=cut_gal, NOCOS=nocos

  if not keyword_set(psfile) then psfile = 'fig_all_lya_1mpc.ps'
  if keyword_set(NOCOS) then psfile = 'fig_all_lya_1mpc_nocos.ps'

  ;; HISTOGRAM
  if not keyword_set(binsL) then binsL = 0.2
  if not keyword_set(lsz) then lsz = 1.4 
  if not keyword_set(csz) then csz = 1.5 

  ;; Plot
  x_psopen,psfile,/maxs
  !p.multi=[0,1,1]

  POS=[0.10, 0.4, 0.48, 0.97]
  fig_all_rho_1mpc, CSZ=csz, LSZ=lsz, /nonhi,  POS=POS, NOCOS=nocos
  POS[0] = 0.57
  POS[2] = 0.97
  fig_lya_incidence, /nops, CSZ=csz, LSZ=lsz,  POS=POS

  x_psclose
  !p.multi=[0,1,1]

  print, 'fig_lya_1mpc: All done'

  return
end
