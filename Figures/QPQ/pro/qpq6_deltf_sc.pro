;;  Compares QPQ6 deltF values against SC simulation output
pro qpq6_deltf_sc, ORIG=orig, RANDOM=random, MEDIAN=median, NODLA=nodla, $
                      PROX=prox, ZLYA=zlya, NEW=new, KLUDGE=kludge

  compile_opt strictarr

;  
  if not keyword_set(PSFILE) then psfile = 'qpq6_deltf_sc.ps'
  if not keyword_set( NTOT ) then ntot = 12L
  if not keyword_set( LSZ ) then lsz = 1.8
  if not keyword_set( LTHICK ) then lthick = 4.
  if not keyword_set( PTHICK ) then pthick = 2.
  if not keyword_set( BLSZ ) then blsz = 1.7
  if not keyword_set( CSIZE ) then csize = 2.3
  if not keyword_set( KLUDGE ) then kludge = 1.

  if not keyword_set( CONTI ) then conti = 0.8
  if not keyword_set(LIT_H) then LIT_H = 0.70  ;; (Hubble)

  vmnx = [-1000, 1000]

  qpq_fil = '~/Dropbox/QSOPairs/qpq6_final_cut.fits'
  qpq_strct = xmrdfits(qpq_fil, 1)
  gdz = where(qpq_strct.z_fg GT 2.0 and qpq_strct.z_fg LT 3.); and $
  gd_strct = qpq_strct[gdz]

  print, 'Median z: ', median(gd_strct[where(gd_strct.R_phys GT 500.)].z_fg)

  VINT = 1000.  ; km/s
  ewi = round( median( where(abs(gd_strct.vEW - VINT) LT 1.)  $
                       MOD n_elements(gd_strct[0].vEW) ) )

  Rbins = [ [30., 100.], [100, 200], [200, 300], [300, 500], [500, 1000]] ;; kpc
  szR = size(Rbins, /dimen)

  ;; Plot
  if keyword_set( PSFILE ) then x_psopen, psfile, /maxs
  !p.multi=[0,1,1]
  clr = getcolor(/load)
  fclr = clr.white

  thisletter = byte(94)
  perpletter = '!9' + string(thisletter) + '!X'

  ;; SC
  root = getenv('QSO_DIR')+'tex/QPQ6/'
  readcol, root+'/Analysis/Mass/SimDelta_Standard.dat', R_SCs, deltF_SCs
  R_SCs = 10.^R_SCs * 1e3
  readcol, root+'/Analysis/Mass/SimDelta_ColdHalo.dat', R_SCc, deltF_SCc
  R_SCc = 10.^R_SCc * 1e3

  yrng = [0., 0.5]
  plot, [0], [0], $
        xrange=[20., 2000], xtit='R!d'+perpletter+'!N (kpc)',$
        ytit='!9d!X!d<F>!N', $
        yrange=yrng, xtickn=xspaces, xmargin=[7,13], $
        ymargin=[4,1], NODATA=nblnd, $ ;ytickn=yspaces, $
        charsize=csize, psym=10, background=clr.white, color=fclr, $
        xstyle=1, ystyle=1, thick=pthick, ytickinterval=ytint, /xlog

  ;; LOOP
  for ii=0L,szR[1]-1 do begin

     ;; Data
     idxR  = where(gd_strct.R_phys GT  Rbins[0, ii] and $
                   gd_strct.R_phys LT Rbins[1,ii], npair)
     deltF = gd_strct[idxR].deltFLya[ewi]
     stat = moment(deltF)

     ;; Mean R
     plotsym, 8, 1.5, /fill 
     meanR = mean( gd_strct[idxR].R_phys )
     oploterror, [meanR], [stat[0]], [Rbins[1,ii]-meanR], $
                 [sqrt(stat[1])/npair], color=fclr, errcolo=fclr, psym=8, /hib
     oploterror, [meanR], [stat[0]], [abs(Rbins[0,ii]-meanR)], $
                 [sqrt(stat[1])/npair], color=fclr, errcolo=fclr, psym=8, /lob

  endfor

  ;; Models
  plotsym, 0, 1.5, /fill 
  oplot, R_SCs, deltF_SCs, color=clr.yellow, psym=-8
  plotsym, 4, 1.5, /fill 
  oplot, R_SCc, deltF_SCc, color=clr.cyan, psym=-8, linesty=1

  ;; Label
  xlbl = 300.
  xpt = 250.
  ssz=1.3
  plotsym, 8, ssz, /fill 
  oplot, [xpt], [0.435], psym=8, color=fclr
  xyouts, xlbl, 0.43, 'QPQ6 data', color=fclr, charsi=lsz
  plotsym, 0, ssz, /fill 
  oplot, [xpt], [0.405], psym=8, color=clr.yellow
  xyouts, xlbl, 0.4, 'Standard simulation', color=clr.yellow, charsi=lsz
  plotsym, 4, ssz, /fill 
  oplot, [xpt], [0.375], psym=8, color=clr.cyan
  xyouts, xlbl, 0.37, 'Cold model (T<10!u4!N K)', color=clr.cyan, charsi=lsz

  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]
  
  return
end
