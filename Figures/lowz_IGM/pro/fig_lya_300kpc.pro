;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro fig_lya_300kpc, summ_fil, STRUCT=struct, PSFILE=psfile, NOPS=nops

  if not keyword_set(psfile) then psfile = 'fig_lya_300kpc.ps'
  if not keyword_set(binsz) then binsz = 0.01
  if not keyword_set(lsz) then lsz = 2.0 
  if not keyword_set(csz) then csz = 1.7 
  if not keyword_set(ssiz) then ssiz=2.
  if not keyword_set(dlim) then dlim = 5. ;arcmin
  if not keyword_set(RHO_GAL) then rho_gal = 150. ; kpc

  root = '~/paper/OVI/Galaxies/'
  all_gal = xmrdfits(root+'Analysis/all_galabs_strct.fits',1)
  summ_fil = ['Analysis/dwarf_galabs_strct.fits', $
              'Analysis/subls_galabs_strct.fits', $
              'Analysis/lstar_galabs_strct.fits']

  ;; Plot
  if keyword_set(psfile) then x_psopen,psfile,/maxs
  !p.multi=[0,2,1]
  clr = getcolor(/load)
  nclr = clr.darkgray
  lclr = clr.lightgray

  ;; Input
  for tt=0,2 do begin
     struct = xmrdfits(root+summ_fil[tt],1)
     ngal = n_elements(struct)
     case tt of
        0: begin
           pclr = clr.yellow
;           pclr = clr.darkgreen
           glbl='Dwarfs (L<0.1L*)' 
        end
        1: begin 
;           pclr = clr.blue
           pclr = clr.cyan
           glbl = 'Sub-L* (0.1L* < L < L*)'
        end
        2: begin
;           pclr = clr.red
           pclr = clr.tomato
           glbl='L* Galaxies (L>L*)' 
        end
     endcase


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Lya
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; EW
  xrng = [10., 400.]
  yrng = [30., 3000] ; mA
  xmrg = [8,0.5]
  ymrg = [12,0.5]
  ppsym = 1
  
  plot,[0.],[0.],xrange=xrng,yrange=yrng,color=lclr, $
       background=clr.white, charsize=csz, thick=9, $
       xmargin=xmrg, ymargin=ymrg,ytitle='W!uLy!9a!X!N (mA)',$
       xtitle='!9r!X (kpc)',/nodata,xstyle=1,ystyle=1, /xlog, /ylog
  
;  xyouts, 12., 40.*(2.^tt), glbl, color=pclr, charsiz=lsz, align=0.

  ;; Limit?
  lim = where(struct.mag[2] GT 0 AND struct.mag[5] GT 99. $
              and strmatch(strtrim(struct.gal_type,2), 'Early'), nlim) 
  plotsym, 1, 2.9, thick=7
  if nlim GT 0 then $
     oplot, [struct[lim].dra], [struct[lim].mag[4]], color=pclr, psym=8
  lim = where(struct.mag[2] GT 0 AND struct.mag[5] GT 99. $
              and strmatch(strtrim(struct.gal_type,2), 'Late'), nlim) 
  plotsym, 1, 2.3, thick=7
  if nlim GT 0 then $
     oplot, [struct[lim].dra], [struct[lim].mag[4]], color=pclr, psym=8

  ;; Value
  val = where(struct.mag[2] GT 0 AND struct.mag[5] LT 99. $
              and strmatch(strtrim(struct.gal_type,2), 'Early'), nval) 
  plotsym, 0, 2.0, thick=7
  if nval GT 0 then oplot, [struct[val].dra], [struct[val].mag[4]], color=pclr, psym=8
  val = where(struct.mag[2] GT 0 AND struct.mag[5] LT 99. $
              and strmatch(strtrim(struct.gal_type,2), 'Late'), nval) 
  if nval GT 0 then oplot, [struct[val].dra], [struct[val].mag[4]], color=pclr, psym=2, symsiz=ssiz

  ;; Nearby galaxies
  dwarfs = where(struct.mag[2] GT 0, ndwarf)
  for ii=0L,ndwarf-1 do begin
     idx = dwarfs[ii]
     dv = x_relvel(struct[idx].z, all_gal.z, /rev)
     mtch = where(all_gal.dra LT RHO_GAL AND $
                  strmatch(strtrim(all_gal.field,2),strtrim(struct[idx].field,2)) AND $
                  (all_gal.ddec GT 0.1 AND all_gal.ddec LT 1.) AND $
                  abs(dv) LT 300., nmtch)
     if nmtch GT 0 then begin
        plotsym, 8, 1.5, thick=4
        ;oplot, [struct[idx].dra], [struct[idx].mag[4]], color=nclr, psym=8
        struct[idx].area = -1
     endif
     mtch = where(all_gal.dra LT RHO_GAL AND $
                  strmatch(strtrim(all_gal.field,2),strtrim(struct[idx].field,2)) AND $
                  all_gal.ddec GT 1.  AND $
                  abs(dv) LT 300., nmtch)
     if nmtch GT 0 then begin ;; L*
;        plotsym, 4, 1.3, thick=4
        ;oplot, [struct[idx].dra], [struct[idx].mag[4]], color=nclr, psym=4, symsi=1.7
        struct[idx].area = -1
     endif
  endfor

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; NHI
  xrng = [10., 400.]
  yrng = [12.7, 16] 
  xmrg = [8,0.5]
  
  plot,[0.],[0.],xrange=xrng,yrange=yrng,color=lclr, $
       background=clr.white, charsize=csz, thick=9, $
       xmargin=xmrg, ymargin=ymrg,ytitle='log N!dHI!N', $
       xtitle='!9r!X (kpc)',/nodata,xstyle=1,ystyle=1, /xlog

  ;; Upper Limit?
  lim = where(struct.mag[2] GT 0 AND struct.mag[9] LE 0. $
              and strmatch(strtrim(struct.gal_type,2), 'Early'), nlim) 
  plotsym, 1, 2.9, thick=5
  if nlim GT 0 then oplot, [struct[lim].dra], [struct[lim].mag[8]], color=pclr, psym=8
  lim = where(struct.mag[2] GT 0 AND struct.mag[9] LE 0. $
              and strmatch(strtrim(struct.gal_type,2), 'Late'), nlim) 
  plotsym, 1, 2.3, thick=5
  if nlim GT 0 then oplot, [struct[lim].dra], [struct[lim].mag[8]], color=pclr, psym=8

  ;; Lower Limit?
  lim = where(struct.mag[2] GT 0 AND (struct.mag[9] GT 9. OR $
                                      struct.mag[8] GT yrng[1]), nlim) 
  plotsym, 2, 2.9, thick=5
  if nlim GT 0 then oplot, [struct[lim].dra], [struct[lim].mag[8] < (yrng[1]-0.2)], $
                           color=pclr, psym=8

  ;; Value
  val = where(struct.mag[2] GT 0 AND (struct.mag[9] LT 9. and struct.mag[9] GT 0) $
              and strmatch(strtrim(struct.gal_type,2), 'Early'), nval) 
  plotsym, 0, 2.0, thick=4
  if nval GT 0 then oplot, [struct[val].dra], [struct[val].mag[8]], color=pclr, psym=8
  val = where(struct.mag[2] GT 0 AND (struct.mag[9] LT 9. and struct.mag[9] GT 0) $
              and strmatch(strtrim(struct.gal_type,2), 'Late'), nval) 
  if nval GT 0 then oplot, [struct[val].dra], [struct[val].mag[8]], color=pclr, psym=2, $
                           symsiz=ssiz
  val = where(struct.mag[2] GT 0 AND (struct.mag[9] LT 9. and struct.mag[9] GT 0) $
              and strmatch(strtrim(struct.gal_type,2), 'Unkn'), nval) 
;  if nval GT 0 then oplot, [struct[val].dra], [struct[val].mag[8]], color=pclr, psym=5

  ;; Nearby galaxies
  dwarfs = where(struct.mag[9] GT 0, ndwarf)
  for ii=0L,ndwarf-1 do begin
     idx = dwarfs[ii]
     dv = x_relvel(struct[idx].z, all_gal.z, /rev)
     mtch = where(all_gal.dra LT RHO_GAL AND $
                  all_gal.ddec GT 1.  AND $
                  abs(dv) LT 300., nmtch)
     if nmtch GT 0 then begin
;        plotsym, 0, 1.3, thick=4
;        oplot, [struct[idx].dra], [struct[idx].mag[8]], color=nclr, psym=4, symsi=1.5
        struct[idx].area = -1
     endif
  endfor


endfor

  x_psclose
  !p.multi=[0,1,1]
  return

end
