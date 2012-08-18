;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro fig_lstar_ovi_1mpc, summ_fil, STRUCT=struct, PSFILE=psfile, NEAR=near, NOPS=nops

  if not keyword_set(psfile) then psfile = 'fig_lstar_1mpc_ovi.ps'
  if not keyword_set(binsz) then binsz = 0.01
  if not keyword_set(lsz) then lsz = 2.2 
  if not keyword_set(csz) then csz = 2.4 
  if not keyword_set(dlim) then dlim = 5. ;arcmin
  if not keyword_set(RHO_GAL) then rho_gal = 300. ; kpc

  ;; Input
  if not keyword_set(SUMM_FIL) then $
     summ_fil = '~/paper/OVI/Galaxies/Analysis/lstar_galabs_1Mpc_strct.fits'
  struct = xmrdfits(summ_fil,1)
  ngal = n_elements(struct)
  all_gal = xmrdfits('~/paper/OVI/Galaxies/Analysis/all_galabs_strct.fits',1)


  ;; Plot
  if keyword_set(psfile) then x_psopen,psfile,/maxs
  !p.multi=[0,1,1]
  clr = getcolor(/load)
  nclr = clr.green
  pclr = clr.tomato
  ;lclr = clr.white
  lclr = clr.black
  ppsym = 4

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; OVI
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; EW
  xrng = [10., 1100.]
  yrng = [10., 1000] ; mA
  xmrg = [8,2.0]
  ymrg = [4,0.5]
  
  plot,[0.],[0.],xrange=xrng,yrange=yrng,color=lclr,$
       background=clr.white, charsize=csz, thick=9, $
       xmargin=xmrg, ymargin=ymrg,ytitle='W!u1031!N (mA)',$
       xtitle='!9r!X (kpc)',/nodata,xstyle=1,ystyle=1, /xlog, /ylog

;  xyouts, xrng[0]*1.1, 15., '(b)', color=lclr, charsiz=lsz, align=0.

  ;; Limit?
  lim = where(struct.magerr[2] GT 0 AND struct.magerr[5] GT 99. $
              and strmatch(strtrim(struct.gal_type,2), 'Early'), nlim) 
  plotsym, 1, 2.9, thick=7
  if nlim GT 0 then oplot, [struct[lim].dra], [struct[lim].magerr[4]], color=pclr, psym=8
  lim = where(struct.magerr[2] GT 0 AND struct.magerr[5] GT 99. $
              and strmatch(strtrim(struct.gal_type,2), 'Late'), nlim) 
  plotsym, 1, 2.3, thick=7
  if nlim GT 0 then oplot, [struct[lim].dra], [struct[lim].magerr[4]], color=pclr, psym=8

  ;; Value
  val = where(struct.magerr[2] GT 0 AND struct.magerr[5] LT 99. $
              and strmatch(strtrim(struct.gal_type,2), 'Early'), nval) 
  plotsym, 0, 2.0, thick=7
  if nval GT 0 then oplot, [struct[val].dra], [struct[val].magerr[4]], color=pclr, psym=8
  val = where(struct.magerr[2] GT 0 AND struct.magerr[5] LT 99. $
              and strmatch(strtrim(struct.gal_type,2), 'Late'), nval) 
  if nval GT 0 then oplot, [struct[val].dra], [struct[val].magerr[4]], color=pclr, psym=2, $
                           symsiz=2.3

  ;; Nearby galaxies
  if keyword_set(NEAR) then begin
     dwarfs = where(struct.mag[2] GT 0, ndwarf)
     for ii=0L,ndwarf-1 do begin
        idx = dwarfs[ii]
        dv = x_relvel(struct[idx].z, all_gal.z, /rev)
        ;; Sub-L*
        mtch = where(all_gal.dra LT struct[idx].dra AND $
                     strmatch(strtrim(all_gal.field,2),strtrim(struct[idx].field,2)) AND $
                     (all_gal.ddec GT 0.1 AND all_gal.ddec LT 1.) AND $
                     abs(dv) LT 300., nmtch)
        if nmtch GT 0 then begin
           plotsym, 8, 1.7, thick=4
           oplot, [struct[idx].dra], [struct[idx].magerr[4]], color=nclr, psym=8
           struct[idx].area = -1
        endif
        ;; Dwarfs
        mtch = where(all_gal.dra LT struct[idx].dra AND $
                     strmatch(strtrim(all_gal.field,2),strtrim(struct[idx].field,2)) AND $
                     all_gal.ddec LT 0.1  AND $
                     abs(dv) LT 300., nmtch)
        if nmtch GT 0 then begin
           plotsym, 4, 1.7, thick=4
           oplot, [struct[idx].dra], [struct[idx].magerr[4]], color=nclr, psym=8
           struct[idx].area = -1
        endif
     endfor
  endif

  ;; Input COS-Halos
  fmt = 'a,a,f,i,f,f,f,f,f,f,f,f,f,i,a,i,i,f,f,i,i,a'
  readcol, '~/Dropbox/COS-Halos/Analysis/OVI_EW_rho_color',f=fmt,comment='#',qso,galname, $
           zgal,vis,u,g,r,mu,mg,mr,umr,mstar,zqso,rho,flg,ew,ewerr, n, nerr, velwidth, line, tbflg

  ;; Hits
  detect = where((strmatch(flg,'T') OR strmatch(flg,'H')) AND $
                line EQ 1031)
  oplot, rho[detect], EW[detect], color=clr.tomato, psym=4, symsiz=2.5

  ;; Misses
  miss = where((strmatch(flg,'M') OR strmatch(flg,'P')) AND $
                line EQ 1031, nmiss)
  plotsym, 1, 2.3, thick=7
  oplot, rho[miss], 2*EWerr[miss], color=clr.tomato, psym=8

;  xyouts, 15., 50., 'L* Galaxies', color=clr.tomato, charsi=lsz, align=0.

  x_psclose
  !p.multi=[0,1,1]

  return

end
