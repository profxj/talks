;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; coshalo_gallery, 'J0042-1037', '358_9', IMG=img, STRCT=strct
;; coshalo_gallery, 'J0935+0204', '15_28' -- MagE!
;; coshalo_gallery, 'J1342-0053', '304_29' -- MagE!
;; coshalo_gallery, 'J1009+0713',  '204_17'
;; coshalo_gallery, 'J1555+3628', '88_11'
;; coshalo_gallery, 'J1419+4207', '132_30'
pro coshalo_gallery, qso, gal_id, IMG=img, STRCT=strct, PSFILE=psfile, SMAX=smax

  if not keyword_set(Ipsfile) then psfile = qso+'_'+gal_id+'_gallery.ps'

  ;; HISTOGRAM
  if not keyword_set(binsL) then binsL = 0.2
  if not keyword_set(lsz) then lsz = 1.4 
  if not keyword_set(lsz2) then lsz2 = 1.6 
  if not keyword_set(csz) then csz = 1.5 
  if not keyword_set(isz) then isiz = 3.4

  ;; Mega structure
  ldir = getenv('DROPBOX_DIR')+'/COS-Halos/lowions/'
  restore, ldir+'/cosmetals_megastructure.sav'

  ;; Find the galaxy
  if not keyword_set(STRCT) then begin
     idx = where(megastruct.galaxy.field EQ qso and megastruct.galaxy.galid EQ gal_id, nidx)
     if nidx NE 1 then stop
     strct = megastruct[idx]
  endif

  ;; Grab the SDSS image
  imsize = 1.
  npix = round(imsize*60./0.39612)
  x_radec, strct.galaxy.ra, strct.galaxy.dec, rad, decd
  if not keyword_set(IMG) then img = sdss_queryimage(RAD, DECD, xs=npix, ys=npix) 

  
  ;; Plot
  x_psopen,psfile,/maxs
  !p.multi=[0,1,1]

  ;; Image
  loadct, 0, /silent
  szi = size(img, /dimen)
  img[*,*,0] = 255
  img[*,0,*] = 255
  img[*,*,szi[2]-1] = 255
  img[*,szi[1]-1,*] = 255
  tv, img, 0.5, 4.5, /true, xsiz=isiz, ysiz=isiz, /inches

  clr = getcolor(/load)
  bclr = clr.black
  fclr = clr.lightgray
  lclr = clr.yellow

  xyouts, 0.20, 0.93, qso+'  '+gal_id, color=lclr, charsiz=lsz2, /normal, align=0.5
  xyouts, 0.06, 0.63, '!9r!X = '+strtrim(string(round(strct.rhofinal),format='(i3)'),2)+'kpc', $
          color=lclr, charsiz=lsz2, /normal, align=0.
  xyouts, 0.06, 0.60, 'z = '+string(strct.zfinal,format='(f5.3)'), $
          color=lclr, charsiz=lsz2, /normal, align=0.
  xyouts, 0.06, 0.57, 'M!dstar!N = 10!u'+strtrim(string(strct.logMfinal,format='(f4.1)'),2)+'!N M!dSun!N', $
          color=lclr, charsiz=lsz2, /normal, align=0.


  ;; Keck spectrum

  red_spec_fil = getenv('DROPBOX_DIR')+ $
                 '/coshaloanalysis/fields/'+qso+'/'+gal_id+'/spec1d/'+strct.galaxy.spec1d_redcorr
  blue_spec_fil = getenv('DROPBOX_DIR')+ $
                 '/coshaloanalysis/fields/'+qso+'/'+gal_id+'/spec1d/'+strct.galaxy.spec1d_bluecorr

  blue_spec = x_readspec(blue_spec_fil, infl=2, wav=blue_wv, npix=blue_npix)
  red_spec = x_readspec(red_spec_fil, infl=2, wav=red_wv, npix=red_npix)

  xrng = [3500., 8000]
  yrng = [0.,  max( [blue_spec, red_spec]) ] * 1.05
  if keyword_set(SMAX) then yrng = [0., smax]
  xtitle='Wavelength (Ang)'
  
  POS=[0.45, 0.6, 0.97, 0.97]

  plot,[0.],[0.],xrange=xrng,yrange=yrng,color=fclr,$
       background=bclr, charsize=csz,  xtickn=spcs, $
       POS=pos, xtitl=xtitle, $
       xmargin=xmrg, ymargin=ymrg, ytitle='Flux', $
       /nodata,xstyle=1,ystyle=1, /noeras

  oplot, blue_wv, blue_spec, color=clr.cyan, thick=3, psym=10
  oplot, red_wv, red_spec, color=clr.tomato, thick=3, psym=10

  xyouts, 5000., yrng[1]*0.9, 'Keck/LRIS ', color=lclr, charsi=lsz2, align=0.
  xyouts, 5040., yrng[1]*0.82, '(Werk+11)', color=lclr, charsi=lsz, align=0.

  ;; HST/COS
  ;; Lya line 
  COS_spec_fil = getenv('DROPBOX_DIR')+ $
                 '/COS-Halos/Targets/'+qso+'/Data/'+qso+'_nbin3_coadd.fits'

  COS_spec = x_readspec(COS_spec_fil, infl=3, wav=COS_wv, npix=COS_npix)
  vmnx = [-1000., 1000]

  x_pixminmax, COS_wv, 1215.6701, strct.zfinal, vmnx[0], vmnx[1], PIXMIN=pmn, $
                       PIXMAX=pmx, VELO=velo
  

  mx =  max( COS_spec[pmn:pmx] ) 
  yrng = [-0.1*mx,  mx*1.1]
  xtitle='Velocity (km/s)'
  bclr = clr.black
  fclr = clr.lightgray
  
  POS=[0.10, 0.1, 0.97, 0.5]

  plot,[0.],[0.],xrange=vmnx,yrange=yrng,color=fclr,$
       background=bclr, charsize=csz,  xtickn=spcs, $
       POS=pos, xtitl=xtitle, $
       xmargin=xmrg, ymargin=ymrg, ytitle='Flux', $
       /nodata,xstyle=1,ystyle=1, /noeras

  oplot, velo, COS_spec, color=fclr, thick=4, psym=10

  oplot, [0., 0], yrng, color=clr.yellow, linesty=1

  xyouts, vmnx[0]+ 0.05*(vmnx[1]-vmnx[0]),  yrng[0]+0.15*(yrng[1]-yrng[0]), $
          'HST/COS', color=lclr, charsi=lsz2, align=0.
  xyouts, vmnx[0]+ 0.06*(vmnx[1]-vmnx[0]),  yrng[0]+0.07*(yrng[1]-yrng[0]), $
          ' (Thom/Tumlinson+12)', color=lclr, charsi=lsz, align=0.

  x_psclose
  !p.multi=[0,1,1]

  print, 'coshalo_gallery: All done'

  return
end
