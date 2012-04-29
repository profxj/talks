;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; coshalo_gallery, 'J0042-1037', '358_9', IMG=img, STRCT=strct
pro coshalo_gallery, qso, gal_id, IMG=img, STRCT=strct

  if not keyword_set(psfile) then psfile = qso+'_'+gal_id+'_gallery.ps'

  ;; HISTOGRAM
  if not keyword_set(binsL) then binsL = 0.2
  if not keyword_set(lsz) then lsz = 1.4 
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
  tv, img, 0.5, 4.5, /true, xsiz=isiz, ysiz=isiz, /inches

  clr = getcolor(/load)

  ;; Keck spectrum

  red_spec_fil = getenv('DROPBOX_DIR')+ $
                 '/coshaloanalysis/fields/'+qso+'/'+gal_id+'/spec1d/'+strct.galaxy.spec1d_redcorr
  blue_spec_fil = getenv('DROPBOX_DIR')+ $
                 '/coshaloanalysis/fields/'+qso+'/'+gal_id+'/spec1d/'+strct.galaxy.spec1d_bluecorr

  blue_spec = x_readspec(blue_spec_fil, infl=2, wav=blue_wv, npix=blue_npix)
  red_spec = x_readspec(red_spec_fil, infl=2, wav=red_wv, npix=red_npix)

  xrng = [3500., 8000]
  yrng = [0.,  max( [blue_spec, red_spec]) ] * 1.05
  xtitle='Wavelength (Ang)'
  bclr = clr.black
  fclr = clr.lightgray
  
  POS=[0.45, 0.6, 0.97, 0.97]

  plot,[0.],[0.],xrange=xrng,yrange=yrng,color=fclr,$
       background=bclr, charsize=csz,  xtickn=spcs, $
       POS=pos, xtitl=xtitle, $
       xmargin=xmrg, ymargin=ymrg, ytitle='Flux', $
       /nodata,xstyle=1,ystyle=1, /noeras

  oplot, blue_wv, blue_spec, color=clr.blue, thick=3, psym=10
  oplot, red_wv, red_spec, color=clr.red, thick=3, psym=10

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
  
  POS=[0.15, 0.1, 0.97, 0.5]

  plot,[0.],[0.],xrange=vmnx,yrange=yrng,color=fclr,$
       background=bclr, charsize=csz,  xtickn=spcs, $
       POS=pos, xtitl=xtitle, $
       xmargin=xmrg, ymargin=ymrg, ytitle='Flux', $
       /nodata,xstyle=1,ystyle=1, /noeras

  oplot, velo, COS_spec, color=clr.black, thick=3, psym=10


  x_psclose
  !p.multi=[0,1,1]

  print, 'coshalo_gallery: All done'

  return
end
