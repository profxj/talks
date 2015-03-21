;+
;procedure che fa il ps di una simulation e aggiunge color bar
;range--> range of colorbar
;zoom is the the factor to zoom in a region 
;eps for encapsulated
;ps is the plate scale
;print save output in ps
;log for log scale
;zscale --> set zscale, if not set defualt if min max
;text allows label
;tlocation --> where to put the label
;colormap  allows to change the colormap
;grey     -->set to colormap 0
;reverse  reverse color map
;format  set the format for the colorbar ('F4.3')
;title -- display title
;coltitle --> title for the color bar
;colpos --> position for the color bar
;rotate --> 0-7 angle parameter for  rotation
;smooth --> the width of the smmothing function
;fits--> if set, im1 is passed as an array
;nocolbar  --> disable colorbar 
;x/ytitle   --> put a title in the bottom/left of the plot
;white     --> if set, the text is printed in black and background is
;              in white
;left      --> move the label to the left
;black     --> force black background
;circ      --> circle or radius pix, centered
;
;-
;; EW maps for CII in QPQ, LBG, COS-Halos
PRO cos_lya_map, wrest, RANGE=range, EPS=eps,$
                  PS=ps, print=print, log=log, zscale=zscale,  units=units, text=text, tlocation=tlocation,$
                  colormap=colormap,grey=grey,  reverse=reverse, format=format, title=title,$
                  coltitle=coltitle, rotate=rotate, PSFILE=psfile, $
                  smooth=smooth, fits=fits, nocolbar=nocolbar, xtitle=xtitle, ytitle=ytitle,$
                  white=white, left=left, black=black, _extra=extra, annotatecolor=annotatecolor, circ=circ, $
                  SIG_CLM=sig_clm, SIG_FINE=sig_fine, MEGASTRUCT=megastruct

  compile_opt strictarr

  ;; Initialize
  wrest = 1215.6701d
  range = [0., 2]
  sig_fine = 0.05
  BLACK = 1
  getion, wrest, ioni, elm, Z=ionz, NM=strioni
  transc = elm+strioni+strtrim(fix(wrest),2)
  if not keyword_set( PSFILE ) then psfile = 'cos_lya_map.ps'

  if not keyword_set(RANGE) then range = [0., 5]

  if not keyword_set(LOW_BOOST) then low_boost = 1.  ;; Boost saturated
  if not keyword_set(UPPER_BOOST) then upper_boost = 1.  ;; Drop lower limit

  if not keyword_set(SIG_CLM) then sig_clm = 0.05 ; Ang
  if not keyword_set(SIG_FINE) then sig_fine = 0.05 ; dex

  if not keyword_set(SAMP_LIM) then samp_lim = 6L  ;; Minimum number for stats

  ;; Dimensions
  grid_sz = 200. ; kpc
  if not keyword_set(PIX_SZ) then pix_sz = 5. ; kpc


  ;; set color map
  if ~keyword_set(colormap) then colormap=33
  if keyword_set(grey) then colormap=0
  if ~keyword_set(format) then format='(F5.2)'
  if keyword_set(white) then colormap=39
  if ~keyword_set(annotatecolor) then annotatecolor='white'

  if ~keyword_set(csz) then csz = 1.3
  if ~keyword_set(reverse) then ctload, colormap
  if keyword_set(reverse) then ctload, colormap, /reverse

  IF keyword_set(PSFILE) THEN x_psopen, psfile,/maxs, /encapsulated,/color
  ;;set scales
  IF KEYWORD_SET(range) THEN scale=range ELSE stop ; scale=[MIN(fits),MAX(fits)]
  IF KEYWORD_SET(zscale) THEN scale=ZSCALE_RANGE(fits)

  ;; Positions
  seed = -32188L

  xwid = 0.29
  x0 = 0.02
  for ss=0,1 do begin

     ;pos = [x0 + xwid*(ss+1) + 0.01*(ss+1), 0.50, x0 + 0.01*(ss+1) + xwid*((ss+1)+1),0.92]
     pos = [x0 + xwid*(ss+2) + 0.01*(ss+2), 0.50, x0 + 0.01*(ss+2) + xwid*((ss+2)+1),0.92]
     if ss GT 0 then begin
        pos = [x0 + xwid*(ss+1) + 0.01*(ss+1), 0.13, x0 + 0.01*(ss+1) + xwid*((ss+1)+1),0.50]
        ;pos = [x0 + xwid*(ss-1) + 0.01*(ss-1), 0.13, x0 + 0.01*(ss-1) + xwid*((ss-1)+1),0.50]
     endif
     case ss of 
        0: begin ;; COS-Halos
           lbl = 'COS-Halos'
           rvir = 160. ; kpc
           rsty = 1
           ;; Initialize COS-Halos
           cd, getenv('PCOSHALOS')+'/Analysis/pro/', curr=curr
           RESOLVE_ROUTINE, 'coshalos_lowmetals_initparm', /compile_full, /either
           RESOLVE_ROUTINE, 'coshalos_cut_mega', /compile_full, /either
           RESOLVE_ROUTINE, 'coshalos_calc_lowcolm', /compile_full, /either
           coshalos_lowmetals_initparm, coshalos_init
           cd, curr
           
           ldir = getenv('DROPBOX_DIR')+'/COS-Halos/lowions/'
           tdir = getenv('DROPBOX_DIR')+'/COS-Halos/Targets/'
           
           all_mega_fil = ldir+coshalos_init.all_mega_fil 
           if not keyword_set(MEGASTRUCT) then begin
              restore, all_mega_fil 
              ;; Cut on luminosity and R
              megastruct = coshalos_cut_mega(megastruct)
           endif
           
           ;; Parse
           cos_halos_redblue, megastruct, red, blue, NRED=nred, NBLUE=nblue, FLAG=0
           cos_ew = cos_mega_parse(wrest, mega=megastruct)

           ;; Cut
           gdcos = where(cos_ew.flgEW GT 0, ngd)
           rho_obs = megastruct[gdcos].rhofinal
           ALL_EW = cos_ew[gdcos].trueEW/1e3
           flg_obs = replicate(1,ngd)
           ;; Binning
           rbins = [ [20., 50], $
                     [50, 100], $ 
                     [100, 160]] 
        end
        1: begin ;; COS-Halos
           lbl = 'COS-Dwarfs'
           rvir = 190. ; kpc
           rsty = 0
           
           ;; Read in data
           ldir = getenv('DROPBOX_DIR')+'/COS-Dwarfs/lowions/'
           restore, ldir+'/cosmetals_megastructure.sav' 
           
           ;; Parse
           dwf_ew = cos_mega_parse(wrest, mega=megastruct)

           ;; Cut
           gddwf = where(dwf_ew.flgEW GT 0, ngd)
           rho_obs = megastruct[gddwf].rhofinal
           ALL_EW = dwf_ew[gddwf].trueEW/1e3
           flg_obs = replicate(1,ngd)

           ;; Binning
           rbins = [ [20., 50], $
                     [50, 80], $ 
                     [80, 110], $ 
                     [110, 160]] 
        end
        2: begin ;; LBG
           lbl = 'LBG'
           rvir = 100. ; kpc
           rsty = 0
           ;; Kludging from Steidel+10 and Rakic+12
           R_LBG = [0.031, 0.063, 0.0905958,     0.152221,     0.215018,     0.30372,  0.429, $ 
                    0.606002,     0.856001,      1.20913] * 1e3 ; kpc
           EW_LBG = [2.01, 1.23, 1.0730036,      0.41661176,      0.65806370,      0.44781657, $ ; Ang
                     0.42758677,      0.39706513,      0.45419831,      0.40758605]
           sig_LBG = 50. ;; Percent!
           nlbg = 200L
           rho_obs = 30. + (200-30.)*findgen(nlbg)/(nlbg-1)

           ;; Spline
           splin = spl_init(R_LBG, EW_LBG, /double)
           ALL_EW = spl_interp(R_LBG, EW_LBG, splin, rho_obs)
           SIG_EW = ALL_EW * (sig_LBG/100)
           ALL_EW = ALL_EW + SIG_EW * randomn(seed, nlbg)
           flg_obs = replicate(1, nlbg)

           ;; Binning
           rbins = [ [30., 50], $
                     [50, 70], $ 
                     [70, 90], $ 
                     [90, 110], $ 
                     [110, 130], $ 
                     [130, 150], $ 
                     [150, 170], $ 
                     [170, 200]] 
        end
        3: begin ;; QPQ
           lbl = 'QPQ'
           rvir = 150. ; kpc
           rsty = 0
           ;;
           ;; Read in pair structure
           qpq_fil = '~/Dropbox/QSOPairs/qpq6_final_cut.fits'
           qpq_strct = xmrdfits(qpq_fil, 1)
           
           ;; Redshift cut
           gdz = where(qpq_strct.z_fg LT 3. and qpq_strct.flg_EWlya GT 0)
           gd_strct = qpq_strct[gdz]
           npair = n_elements(gd_strct)
           print, 'Npair = ', npair
           
           lowR = where(gd_strct.R_phys LE 200., nlow)
           print, 'N (R<200) = ', nlow
           
           ;; Push into this programs format
           rho_obs = gd_strct[lowR].R_phys ; kpc
           ALL_EW = gd_strct[lowR].EWLya
           flg_obs = replicate(1,nlow)
           ;stop

           ;; Binning
           rbins = [ [30., 100], $
                     [100, 150], $ 
                     [100, 200]] 
        end
     endcase

     ;; Generate the "Image"
     npix = round(grid_sz/pix_sz) * 2 + 1
     low_fits = fltarr(npix,npix)
     xval = pix_sz*findgen(npix) - grid_sz + 0.5 * pix_sz
     yval = pix_sz*findgen(npix) - grid_sz + 0.5 * pix_sz
     
     xgrid = xval # replicate(1., npix)
     ygrid = replicate(1., npix) # yval
     rgrid = sqrt(xgrid^2 + ygrid^2)

     nbins = (size(rbins,/dimen))[1]
     
     for qq=0L,nbins-1 do begin
        ;; Grid
        gdr = where( rgrid GE rbins[0,qq] and rgrid LT rbins[1,qq], ngdr )
        if ngdr EQ 0 then continue
        
        ;; Observation
        gdo = where( rho_obs GE rbins[0,qq] and rho_obs LT rbins[1,qq], ngdo )
        
        if ngdo LT SAMP_LIM then stop ;; Too few data
        print, 'N samp = ', ngdo
        sub_colm = ALL_EW[gdo]
        sub_flg = flg_obs[gdo]
        
        ;; Random drawing
        rani = fix( randomu(seed, ngdr) * ngdo ) 
        ran_clm = sub_colm[ rani ]
        ran_flg = sub_flg[ rani ]
        
        ;; Lower limits
        low = where(ran_flg EQ 2, nlow)
        if nlow GT 0 then ran_clm[low] = ran_clm[low] + LOW_BOOST * randomu(seed, nlow)
        
        ;; Upper limits
        upper = where(ran_flg EQ 3, nupper)
        if nupper GT 0 then ran_clm[upper] = ran_clm[upper] - UPPER_BOOST * randomu(seed, nupper)
        
        ;; Add an error
        ran_clm = ran_clm + SIG_CLM * randomn(seed, ngdr)
        
        ;; Fill in
        low_fits[gdr] = ran_clm
     endfor
     
     ;; Remap to a 1kpc grid and add scatter
     new_bin = round(pix_sz)
     fits = rebin( low_fits, npix*new_bin, npix*new_bin)
     
     ;; Add scatter
     nfine = n_elements(fits)
     fits = fits + SIG_FINE * randomn(seed, nfine)
     
     ;; Smooth?
     if KEYWORD_SET(black) THEN TVLCT, [[0], [0], [0]], 0
      
     !x.style=1
     !y.style=1
  
     IMDISP, fits, RANGE=scale, /axis, XTickformat='(A1)', YTickformat='(A1)',$
             XTICKLEN=1D-5,YTICKLEN=1D-5, title=title, xtitle=xtitle , ytitle=ytitle,$
             charsize=2.5, position=pos, margin=0.
  
     ;; rvir (or 1/2 rvir)
     if keyword_set(SHOWCIRC) then begin
        x_oplotcirc, rvir, x0=((size(fits))[1])*0.5, y0=((size(fits))[2])*0.5, $
                     color=fsc_color("white"), linesty=rsty
     endif
     ;; Label
     xlbl = 10.
     xyouts, xlbl, 380., lbl, CHARSIZE=1.5,CHARTHICK=2.5 , color=fsc_color("white") 
     xyouts, xlbl, 20., 'HI', CHARSIZE=1.5,CHARTHICK=2.5 , color=fsc_color("white") 
  endfor

  ;; Color bar
  if ~keyword_set(colpos) then colpos=[0.96,0.13,0.99,0.83]
  if ~keyword_set(nocolbar) then fsc_colorbar, range=scale,FORMAT=format, /vertical,$
     title=coltitle, POSITION=colpos, charsize=csz, annotatecolor=annotatecolor, _extra=extra
  xyouts, 26500., 18300, 'EW', color=fsc_color("white"), /device, charsi=1.5, align=0.5
  xyouts, 26500., 17700, '(Ang)', color=fsc_color("white"), /device, charsi=1.5, align=0.5
  IF keyword_set(PSFILE) THEN x_psclose




  ;;RESTORE COLOR BAR
  ctload, 39


END
