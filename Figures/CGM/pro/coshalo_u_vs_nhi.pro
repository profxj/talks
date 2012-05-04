pro coshalo_u_vs_nhi, grid, struct
;program to loop through all galaxies and make cloudy "solution" plots
;for metallicity and ionization parameter. 
;; jw_solvecloudy uses 2d krigging to get the interpolation of the
;; cloudy grid. KRIG2D is very very slow. Since it has to be called
;; multiple times, it makes running the program take approximately 1
;; -2 hours on the finest cloudy grid with a mesh grid of 100 X 100. 
;;If that time is unacceptable, try jw_solvecloudy faster!
;;; This uses a triangulation method which is supposed to be faster. 
;; It is way faster! Holy Shit! Use this version. 
;;; both jw_solvecloudy and jw_solvecloudyfaster should be run on a
;;; machine with power!! (i.e. a macbookpro will be painfully slow)
;; jw 2/4/12
; 
ldir = getenv('DROPBOX_DIR')+'/COS-Halos/lowions/'
if not keyword_set(grid) then grid = ldir+'finercldygrid_logn_coshalos.fits.gz'
if not keyword_set(struct) then struct = ldir+'cosmetals_megastructure.sav'
tdir = getenv('DROPBOX_DIR')+'/COS-Halos/Targets/'
hidir =  getenv('DROPBOX_DIR')+'/COS-Halos/Analysis/HI/'
restore, file = ldir+struct ; METASTRUCTURE, called megastruct
grid = mrdfits(ldir+grid, 1) ; the grid must be in the working directory. 

ngals = n_elements(megastruct.galaxy.galid) ; proper variable name
;;  Only do Mg if it's detected.
elementarr = [14,6,7,12,8] ;; Si, C, N, Mg, and Oxygen
nele = n_elements(elementarr)
nios = 3 ; max number of states per element. should be fine. 
 uniqhivals = uniq(grid.nhi, sort(grid.nhi))
 numplots = n_elements(uniqhivals)


;;; the fine spacing of the interpolated grid. 

meshx = 100L ; u
meshy = 100L ; z

;; the u and z values

ukrig = fltarr(meshx) + min(grid.u)+ (findgen(meshx) * (max(grid.u) - min(grid.u)) / (meshx - 1) ) 
zkrig =  fltarr(meshy) + min(grid.feh)+ (findgen(meshy) * (max(grid.feh) - min(grid.feh)) / (meshy - 1) ) 

;;; set the mask values to 1 by default. Will be 0 for unallowed
;;; areas. 
;; create the final mask array for storing the UZ solutions
mask_uz = bytarr(ngals,numplots,nele, meshx, meshy)
mask_uz[*] = 1
print, 'finished with mesh.'

;; the final output from the loop. allowable u and z ranges for each
;; HI level and for each element (in theory you wouldn't have
;; to do each element, but it's too restrictive otherwise.  

umax = fltarr(ngals,numplots,nele)
umin = fltarr(ngals,numplots,nele)
zmax = fltarr(ngals,numplots,nele)
zmin = fltarr(ngals,numplots,nele)



FOR i = 28, ngals - 1. do begin
if i eq 29 then begin
STOP
endif


psfil1 = tdir+megastruct[i].galaxy.qso+'/'+$
         megastruct[i].galaxy.sysnm+'/plots/'+megastruct[i].galaxy.sysnm+$
         '_fastcloudycolumns.ps'
print, 'Will write: '+psfil1

;;;;;;;;;;;;;;; Log U, Log Z contour plots

set_plot, 'ps'
device, filename= psfil1, /color,$
        /tt_font,set_font = 'Times',xsize=22.,ysize=18.
         device, decomposed = 0,  /ISOLATIN1, bits_per_pixel = 8
         !p.font = 0
         !p.thick = 3
         !x.thick = 3
         !y.thick = 3
         !p.charthick = 6
         !P.CHARSIZE = 2.0
         !p.multi = [0, 2, 2, 0, 0]
         

         clr = getcolor(/load)
         xclr = jw_setclrs(/white)
         
        
          xmin = -5.0
          xmax = 0.0
          ymin = -3.0
          ymax = 0.0

         ;;;;; Now get the measured HI for each system from
         ;;;;; Chris Thom's analysis. Some lower limits have
         ;;;;; upper limits based on the absence of Damping wings from
         ;;;;; their HI profiles. It seems that Chris has not
         ;;;;; accounted for this, so some may need to be revisited
         ;;;;; still. 
              hifile = mrdfits(hidir+'HI_results.fits', 1)
              ;;; hiflags are 1 = good, 4 = non detection, upper
              ;;; limit, 8 = lower limit (saturated) 
              ;; indices are probably the same, since they come from
              ;; same "systems to calculate" file. BUT, just to be
              ;; sure the indexing is correct. 

              calcgal = where(hifile.qso eq megastruct[i].galaxy.field AND strtrim(hifile.galname, 2) eq megastruct[i].galaxy.galid)

              lognvalhi =  hifile[calcgal].nhi
              uplognvalhi = lognvalhi + (1.0 * hifile[calcgal].nhierr)
              lowlognvalhi = lognvalhi - (1.0 * hifile[calcgal].nhierr )
              labuphi = strtrim(string(uplognvalhi, format = '(f5.1)'), 2)
              lablowhi = strtrim(string(lowlognvalhi,  format = '(f5.1)'), 2)
              labvalhi = strtrim(string(lognvalhi,  format = '(f5.1)'), 2)

;;;; Start with Si Page. 
;;;; Make 8 plots for different NHI values. 

     
      
      ;; Make the solution arrays storeable

     

      FOR ii = 0, numplots - 1. DO BEGIN

         print, 'HI plot: ', ii

         FOR zx = 0., nele - 1. DO BEGIN

            print, 'element: ', zx
      
            plot, [-10],[-10], xrange = [xmin, xmax], yrange = [ymin, ymax], $
                  xsty = 1, ysty = 1, xtitle = 'Log U', $
                  ytitle = 'Log [M/H]', chars = 1.2, xmarg = [4,4], $
                  ymarg = [2,2]

            xyouts, -4.75, -2.85, 'log N!dHI!n = '$
                    +strtrim(string(grid[uniqhivals[ii]].nhi, $
                                    format = '(f5.1)'), 2), $
                    chars = 1.0

            colindex = where(grid.nhi eq grid[uniqhivals[ii]].nhi)

            ;;;Get the Actual Data on there. 
            ; First for Si
            
           

            zoi = elementarr[zx]
            ioi = [1,2,3,4] ; 6 for OVI in specialized cases. See if this works. 
            
            nions = n_elements(ioi)

            strioi = [ 'I', 'II','III', 'IV']

            IF (i eq 26 OR i eq 29 or i eq 33 or i eq 35 or i eq 36 or i eq 38 or i eq 45) THEN BEGIN

               ioi = [1,2,3] ; SiIV is not consistent here.  
            
               nions = n_elements(ioi)

            strioi = [ 'I', 'II','III']

            ENDIF


            getabnd, elm, zoi, abnd, flag = 1

            xyouts, -4.75, -2.55, elm, chars = 1.0
            xyouts, -4.75, -2.65, megastruct[i].galaxy.qso + ' '+$
                    megastruct[i].galaxy.sysnm, chars = 0.8

            FOR zz = 0., n_elements(ioi) - 1. DO BEGIN

           
             ;; unusable = 0, uplim = 1, lowlim = 2, good = 3
            ; default is nothing

               ;;; The Values

              lognval =  megastruct[i].ion[zoi,ioi[zz]].lognion
              uplognval = lognval + (1.0 * megastruct[i].ion[zoi,ioi[zz]].siglognion)
              lowlognval = lognval - (1.0 * megastruct[i].ion[zoi,ioi[zz]].siglognion)
              labup = strtrim(string(uplognval, format = '(f5.1)'), 2)
              lablow = strtrim(string(lowlognval,  format = '(f5.1)'), 2)
              labval = strtrim(string(lognval,  format = '(f5.1)'), 2)

             

              nm = elm+strioi[zz]

               ;;;;;;

            IF (megastruct[i].ion[zoi, ioi[zz]].nflg EQ 0 OR (zoi eq 12 AND ioi[zz] eq 1)) THEN BEGIN
              print, 'No Coverage or Bad:', zoi, ioi[zz]
              megastruct[i].ion[zoi, ioi[zz]].nflg = 0
              ;; mask values are unchanged (all 1's) everything is allowed if
              ;; there are no observational constraints, after-all. 
              ;;; Mg I is too hard to model? It is not included on
              ;;; these plots since it is often quite discrepant. 

            ENDIF

          
             IF (megastruct[i].ion[zoi, ioi[zz]].nflg EQ 1) THEN BEGIN
              print, 'Good Detection.', zoi, ioi[zz]
             

              contour, grid[colindex].lognion[zoi,ioi[zz]], $
                       grid[colindex].u, grid[colindex].feh, $
                       /irreg, /over, c_label = [1, 1], color = xclr[zz + 1], $
                       level = [lowlognval, uplognval], $
                       c_annot = [nm+' = '+lablow, nm+' = '+labup], $
                       c_chars = 0.5, /follow
              ; find the allowed z and u values by triangulation

              sortu = sort(grid[colindex].u)
              triangulate, grid[colindex[sortu]].u, grid[colindex[sortu]].feh, $
                           trivert
              thistri = TRIGRID(grid[colindex[sortu]].u, grid[colindex[sortu]].feh, $
                                grid[colindex[sortu]].lognion[zoi,ioi[zz]], $
                                trivert, nx = meshx, ny = meshy)

            
              print, 'done triangulating for:',  zoi, ioi[zz]

              ;; now find where the interpolated surface is outside
              ;; measured values (i.e. not allowed) and set to 0s. 

              forbid = where(thistri lt lowlognval OR thistri gt uplognval)
              
              ; get the indices 

              IF (forbid[0] ge 0) THEN BEGIN
                 thistri[*] = 1 ;make thiskrig mirror the mask to save time. 
                  thistri[forbid] = 0 ; now set the forbidden values to 0
                  forbidu = forbid MOD meshx
                  forbidz = forbid/meshx
                  ; older 0's will multiply out
                  mask_uz[i,ii,zx,*, *] = mask_uz[i,ii,zx,*,*] * thistri  

              ENDIF

              

        

            
           ENDIF
             

               IF (megastruct[i].ion[zoi, ioi[zz]].nflg EQ 2) THEN BEGIN
              print, 'Lower Limit.', zoi, ioi[zz]
             

              contour, grid[colindex].lognion[zoi,ioi[zz]], $
                       grid[colindex].u, grid[colindex].feh, $
                       /irreg, /over, c_label = [1], color = xclr[zz + 1], $
                       level = [lognval], $
                       c_annot = [nm+' > '+labval], $
                       c_chars = 0.5, /follow

               ; find the allowed z and u values by triangulation

              sortu = sort(grid[colindex].u)
              triangulate, grid[colindex[sortu]].u, grid[colindex[sortu]].feh, $
                           trivert
              thistri = TRIGRID(grid[colindex[sortu]].u, grid[colindex[sortu]].feh, $
                                grid[colindex[sortu]].lognion[zoi,ioi[zz]], $
                                trivert, nx = meshx, ny = meshy)

            
              print, 'done triangulating for:',  zoi, ioi[zz]

              ;; now find where the interpolated surface is outside
              ;; measured values (i.e. not allowed) and set to 0s. 

              forbid = where(thistri lt lognval)
              
              ; get the indices 

              IF (forbid[0] ge 0) THEN BEGIN
                 thistri[*] = 1 ;make thiskrig mirror the mask to save time. 
                  thistri[forbid] = 0 ; now set the forbidden values to 0
                  forbidu = forbid MOD meshx
                  forbidz = forbid/meshx
                  ; older 0's will multiply out
                  mask_uz[i,ii,zx,*, *] = mask_uz[i,ii,zx,*,*] * thistri  

              ENDIF

     


           ENDIF

                IF (megastruct[i].ion[zoi, ioi[zz]].nflg EQ 3) THEN BEGIN
              print, 'Upper Limit.', zoi, ioi[zz]
             

              contour, grid[colindex].lognion[zoi,ioi[zz]], $
                       grid[colindex].u, grid[colindex].feh, $
                       /irreg, /over, c_label = [1], color = xclr[zz + 1], $
                       level = [lognval], $
                       c_annot = [nm+' < '+labval], $
                       c_chars = 0.5, /follow

               ; find the allowed z and u values by triangulation

              sortu = sort(grid[colindex].u)
              triangulate, grid[colindex[sortu]].u, grid[colindex[sortu]].feh, $
                           trivert
              thistri = TRIGRID(grid[colindex[sortu]].u, grid[colindex[sortu]].feh, $
                                grid[colindex[sortu]].lognion[zoi,ioi[zz]], $
                                trivert, nx = meshx, ny = meshy)

            
              print, 'done triangulating for:',  zoi, ioi[zz]

              ;; now find where the interpolated surface is outside
              ;; measured values (i.e. not allowed) and set to 0s. 

              forbid = where(thistri gt lognval)
              
              ; get the indices 

              IF (forbid[0] ge 0) THEN BEGIN
                 thistri[*] = 1 ;make thiskrig mirror the mask to save time. 
                  thistri[forbid] = 0 ; now set the forbidden values to 0
                  forbidu = forbid MOD meshx
                  forbidz = forbid/meshx
                  ; older 0's will multiply out
                  mask_uz[i,ii,zx,*, *] = mask_uz[i,ii,zx,*,*] * thistri  

              ENDIF

     
             

            ENDIF

          
              

             ENDFOR ; close ion for loop
          
            

        

       
         ;;; now the min and max values for every HI value for every element
         
         goodvals = where(mask_uz[i,ii,zx,*,*] eq 1) ; maybe there are some left...
        
        
         
         IF (goodvals[0] GE 0) THEN BEGIN
         
            print, 'we have a solution for this HI level! And this element'
            print, strtrim(string(grid[uniqhivals[ii]].nhi, format = '(f5.1)'), 2)

         ;first find u and z indexes

         goodu = goodvals MOD meshx
         goodz = goodvals/meshx

                                ; since they are the same size from
                                ; the grid values, take the max index
                                ; number and use the u value from
                                ; that. 
         maxgoodu = max(goodu)
         mingoodu = min(goodu)
         maxgoodz = max(goodz)
         mingoodz = min(goodz)
         
         ;; one for each HI level. And, for each galaxy. And each
         ;; element. 

         umax[i,ii,zx] = ukrig[maxgoodu]
         umin[i,ii,zx] = ukrig[mingoodu]
         ;;; Add a dex to each for non solar abundance ratios. 
         zmax[i,ii,zx] = zkrig[maxgoodz] + 0.1
         zmin[i,ii,zx] = zkrig[mingoodz] - 0.1
         
      ENDIF ELSE BEGIN
            print, 'This HI level is not allowed!!!!'
            print, strtrim(string(grid[uniqhivals[ii]].nhi, format = '(f5.1)'), 2)

            ;;; So, in some cases the mask will be all zeros... just
            ;;; set umax, umin, zmax, zmin equal to -999.9

            umax[i,ii,zx] = -999.9
            umin[i,ii,zx] = -999.9
            zmax[i,ii,zx] = -999.9
            zmin[i,ii,zx] = -999.9


      ENDELSE

      

   ENDFOR ; close element for loop      
   ENDFOR                       ; close HI for loop

      device, /close
      ; make into a pdf file. 
      spawn, 'idlepstopdf.sh '+psfil1

    ;;solve

;; minzmax and maxzmin for every hivalue. 
minzmax = fltarr(ngals, numplots)
maxzmin = fltarr(ngals, numplots)
minumax = fltarr(ngals, numplots)
maxumin = fltarr(ngals, numplots)
hisolz = fltarr(ngals, numplots)
hisolu = fltarr(ngals, numplots)

FOR iv = 0., numplots - 1. DO BEGIN

minzmax[i,iv] = min(zmax[i,iv,*], minzm)
maxzmin[i,iv] = max(zmin[i,iv,*], maxzm)
minumax[i,iv] = min(umax[i,iv,*], minum)
maxumin[i,iv] = max(umin[i,iv,*], maxum)

ENDFOR


;; now find the HI values where they cross. This is now mostly
;; irrelevant. 

FOR iv = 0., numplots - 1. DO BEGIN
IF (minzmax[i,iv] gt maxzmin[i,iv]) THEN BEGIN
hisolz[i,iv] = 1
ENDIF ELSE BEGIN
hisolz[i,iv] = 0
ENDELSE

IF (minumax[i,iv] gt maxumin[i,iv]) THEN BEGIN
hisolu[i,iv] = 1
ENDIF ELSE BEGIN
hisolu[i,iv] = 0
ENDELSE

ENDFOR

crossz = where(hisolz[i,*] eq 1)
IF (crossz[0] ge 0) THEN BEGIN
allowedhiz = grid[uniqhivals[crossz]].nhi
ENDIF ELSE BEGIN
;would be a little weird. 
allowedhiz = [15.0, 19.0]
ENDELSE

crossu = where(hisolu[i,*] eq 1)
IF (crossu[0] ge 0) THEN BEGIN
allowedhiu = grid[uniqhivals[crossu]].nhi
ENDIF ELSE BEGIN
allowedhiu = [15.0, 19.0]
ENDELSE
nallow = n_elements(allowedhiu)
allowedhi = fltarr(nallow)


FOR iv = 0., nallow - 1. DO BEGIN
allowedhi[iv] = where(allowedhiz eq allowedhiu[iv])
ENDFOR

goodhi = where(allowedhi ge 0)

;;;Get range of Allowed HI!Finally!!

;;;;;;;;;;;;;;;;;;OLDER VERSION
;IF (goodhi[0] GE 0) THEN BEGIN
;hisol = grid[uniqhivals[crossz[allowedhi[goodhi]]]].nhi

;maxhisol = max(hisol)
;minhisol = min(hisol)
;imaxhi = where(grid[uniqhivals].nhi eq maxhisol)
;iminhi = where(grid[uniqhivals].nhi eq minhisol)
;ENDIF ELSE BEGIN

;imaxhi = 40 ;; just make it 19.0. will be notable on plot. 
;iminhi = 15 ;; make it 16.5

;ENDELSE
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;INSTEAD USE MEASUREMENTS FOR SHOWCASED RANGE. 
IF (hifile[calcgal].nhiflg EQ 0 AND goodhi[0] GE 0) THEN BEGIN

print, 'Bad HI solution.'

hisol = grid[uniqhivals[crossz[allowedhi[goodhi]]]].nhi
maxhisol = max(hisol)
minhisol = min(hisol)
imaxhi = where(grid[uniqhivals].nhi eq maxhisol)
iminhi = where(grid[uniqhivals].nhi eq minhisol)
hiclr = xclr[2]
hiflg = 0
ENDIF ELSE BEGIN

imaxhi = 40 ;; just make it 19.0. will be notable on plot. 
iminhi = 15 ;; make it 16.5
hiflg = 0
ENDELSE

IF (hifile[calcgal].nhiflg EQ 1) THEN BEGIN
print, 'Good HI range'
hisol = lognvalhi
maxhisol = uplognvalhi
minhisol = lowlognvalhi
imaxhi = where(grid[uniqhivals].nhi ge maxhisol AND grid[uniqhivals].nhi lt maxhisol + 0.1)
iminhi = where(grid[uniqhivals].nhi le minhisol AND grid[uniqhivals].nhi gt minhisol - 0.1 )
IF (maxhisol gt 19.0) THEN BEGIN
imaxhi = 40
ENDIF
IF (maxhisol lt 15.0) THEN BEGIN
imaxhi = 0  ; sigh. this happens apparently. 

ENDIF
IF (minhisol lt 15.0) THEN BEGIN
iminhi = 0
ENDIF

IF (minhisol gt 19.0) THEN BEGIN
iminhi = 40
ENDIF

hiclr = xclr[3]
hiflg = 1
ENDIF

IF (hifile[calcgal].nhiflg EQ 8) THEN BEGIN ;; lower limit. most cases. 
print, 'Good HI range'
hisol = lognvalhi
maxhisol = 19.0
minhisol = lognvalhi
imaxhi = where(grid[uniqhivals].nhi ge maxhisol AND grid[uniqhivals].nhi lt maxhisol + 0.1)
iminhi = where(grid[uniqhivals].nhi le minhisol AND grid[uniqhivals].nhi gt minhisol - 0.1 )


IF(n_elements(iminhi) gt 1) THEN BEGIN
;;there is the one case of 18.00000 sigh. 
iminhi = max(iminhi)
ENDIF

IF (maxhisol gt 19.0) THEN BEGIN
imaxhi = 40
ENDIF
IF (maxhisol lt 15.0) THEN BEGIN
imaxhi = 0
ENDIF
IF (minhisol lt 15.0) THEN BEGIN
iminhi = 0
ENDIF
IF (minhisol gt 19.0) THEN BEGIN
iminhi = 40
ENDIF

;;;;;Added by hand for the hi by-eye where necessary. 

IF (i eq 0) THEN BEGIN
iminhi = 20 ;; want to see about 17. plotted here. 
ENDIF
IF (i eq 1) THEN BEGIN
iminhi = 0 ;; want to see about 17. plotted here. 
imaxhi = 30
ENDIF
IF (i eq 5) THEN BEGIN
iminhi = 20 ;; range from fit 
imaxhi = 30
ENDIF
IF (i eq 6) THEN BEGIN
iminhi = 15 ;; want to see about 17. plotted here. 
imaxhi = 40
ENDIF
IF (i eq 13) THEN BEGIN
iminhi = 0 ;;
imaxhi = 20 ;; strict max of 17. 
ENDIF
IF (i eq 15) THEN BEGIN
iminhi = 20 ;; want to see about 17.5 plotted here. 
imaxhi = 35
ENDIF
IF (i eq 16) THEN BEGIN
iminhi = 20 ;; want to see about 17.5 plotted here. 
imaxhi = 40
ENDIF
IF (i eq 17) THEN BEGIN
iminhi = 30 ;; want to see about 17.5 plotted here. 
imaxhi = 35
ENDIF
IF (i eq 18) THEN BEGIN
iminhi = 30 ;; want to see about 18.1 plotted here. 
imaxhi = 40
ENDIF
IF (i eq 19) THEN BEGIN
iminhi = 20 ;; 17. 
imaxhi = 35
ENDIF

IF (i eq 21) THEN BEGIN
iminhi = 14 ;; 16.4 
imaxhi = 25
ENDIF

IF (i eq 22) THEN BEGIN
iminhi = 30 ;; 18. 
imaxhi = 40
ENDIF

IF (i eq 23) THEN BEGIN
iminhi = 10 ;; 16. 
imaxhi = 30
ENDIF

IF (i eq 24) THEN BEGIN
iminhi = 10 ;; 16. 
imaxhi = 27
ENDIF

IF (i eq 26) THEN BEGIN
iminhi = 19 ;; 16.9 
imaxhi = 33
ENDIF

IF (i eq 28) THEN BEGIN
iminhi = 21 ;; 16.9 
imaxhi = 33
ENDIF

IF (i eq 31) THEN BEGIN
iminhi = 16 ;; 16.6 
imaxhi = 35  ;; 18.5
ENDIF

IF (i eq 32) THEN BEGIN
iminhi = 25 ;; 17.5 
imaxhi = 33  ;; 18.3
ENDIF

IF (i eq 33) THEN BEGIN
iminhi = 35 ;; 18.5 
imaxhi = 40  ;; 19
ENDIF

IF (i eq 35) THEN BEGIN
iminhi = 25 ;; 17.5 
imaxhi = 30  ;; 18.3
ENDIF

IF (i eq 36) THEN BEGIN
iminhi = 17 ;; 17.5 
imaxhi = 31  ;; 18.3
ENDIF

IF (i eq 39) THEN BEGIN
iminhi = 17 ;; 16.7 
imaxhi = 31  ;; 18.1
ENDIF

IF (i eq 42) THEN BEGIN
iminhi = 20 ;; 17 
imaxhi = 30  ;; 18
ENDIF

IF (i eq 45) THEN BEGIN
iminhi = 15 ;; 16.5 
imaxhi = 30  ;; 18
ENDIF

IF (i eq 46) THEN BEGIN
iminhi = 5 ;; 15.5 
imaxhi = 25  ;; 17.5
ENDIF

hiclr = xclr[1]
hiflg = 2
ENDIF

IF (hifile[calcgal].nhiflg EQ 4) THEN BEGIN ;; upper limit. (non detection)
 
print, 'Good HI range'
hisol = lognvalhi
maxhisol = lognvalhi
minhisol = 15.0
imaxhi = where(grid[uniqhivals].nhi ge maxhisol AND grid[uniqhivals].nhi lt maxhisol + 0.1)
iminhi = where(grid[uniqhivals].nhi le minhisol AND grid[uniqhivals].nhi gt minhisol - 0.1 )
IF (maxhisol gt 19.0) THEN BEGIN
imaxhi = 40
ENDIF
IF (maxhisol lt 15.0) THEN BEGIN
imaxhi = 0
ENDIF
IF (minhisol lt 15.0) THEN BEGIN
iminhi = 0
ENDIF
IF (minhisol gt 19.0) THEN BEGIN
iminhi = 40
ENDIF
hiclr = xclr[1]
hiflg = 3
ENDIF

;;;; print out a text file with some relevant information. 



  txfil1 = tdir+megastruct[i].galaxy.qso+'/'+$
         megastruct[i].galaxy.sysnm+'/plots/'+megastruct[i].galaxy.sysnm+$
         '_fastcloudysol.txt'
print, 'Will write: '+txfil1
openw, 3, txfil1, width = 200

       ;;; write a row for each column of HI. 
printf, 3, megastruct[i].galaxy.qso + ' ' + megastruct[i].galaxy.sysnm
printf, 3, 'Log NHI, Element, Umax, Umin, Zmax, Zmin'

FOR zx = 0., nele - 1. DO BEGIN

FOR ivv = 0, numplots - 1. DO BEGIN
lin = string(grid[uniqhivals[ivv]].nhi, format = '(f5.1)') + ' '
lin = lin + string(elementarr[zx], format = '(i2)') + ' ' 
lin = lin + string(umax[i,ivv,zx], format = '(f5.1)') + ' '
lin = lin + string(umin[i,ivv,zx], format = '(f5.1)') + ' '
lin = lin + string(zmax[i,ivv,zx], format = '(f5.1)') + ' '
lin = lin + string(zmin[i,ivv,zx], format = '(f5.1)') + ' '
printf, 3, lin

ENDFOR

ENDFOR

close, /all

      ;;; now make a nice solution plot for each galaxy. 
      psfil2 = tdir+megastruct[i].galaxy.qso+'/'+$
         megastruct[i].galaxy.sysnm+'/plots/'+megastruct[i].galaxy.sysnm+$
         '_fastcloudysol.ps'
print, 'Will write: '+psfil2

;;;;;;;;;;;;;;; Log U, Log Z contour plots

set_plot, 'ps'
device, filename= psfil2, /color,$
        /tt_font,set_font = 'Times',xsize=22.,ysize=18.
device, decomposed = 0,  /ISOLATIN1, bits_per_pixel = 8

         !p.font = 0
         !p.thick = 3
         !x.thick = 3
         !y.thick = 3
         !p.charthick = 6
         !P.CHARSIZE = 2.0
         !p.multi = [0, 2, 2]
         

         clr = getcolor(/load)
         xclr = jw_setclrs(/white)

        ;; First metallicity

          xmin = 15.0
          xmax = 19.0
          ymin = -3.0
          ymax = 0.5  ;; to put HI values on top. 
     
          plot, [-10],[-10], xrange = [xmin, xmax], yrange = [ymin, ymax], $
                  xsty = 1, ysty = 1, xtitle = 'Log N!dHI!n', $
                  ytitle = 'Log [M/H]', chars = 1.2, xmarg = [8,4], $
                  ymarg = [2,4]
          
           xyouts, 18.35, 0.79, megastruct[i].galaxy.qso + ' '+$
                    megastruct[i].galaxy.sysnm, chars = 1.2

           IF (hiflg eq 0) THEN BEGIN

              xyouts, 16.5, 0.2, 'No HI Measurement', chars = 1.2
              oplot, [15.0, 19.0], [0.2, 0.2], linestyle = 0, color = hiclr 

           ENDIF

            IF (hiflg eq 1) THEN BEGIN
              ; oplot, [lowlognvalhi, uplognvalhi], [0.2, 0.2],linestyle = 1, color = hiclr 
               plotsym, 1, 2
               oplot, [lognvalhi], [0.4], color = hiclr, psym = 8, thick = 5
                xyouts, uplognvalhi + 0.5, 0.3, 'HI Measurement = '+strtrim(string(lognvalhi, format = '(f5.1)'), 2) + ' +/- '+$
                      strtrim(string(hifile[calcgal].nhierr, format = '(f5.2)'), 2) , chars = 0.8
            ENDIF

             IF (hiflg eq 2) THEN BEGIN

                oplot, [lognvalhi, 18.8], [0.2,0.2], linestyle = 0,  color = hiclr 
                plotsym, 7, 2
                oplot, [18.6], [0.2], psym = 8, thick = 5, color = hiclr
               
                xyouts, 15.4, [0.3], 'HI Measured Lower Limit > '+strtrim(string(lognvalhi, format = '(f5.1)'), 2), chars = 0.8

             ENDIF

              IF (hiflg eq 3) THEN BEGIN
                 oplot, [15.0, lognvalhi], [0.2,0.2], linestyle = 0,  color = hiclr 
                plotsym, 6, 2
                oplot, [15.2], [0.2], psym = 8, thick = 5, color = hiclr
                
                xyouts, 16.0, [0.3], 'HI Measured Upper Limit < '+strtrim(string(lognvalhi, format = '(f5.1)'), 2), chars = 0.8
              ENDIF


           
          FOR zx = 0, nele - 1. DO BEGIN

          

            ;; Loop for lines connecting the points. 
            FOR iv = 0, numplots - 1. DO BEGIN
               FOR zz = 0, nions - 1. DO BEGIN
               IF (megastruct[i].ion[elementarr[zx], ioi[zz]].nflg ne 0) THEN BEGIN
            oplot,[grid[uniqhivals[iv]].nhi + (zx * 0.020), grid[uniqhivals[iv]].nhi + (zx * 0.020)], $
                   [zmax[i,iv,zx], zmin[i,iv,zx]], $
                   linestyle = 0, color = xclr[zx +1], thick = 4
         ENDIF
            ENDFOR



               IF (zmax[i,iv,zx] - zmin[i,iv,zx] lt 3.2) THEN BEGIN
            getabnd, elm, elementarr[zx], abnd, flag = 1
            xyouts, 18.65 + (zx * 0.5), 0.60, elm, chars = 1.0, color = xclr[zx +1], charthick = 5
         ENDIF


         ENDFOR
         ENDFOR


           ;;; now Ionization parameter

            ymin = -5.0
            ymax = 0.7
            plot, [-10],[-10], xrange = [xmin, xmax], yrange = [ymin, ymax], $
                  xsty = 1, ysty = 1, xtitle = 'Log  N!dHI!n', $
                  ytitle = 'Log U', chars = 1.2, xmarg = [4,8], $
                  ymarg = [2,4]  

             IF (hiflg eq 0) THEN BEGIN

              xyouts, 16.5, 0.3, 'No HI Measurement', chars = 1.2
              oplot, [15.0, 19.0], [0.2, 0.2], linestyle = 0, color = hiclr 

           ENDIF

            IF (hiflg eq 1) THEN BEGIN
              ; oplot, [lowlognvalhi, uplognvalhi], [0.2, 0.2],linestyle = 1, color = hiclr 
               plotsym, 1, 2
               oplot, [lognvalhi], [0.4], color = hiclr, psym = 8, thick = 5
               xyouts, uplognvalhi + 0.5, 0.3, 'HI Measurement = '+strtrim(string(lognvalhi, format = '(f5.1)'), 2) + ' +/- '+$
                      strtrim(string(hifile[calcgal].nhierr, format = '(f5.2)'), 2) , chars = 0.8
            ENDIF

             IF (hiflg eq 2) THEN BEGIN

                oplot, [lognvalhi, 18.8], [0.2,0.2], linestyle = 0,  color = hiclr 
                plotsym, 7, 2
                oplot, [18.6], [0.2], psym = 8, thick = 5, color = hiclr
                xyouts, 15.4, [0.3], 'HI Measured Lower Limit > '+strtrim(string(lognvalhi, format = '(f5.1)'), 2), chars = 0.8

             ENDIF

              IF (hiflg eq 3) THEN BEGIN
                 oplot, [15.0, lognvalhi], [0.2,0.2], linestyle = 0,  color = hiclr 
                plotsym, 6, 2
                oplot, [15.2], [0.2], psym = 8, thick = 5, color = hiclr
                xyouts, 16.0, [0.3], 'HI Measured Upper Limit <  '+strtrim(string(lognvalhi, format = '(f5.1)'), 2), chars = 0.8
              ENDIF



             FOR zx = 0, nele - 1. DO BEGIN

          ; oplot, grid[uniqhivals].nhi + (zx * 0.03), umax[i,*,zx], psym = sym(1), color = xclr[zx + 1], $
                                symsize = 0.5

           ; oplot, grid[uniqhivals].nhi + (zx * 0.03), umin[i,*,zx], psym = sym(1), color = xclr[zx + 1], $
                                symsize = 0.5

            ;; Loop for lines connecting the points. 
            FOR iv = 0, numplots - 1. DO BEGIN
                FOR zz = 0, nions - 1. DO BEGIN
               IF (megastruct[i].ion[elementarr[zx], ioi[zz]].nflg ne 0) THEN BEGIN
            oplot, [grid[uniqhivals[iv]].nhi + (zx * 0.020), grid[uniqhivals[iv]].nhi + (zx * 0.020)], $
                    [umax[i,iv,zx], umin[i,iv,zx]], $
                   linestyle = 0, color = xclr[zx +1], thick =4
         ENDIF
            ENDFOR


         ENDFOR
         ENDFOR

         
             xmin = -5.0
             xmax = 0.0
             ymin = -3.0
             ymax = 0.0

              plot, [-10],[-10], xrange = [xmin, xmax], yrange = [ymin, ymax], $
                  xsty = 1, ysty = 1, xtitle = 'Log U', $
                  ytitle = 'Log [M/H]', chars = 1.2, xmarg = [8,4], $
                  ymarg = [4,2]
             
                
            ;;; more nested for loops...

              

              
              FOR zx = 0, nele - 1. DO BEGIN
                 
                 ;; find the right values for the lower column density
                 ;; estimate. 
                    gridsize = ((meshy * meshx) - (meshx/10))
                    goodvals2 = where(mask_uz[i,iminhi,zx,*,*] eq 1)
                    ngood = n_elements(goodvals2)
                    IF (goodvals2[0] ge 0 AND ngood le gridsize) THEN BEGIN
                    goodu2 = goodvals2 MOD meshx
                    goodz2 = goodvals2/meshx
                    FOR iii=0., ngood - 1. DO BEGIN

                       ;;;;A poor woman's attempt at transparency. 
                       IF (gridsize GE 40000) THEN BEGIN
                       IF (iii MOD 3 eq 0) THEN BEGIN
                       oplot, [ukrig[goodu2[iii]]], [zkrig[goodz2[iii]]], psym = sym(1), $
                              symsize = (0.7 /((1.1 * zx) + 1)), color = xclr[zx + 1]
                    ENDIF

                    ENDIF ELSE BEGIN
                       oplot, [ukrig[goodu2[iii]]], [zkrig[goodz2[iii]]], psym = sym(1), $
                              symsize = (0.7 /((1.1 * zx) + 1)), color = xclr[zx + 1]
                    ENDELSE

                      ENDFOR

                ENDIF

              ENDFOR

              xyouts, -4.75, -2.5, 'Log N!dHI!n = '+string(grid[uniqhivals[iminhi]].nhi, format = '(f5.1)'), $
                      chars = 1.2, charthick = 3

               plot, [-10],[-10], xrange = [xmin, xmax], yrange = [ymin, ymax], $
                  xsty = 1, ysty = 1, xtitle = 'Log U', $
                  ytitle = 'Log [M/H]', chars = 1.2, xmarg = [4,8], $
                  ymarg = [4,2]
          
             FOR zx = 0, nele - 1. DO BEGIN
                 
                 ;; find the right values for the lower column density
                 ;; estimate. 

                    goodvals2 = where(mask_uz[i,imaxhi,zx,*,*] eq 1)
                    ngood = n_elements(goodvals2)
                    IF (goodvals2[0] ge 0 AND ngood le gridsize) THEN BEGIN
                    goodu2 = goodvals2 MOD meshx
                    goodz2 = goodvals2/meshx
                    FOR iii=0., ngood - 1. DO BEGIN
                     
                        ;;;;A poor woman's attempt at transparency. 
                       IF (gridsize GE 40000) THEN BEGIN
                       IF (iii MOD 3 eq 0) THEN BEGIN
                       oplot, [ukrig[goodu2[iii]]], [zkrig[goodz2[iii]]], psym = sym(1), $
                              symsize = (0.7 /((1.1 * zx) + 1)), color = xclr[zx + 1]
                    ENDIF

                    ENDIF ELSE BEGIN
                       oplot, [ukrig[goodu2[iii]]], [zkrig[goodz2[iii]]], psym = sym(1), $
                              symsize = (0.7 /((1.1 * zx)  + 1)), color = xclr[zx + 1]
                    ENDELSE

                      ENDFOR

                ENDIF

              ENDFOR

              xyouts, -4.75, -2.5, 'Log N!dHI!n = '+string(grid[uniqhivals[imaxhi]].nhi, format = '(f5.1)'), $
                      chars = 1.2, charthick = 3

          device, /close
      ; make into a pdf file. 
      spawn, 'idlepstopdf.sh '+psfil2


   ENDFOR                       ; close galaxy for loop






end
