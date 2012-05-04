pro coshalo_u_vs_nhi, struct, WHITE = white

;;; Program to take the solutions from jw_lowionsolution.txt
;;; And plot them up and correlate with galaxy properties, etc. !!

;; written 3/4/12 by JW. 
;; updated to make nicer looking plots...5/2/12

;; N.B. jw_lowionsolution.txt is written by hand, thus, 
;;; if anything changes in the solution, it will have to be put in by
;;; hand. Running the metaldriver will not be enough. 

if not keyword_set(struct) then struct = 'cosmetals_megastructure.sav'
ldir = getenv('DROPBOX_DIR')+'/COS-Halos/lowions/'

restore, file = ldir+struct ; METASTRUCTURE, called megastruct
;; read in the solution. HI limits are weird!

lowsol = ldir+'jw_lowionsolution.txt'
readcol, lowsol, field, id, nhi, nhierr, nhiflg, nhi_add_uplim, $
         nhi_add_lowlim, nhimulti, umin, umax, zmin,  zmax, multiphase, $
         goodsol, format = '(a, a, f, f, i, f, f, a, f, f, f, f, a, i)'

;;; goodsol ranges from 1 to 5. 1 is unusuable, essentially. 
;;; nhiflg is as follows: 1 = good, 4 = non-detection, 8 = lower limit
;;; in the above cases, we use the simple nhi value. 
;;; however, if the flag is 13, I have added another upper limit
;;; (nhi_add_uplim) -- and use the lower limit in nhi
;;; if the flag is 15, I have added another upper limit AND another
;;; lower limit (use both nhi_add_uplim and nhi_add_lowlim)

;;;; If the flag is 17, I have added an additional lower limit only,
;;;; and thus use the lower limit in nhi_add_lowlim

nh = -2.5 ;; this is what I use for cloudy

usenhilow = findgen(n_elements(id)) * 0.0 
usenhiup = findgen(n_elements(id)) * 0.0 
usenhival =  findgen(n_elements(id)) * 0.0 

hiclr = usenhilow + 1.0


hicase1 = where(nhiflg eq 1)

usenhilow[hicase1] = nhi[hicase1] - nhierr[hicase1]
usenhiup[hicase1] = nhi[hicase1] + nhierr[hicase1] 
usenhival[hicase1] = nhi[hicase1]
hiclr[hicase1] = 1.0

hicase2 = where(nhiflg eq 4)
usenhilow[hicase2] = -99.
usenhiup[hicase2] = nhi[hicase2]
usenhival[hicase2] = nhi[hicase2]
hiclr[hicase2] = 2.0

;hicase3 = where(nhiflg eq 8)
;usenhilow[hicase3] = nhi[hicase3]
;usenhiup[hicase3] = 19.
;usenhival[hicase3] = nhi[hicase3]
;hiclr[hicase3] = 3.0

hicase4 = where(nhiflg eq 13)
usenhilow[hicase4] = nhi[hicase4]
usenhiup[hicase4] = nhi_add_uplim[hicase4]
FOR i = 0., n_elements(hicase4) -1. do begin
usenhival[hicase4[i]] = alog10(mean([10.^(nhi[hicase4[i]]), 10.^(nhi_add_uplim[hicase4[i]])])) 
ENDFOR

hiclr[hicase4] = 4.0

hicase5 = where(nhiflg eq 17)
usenhilow[hicase5] = nhi_add_lowlim[hicase5]
usenhiup[hicase5] = 19.
usenhival[hicase5] = nhi_add_lowlim[hicase5]
hiclr[hicase5] = 5.0

hicase6 = where(nhiflg eq 15)

usenhilow[hicase6] = nhi_add_lowlim[hicase6]
usenhiup[hicase6] = nhi_add_uplim[hicase6]
FOR i = 0., n_elements(hicase6) -1. do begin
usenhival[hicase6[i]] = alog10(mean([10.^(nhi_add_lowlim[hicase6[i]]), 10.^(nhi_add_uplim[hicase6[i]])])) 
ENDFOR
hiclr[hicase6] = 6.0

;; that should cover it. There should be no zeros left. 
bad = where (usenhilow eq 0.0 OR usenhiup eq 0.0)
IF (bad[0] ge 0) THEN BEGIN

print, 'You did not calculate HI properly here'
STOP

ENDIF



;;; make lowest possible HI value = 15. This is for the cloudy stuff. 

othin = where(usenhilow lt 15.0)
usenhilow[othin] = 15.0
usenhiup[othin] = 15.0

;; make the highest value 10^19. Also for cloudy grid. 
othick = where (usenhilow gt 19.0)
usenhilow[othick] = 19.0
usenhiup[othick] = 19.0

;;; There are lots of comments in the lowionsolution.txt file, but we
;;; should be reading in only 51 valid lines here. 

;;; It should be matched to the order in the megastruct (by design),
;;; but of course we will make sure of this, and return an error if it
;;; is not. 

ngals = n_elements(megastruct.galaxy.galid) ; proper variable name
IF (n_elements(id) NE ngals) THEN BEGIN
print, 'Something is wrong with the solution txt. Not enough entries to match.'
STOP
ENDIF


;;;;;Start with the relations converting M* to mhalo, and a constant
;;;;;CGM mass. 

         mhalot = 10.^(findgen(101) / 100. * 5. + 9.) 

         mstart = mhalo2mstar(mhalot)

         nh1t = 1e19 


         
             
; first, plot the total of Mstar and M_CGM for a fixed CGM mass 
  m_cgm = (3e10/1e19*nh1t)+0.0*mhalot 
; from JXP calculation for NHI = 10^19 and R= 300 kpc fixed  

; now, derive the CGM mass based on the virial radius 
  rvir = mhalot * 0.0 
  ;for in=0, n_elements(rvir)-1 do rvir[in] = (nfw(mhalot[in], 15., 0.0)).r200 
  m_cgmv = 1.4 * nh1t * 1.67d-24 * (rvir * 3.086d21)^2 / 1.989d33 
;  oplot, alog10(mhalot), alog10((mstart+m_cgmv) / mhalot), psym=-3, color=xclr[2]
  
 

         gs = where(goodsol GE 4)
         ngs = n_elements(gs)
  


      goodz = where((zmax[gs] - zmin[gs]) le 1.5)  
      ngoodz = n_elements(goodz)


    us = where(megastruct[gs[goodz]].galaxy.sfr_uplim eq 'yes')
    sf =  where(megastruct[gs[goodz]].galaxy.sfr_uplim eq 'no')

 
      goodz = where((zmax[gs] - zmin[gs]) le 1.5 AND usenhilow[gs] lt 18.5)  
      ngoodz = n_elements(goodz)

      
     


    us = where(megastruct[gs[goodz]].galaxy.sfr_uplim eq 'yes')
    sf =  where(megastruct[gs[goodz]].galaxy.sfr_uplim eq 'no')

  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 

;;;;;;;;;;;;;;;

      goodz = where((umax[gs] - umin[gs]) le 1.5)  
      ngoodz = n_elements(goodz)

      impact = megastruct[gs].rhofinal
     
 

    us = where(megastruct[gs[goodz]].galaxy.sfr_uplim eq 'yes')
    sf =  where(megastruct[gs[goodz]].galaxy.sfr_uplim eq 'no')

         ;;; define some variables. 

         ;mbaryerr = m_cgmmean
         ;dumerr = m_cgmmean * 0.00 + 0.20
         zmed = findgen(ngs) * 0.00 
         umed = zmed
FOR i = 0., ngs - 1. do begin

  
   ;;;; the error for the baryonic mass. To make this easier, use a
   ;;;; mean. you may want to use different sizes for lower and upper
   ;;;; in the end. Probably you do. 

   ;mbaryerrup = abs(alog10(m_cgmmean[i]) - alog10(m_cgm300up[i]))
   ;mbaryerrlow = abs(alog10(m_cgmmean[i]) - alog10(m_cgm300low[i]))
   ;mbaryerr[i] = mean([mbaryerrup, mbaryerrlow])
   ;m_cgmmean[i] = alog10(m_cgmmean[i])
  ; zmed[i] = alog10(mean([10.^(zmin[gs[i]]), 10.^(zmax[gs[i]])]))
  ; umed[i] = alog10(mean([10.^(umin[gs[i]]), 10.^(umax[gs[i]])]))

   zmed[i] = (mean([(zmin[gs[i]]), (zmax[gs[i]])]))
   umed[i] = (mean([(umin[gs[i]]), (umax[gs[i]])]))
   
ENDFOR
  
;;;;;;;;;;;;;;;

;;;;;;Now, U vs. NHI



psfil21 = 'goodsol_hiu.eps'
print, 'Will write: '+psfil21

	x_psopen, psfil21, /maxs
;ct_load_colours
 clr = getcolor(/load)

;;Set up colour indices. These correspond to colours defined by jw_setclrs
   xclr = jw_setclrs(/dark)
   bluecol = clr.cyan
   redcol = clr.tomato
   fgcol = 255
   bgcol = 0
if ~keyword_set(white) then begin
ct_psfill_black, col = bgcol
endif


         ymin = -4.0
         ymax = -1.0
         xmin = 14
         xmax = 20.5


         yquant = usenhival[gs[goodz]]
         ylow = usenhilow[gs[goodz]]
         yup = usenhiup[gs[goodz]]

         plot, [-10],[-10], xrange = [xmin, xmax], yrange = [ymin, ymax], $
               xmarg=[8,2], ymarg=[4,1], $
               xsty = 1, ysty = 1, ytitle = 'log U', $
               xtitle = 'log N!dHI!n', chars = 2.2, /noerase, color = fgcol

    
       
FOR i = 0., ngoodz - 1. do begin

  

   oplot, [ylow[i], yup[i]],  $
          [umed[goodz[i]], umed[goodz[i]]], $
          linestyle = 0, $
          thick = 6, color = clr.darkgray

   oplot, [yquant[i], yquant[i]], [umin[gs[goodz[i]]], umax[gs[goodz[i]]]], $
          linestyle = 0, $
          thick = 6, color = clr.darkgray
  

  
   ; oplot, [alog10(m_cgm300low[goodz[i]]), alog10(m_cgm300up[goodz[i]])], $
   ;       [zmed, zmed],  linestyle = 0, $
    ;      color = xclr[hiclr[i]], thick = 4

ENDFOR

  us = where(megastruct[gs[goodz]].galaxy.sfr_uplim eq 'yes')
    sf =  where(megastruct[gs[goodz]].galaxy.sfr_uplim eq 'no')
;;;Make the symbol size range from 0.5 to 6, and assign a number based
;;;on the range plotted. 
    rangeu = abs(umin[gs[goodz]] - umax[gs[goodz]])
    rangez = abs(yup - ylow)
    rangeavg = ((rangeu + rangez) / 2.) + 0.1
    ssclosed = 1./(rangeavg * 0.77)
    ssopen = ssclosed + 0.5
nus = n_elements(us)
nsf = n_elements(sf)
FOR ii =0, nus -1 do begin
    oplot, [yquant[us[ii]]], [umed[goodz[us[ii]]]], psym = symcat(14, col = fgcol), syms = ssopen[us[ii]]
     oplot, [yquant[us[ii]]], [umed[goodz[us[ii]]]], psym = symcat(14, col = redcol), syms = ssclosed[us[ii]]
  endfor

For ii=0, nsf -1 do begin
    oplot, [yquant[sf[ii]]], [umed[goodz[sf[ii]]]], psym = symcat(15, col = fgcol), syms = ssopen[sf[ii]]
    oplot, [yquant[sf[ii]]], [umed[goodz[sf[ii]]]], psym = symcat(15, col = bluecol), syms = ssclosed[sf[ii]]

 endfor

   

 x_psclose

      ; make into a pdf file. 
      spawn, 'idlepstopdf.sh '+psfil21


;;;;;;;\

end
