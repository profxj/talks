; program to plot metalicities and redshifts
; The binned points are like Prochaska et al. 2003b
; Error bars are from bootstrap method

function logfunc, X, A
  return,[10^(A[1]*X +A[0])]
end

; /old sets the old binning, using only DLA's before the new ESI survey
; nofit plots things without the lines fitting it

; Bias introduces a reduction of the lowest bin by 0.1 dex based on
; nestor 2008 result. (Since lowest redshift guys are biased). 

; run does a running bin rather than 

pro fig_dla_mtlvsz, print=print, oldbin=oldbin, nofit=nofit, check=check, ihiz=ihiz, $
                    bias=bias, norman=norman, run=run, xavier=xavier, noterdaeme=noterdaeme

;parse_dlalst, ndla, '~/DLA/Lists/metal_MAR_new_lit.lst', /noelm, ROOT='~/DLA/'

if keyword_set(oldbin) then begin
parse_dlalst, sdla, '~/DLA/Lists/metal_MAR_prochaska03.lst', /noelm, ROOT='~/DLA/'
endif else begin
parse_dlalst, sdla, '~/DLA/Lists/metal_MAR.lst', /noelm, ROOT='~/DLA/'
endelse

parse_dlalst, adla, '~/DLA/Lists/metal_MAR_plus_hiz.lst', /noelm, ROOT='~/DLA/'

;parse_dlalst, sdla, '~/DLA/Lists/metal_MAR_all.lst', /noelm, ROOT='~/DLA/'

;parse_dlalst, sdla, '~/DLA/Lists/metal704_DLA.lst', /noelm, ROOT='~/DLA/'
;parse_dlalst, sdla, '~/DLA/Lists/all_mtl.lst', /noelm, ROOT='~/DLA/'
;parse_dlalst, sdla, '~/DLA/Lists/metal_dla.lst', /noelm


; Now lets read in my new hiz DLAs
;parse_dlalst, hdla, '~/DLA/Lists/metal_MAR_hiz_proposal.lst', /noelm, ROOT='~/DLA/'
parse_dlalst, hdla, '~/DLA/Lists/metal_MAR_hiz.lst', /noelm, ROOT='~/DLA/'

; Now lets read in the hires results
;parse_dlalst, rdla, '~/DLA/Lists/metal_MAR_HIRES_proposal.lst', /noelm, ROOT='~/DLA/'
parse_dlalst, rdla, '~/DLA/Lists/metal_MAR_HIRES.lst', /noelm, ROOT='~/DLA/'


 if keyword_set(bias) then begin
    print, 'CORRECTION FOR BIAS IN LOWEST REDSHIFT BIN!'
                              ; This is for the bias in lowest
                              ; redshift bin as described in nestor
                              ; 2008
    abiasmatch = where(adla.zabs le 1.5)
    sbiasmatch = where(sdla.zabs le 1.5)
        
    adla[abiasmatch].mtl = adla[abiasmatch].mtl-0.1
    sdla[sbiasmatch].mtl = sdla[sbiasmatch].mtl-0.1

 endif


;;;;;;; quick look at dispersion

 b = where(adla.flgmtl eq 1 or adla.flgmtl eq 4 or adla.flgmtl eq 14)
ddla = adla[b]
z2 = where(ddla.zabs ge 1.7 and ddla.zabs lt 3.0)
z3 = where(ddla.zabs ge 3.0 and ddla.zabs lt 4.0)
z4 = where(ddla.zabs ge 4.0)

;print, stddev(ddla[z2].mtl), stddev(ddla[z3].mtl),stddev(ddla[z4].mtl)

if not keyword_set(run) then begin
if keyword_set(oldbin) then begin
   zarr = intarr(6)
   bin = dblarr(6)
   binl = dblarr(6)
   binr = dblarr(6)
   
   
; Lets bin the data
; We want to make the bins have equal numbers of data points.
; First bin is everything less than z=1.5
; The rest between 1.5-4.5 we do by equal number of points 
   
   rest = where(sdla.zabs gt 1.5 and sdla.zabs le 4.5, num)
   nzabs = sdla[rest].zabs
   
   snzabs = nzabs[sort(nzabs)]
   numperbin = num/5
   
   bin1_1 = 0
   bin1_2 = 1.5
   bin2_1 = snzabs[0]
   bin2_2 = snzabs[numperbin]
   bin3_1 = snzabs[1+numperbin]
   bin3_2 = snzabs[1+(numperbin*2)]
   bin4_1 = snzabs[2+(numperbin*2)]
   bin4_2 = snzabs[2+ (numperbin*3)]
   bin5_1 = snzabs[3+(numperbin*3)]
   bin5_2 = snzabs[3+ (numperbin*4)]
   bin6_1 = snzabs[4+(numperbin*4)]
   bin6_2 = snzabs[num-1]
   zarr1 = where(sdla.zabs ge bin1_1 and sdla.zabs le bin1_2, n1)
   zarr2 = where(sdla.zabs ge bin2_1 and sdla.zabs le bin2_2, n2)
   zarr3 = where(sdla.zabs ge bin3_1 and sdla.zabs le bin3_2, n3)
   zarr4 = where(sdla.zabs ge bin4_1 and sdla.zabs le bin4_2, n4)
   zarr5 = where(sdla.zabs ge bin5_1 and sdla.zabs le bin5_2, n5)
   zarr6 = where(sdla.zabs ge bin6_1 and sdla.zabs le bin6_2, n6)


   binl[0] = min(sdla.zabs)
   binr[0] = 1.5
   binl[1] = snzabs[0]
   binr[1] = snzabs[numperbin]
   binl[2] = snzabs[1+numperbin]
   binr[2] = snzabs[1+(numperbin*2)]
   binl[3] = snzabs[2+(numperbin*2)]
   binr[3] = snzabs[2+(numperbin*3)]
   binl[4] = snzabs[3+(numperbin*3)]
   binr[4] = snzabs[3+(numperbin*4)]
   binl[5] = snzabs[4+(numperbin*4)]
   binr[5] = snzabs[num-1]
   
   zarr1 = where(sdla.zabs ge binl[0] and sdla.zabs le binr[0], n1)
   zarr2 = where(sdla.zabs ge binl[1] and sdla.zabs le binr[1], n2)
   zarr3 = where(sdla.zabs ge binl[2] and sdla.zabs le binr[2], n3)
   zarr4 = where(sdla.zabs ge binl[3] and sdla.zabs le binr[3], n4)
   zarr5 = where(sdla.zabs ge binl[4] and sdla.zabs le binr[4], n5)
   zarr6 = where(sdla.zabs ge binl[5] and sdla.zabs le binr[5], n6)

   bin[0] = mean(sdla[zarr1].zabs)
   bin[1] = mean(sdla[zarr2].zabs)
   bin[2] = mean(sdla[zarr3].zabs)
   bin[3] = mean(sdla[zarr4].zabs)
   bin[4] = mean(sdla[zarr5].zabs)
   bin[5] = mean(sdla[zarr6].zabs)

endif else begin
   



   zarr = intarr(9)
   bin = dblarr(9)
   binl = dblarr(9)
   binr = dblarr(9)
   
   
   rest = where(adla.zabs gt 1.5, num)
   nzabs = adla[rest].zabs
   snzabs = nzabs[sort(nzabs)]
   
   numperbin = num/8
   
   binl[0] = min(adla.zabs)
   binr[0] = 1.5
   binl[1] = snzabs[0]
   binr[1] = snzabs[numperbin]
   binl[2] = snzabs[1+numperbin]
   binr[2] = snzabs[1+(numperbin*2)]
   binl[3] = snzabs[2+(numperbin*2)]
   binr[3] = snzabs[2+(numperbin*3)]
   binl[4] = snzabs[3+(numperbin*3)]
   binr[4] = snzabs[3+(numperbin*4)]
   binl[5] = snzabs[4+(numperbin*4)]
   binr[5] = snzabs[4+(numperbin*5)]
   binl[6] = snzabs[5+(numperbin*5)]
   binr[6] = snzabs[5+(numperbin*6)]
   binl[7] = snzabs[6+(numperbin*6)]
   binr[7] = snzabs[6+(numperbin*7)]
   binl[8] = snzabs[7+(numperbin*7)]
   binr[8] = snzabs[num-1]
   
   
   zarr1 = where(adla.zabs ge binl[0] and adla.zabs le binr[0], n1)
   zarr2 = where(adla.zabs ge binl[1] and adla.zabs le binr[1], n2)
   zarr3 = where(adla.zabs ge binl[2] and adla.zabs le binr[2], n3)
   zarr4 = where(adla.zabs ge binl[3] and adla.zabs le binr[3], n4)
   zarr5 = where(adla.zabs ge binl[4] and adla.zabs le binr[4], n5)
   zarr6 = where(adla.zabs ge binl[5] and adla.zabs le binr[5], n6)
   zarr7 = where(adla.zabs ge binl[6] and adla.zabs le binr[6], n7)
   zarr8 = where(adla.zabs ge binl[7] and adla.zabs le binr[7], n8)
   zarr9 = where(adla.zabs ge binl[8] and adla.zabs le binr[8], n9)

   bin[0] = mean(adla[zarr1].zabs)
   bin[1] = mean(adla[zarr2].zabs)
   bin[2] = mean(adla[zarr3].zabs)
   bin[3] = mean(adla[zarr4].zabs)
   bin[4] = mean(adla[zarr5].zabs)
   bin[5] = mean(adla[zarr6].zabs)
   bin[6] = mean(adla[zarr7].zabs)
   bin[7] = mean(adla[zarr8].zabs)
   bin[8] = mean(adla[zarr9].zabs)


;  zarr = intarr(8)
;  bin = dblarr(8)
;  binl = dblarr(8)
;  binr = dblarr(8)
;  
;  
;  rest = where(adla.zabs gt 1.5, num)
;  nzabs = adla[rest].zabs
;  snzabs = nzabs[sort(nzabs)]
;  
;  numperbin = num/7
;  
;  binl[0] = min(adla.zabs)
;  binr[0] = 1.5
;  binl[1] = snzabs[0]
;  binr[1] = snzabs[numperbin]
;  binl[2] = snzabs[1+numperbin]
;  binr[2] = snzabs[1+(numperbin*2)]
;  binl[3] = snzabs[2+(numperbin*2)]
;  binr[3] = snzabs[2+(numperbin*3)]
;  binl[4] = snzabs[3+(numperbin*3)]
;  binr[4] = snzabs[3+(numperbin*4)]
;  binl[5] = snzabs[4+(numperbin*4)]
;  binr[5] = snzabs[4+(numperbin*5)]
;  binl[6] = snzabs[5+(numperbin*5)]
;  binr[6] = snzabs[5+(numperbin*6)]
;  binl[7] = snzabs[6+(numperbin*6)]
;  binr[7] = snzabs[num-1]
;  
;  
;  
;  zarr1 = where(adla.zabs ge binl[0] and adla.zabs le binr[0], n1)
;  zarr2 = where(adla.zabs ge binl[1] and adla.zabs le binr[1], n2)
;  zarr3 = where(adla.zabs ge binl[2] and adla.zabs le binr[2], n3)
;  zarr4 = where(adla.zabs ge binl[3] and adla.zabs le binr[3], n4)
;  zarr5 = where(adla.zabs ge binl[4] and adla.zabs le binr[4], n5)
;  zarr6 = where(adla.zabs ge binl[5] and adla.zabs le binr[5], n6)
;  zarr7 = where(adla.zabs ge binl[6] and adla.zabs le binr[6], n7)
;  zarr8 = where(adla.zabs ge binl[7] and adla.zabs le binr[7], n8)
;
;  bin[0] = mean(adla[zarr1].zabs)
;  bin[1] = mean(adla[zarr2].zabs)
;  bin[2] = mean(adla[zarr3].zabs)
;  bin[3] = mean(adla[zarr4].zabs)
;  bin[4] = mean(adla[zarr5].zabs)
;  bin[5] = mean(adla[zarr6].zabs)
;  bin[6] = mean(adla[zarr7].zabs)
;  bin[7] = mean(adla[zarr8].zabs)


;   zarr = intarr(7)
;   bin = dblarr(7)
;   binl = dblarr(7)
;   binr = dblarr(7)
;   
;   
;   rest = where(adla.zabs gt 1.5, num)
;   nzabs = adla[rest].zabs
;   snzabs = nzabs[sort(nzabs)]
;   
;   numperbin = num/6
;   
;   binl[0] = min(adla.zabs)
;   binr[0] = 1.5
;   binl[1] = snzabs[0]
;   binr[1] = snzabs[numperbin]
;   binl[2] = snzabs[1+numperbin]
;   binr[2] = snzabs[1+(numperbin*2)]
;   binl[3] = snzabs[2+(numperbin*2)]
;   binr[3] = snzabs[2+(numperbin*3)]
;   binl[4] = snzabs[3+(numperbin*3)]
;   binr[4] = snzabs[3+(numperbin*4)]
;   binl[5] = snzabs[4+(numperbin*4)]
;   binr[5] = snzabs[4+(numperbin*5)]
;   binl[6] = snzabs[5+(numperbin*5)]
;   binr[6] = snzabs[num-1]
;   
;   
;   
;   zarr1 = where(adla.zabs ge binl[0] and adla.zabs le binr[0], n1)
;   zarr2 = where(adla.zabs ge binl[1] and adla.zabs le binr[1], n2)
;   zarr3 = where(adla.zabs ge binl[2] and adla.zabs le binr[2], n3)
;   zarr4 = where(adla.zabs ge binl[3] and adla.zabs le binr[3], n4)
;   zarr5 = where(adla.zabs ge binl[4] and adla.zabs le binr[4], n5)
;   zarr6 = where(adla.zabs ge binl[5] and adla.zabs le binr[5], n6)
;   zarr7 = where(adla.zabs ge binl[6] and adla.zabs le binr[6], n7)
;
;   bin[0] = mean(adla[zarr1].zabs)
;   bin[1] = mean(adla[zarr2].zabs)
;   bin[2] = mean(adla[zarr3].zabs)
;   bin[3] = mean(adla[zarr4].zabs)
;   bin[4] = mean(adla[zarr5].zabs)
;   bin[5] = mean(adla[zarr6].zabs)
;   bin[6] = mean(adla[zarr7].zabs)


endelse

endif else begin

   ; this is for a running bin

; number per bin
nbin = 40

; number of points resulting from running bin
abin = n_elements(adla.zabs)-nbin+1

zabs = adla[sort(adla.zabs)].zabs

binl = dblarr(abin)
binr = dblarr(abin)
bin = dblarr(abin)

for idx=0, abin-1 do begin

; subarrays of zabs
szabs = zabs[idx:idx+nbin-1]

bin[idx] = mean(szabs)
binl[idx] = min(szabs)
binr[idx] = max(szabs)

endfor

endelse


; Do the bins from Xavier's NHI paper 2009
; This overrides anything else specified before. 
if keyword_set(xavier) then begin

   binl=[2.2,2.4,2.7,3.0,3.5,4.0]
   binr=[2.4,2.7,3.0,3.5,4.0,5.5]
   binr = binr + 0.0000001

   zarr1 = where(adla.zabs ge binl[0] and adla.zabs le binr[0], n1)
   zarr2 = where(adla.zabs ge binl[1] and adla.zabs le binr[1], n2)
   zarr3 = where(adla.zabs ge binl[2] and adla.zabs le binr[2], n3)
   zarr4 = where(adla.zabs ge binl[3] and adla.zabs le binr[3], n4)
   zarr5 = where(adla.zabs ge binl[4] and adla.zabs le binr[4], n5)
   zarr6 = where(adla.zabs ge binl[5] and adla.zabs le binr[5], n6)

   bin = dblarr(6)

   bin[0] = mean(adla[zarr1].zabs)
   bin[1] = mean(adla[zarr2].zabs)
   bin[2] = mean(adla[zarr3].zabs)
   bin[3] = mean(adla[zarr4].zabs)
   bin[4] = mean(adla[zarr5].zabs)
   bin[5] = mean(adla[zarr6].zabs)

endif


if keyword_set(noterdaeme) then begin

   binl=[1.7,2.23,2.60,2.88,3.20]
   binr=[2.2301,2.6001,2.8801,3.2001,5.19]
   binr = binr + 0.0000001

   zarr1 = where(adla.zabs ge binl[0] and adla.zabs le binr[0], n1)
   zarr2 = where(adla.zabs ge binl[1] and adla.zabs le binr[1], n2)
   zarr3 = where(adla.zabs ge binl[2] and adla.zabs le binr[2], n3)
   zarr4 = where(adla.zabs ge binl[3] and adla.zabs le binr[3], n4)
   zarr5 = where(adla.zabs ge binl[4] and adla.zabs le binr[4], n5)

   bin = dblarr(5)

   bin[0] = mean(adla[zarr1].zabs)
   bin[1] = mean(adla[zarr2].zabs)
   bin[2] = mean(adla[zarr3].zabs)
   bin[3] = mean(adla[zarr4].zabs)
   bin[4] = mean(adla[zarr5].zabs)

endif


lobin = bin-binl
hibin = binr-bin


mtl = dblarr(n_elements(bin))
mtlmedian = dblarr(n_elements(bin))
mtlavg = dblarr(n_elements(bin))
mtlerr = dblarr(n_elements(bin))

if keyword_set(oldbin) then begin
   for idx=0, n_elements(bin)-1 do begin
      
      samp = where(sdla.zabs ge binl[idx] and sdla.zabs le binr[idx], nmatch)
      
                                ; metallicity weighted by the column density
                                ; This is from page L10 in Prochaska 2003b. 
      mtl[idx] = total( (10.^sdla[samp].mtl)*10.^sdla[samp].NHI)/total(10.^sdla[samp].NHI)
      
      if nmatch gt 1 then begin
         mtlerr[idx] = bootstrap(10.^sdla[samp].mtl, weights=10.^sdla[samp].NHI, /umean, check=check)
      endif else mtlerr[idx]=0   
   endfor
endif else begin
   
   for idx=0, n_elements(bin)-1 do begin
      
      samp = where(adla.zabs ge binl[idx] and adla.zabs le binr[idx], nmatch)
      
                                ; metallicity weighted by the column density
                                ; This is from page L10 in Prochaska 2003b. 
      mtl[idx] = total( (10.^adla[samp].mtl)*10.^adla[samp].NHI)/total(10.^adla[samp].NHI)
      mtlmedian[idx] = median(adla[samp].mtl)
      mtlavg[idx] = mean(adla[samp].mtl)
      if nmatch gt 1 then begin
         mtlerr[idx] = bootstrap(10.^adla[samp].mtl, weights=10.^adla[samp].NHI, /umean, check=check)
         
      endif else mtlerr[idx]=0   
   endfor
endelse

logmtl = alog10(mtl)



  ; if keyword_set(bias) then begin
  ;    print, 'CORRECTION FOR BIAS IN LOWEST REDSHIFT BIN!'
  ;                              ; This is for the bias in lowest
  ;                              ; redshift bin as described in nestor
  ;                              ; 2008
  ;    logmtl[0] = logmtl[0] - 0.1
  ;    mtl[0] = 10^logmtl[0]
  ;
  ; endif


;; Doing 95 %, so 2 sigma, so 2 times the error
;logmtl_lo = logmtl-alog10(mtl-2*mtlerr)
;logmtl_hi = alog10(mtl+2*mtlerr)-logmtl

; Just 1 sigma now!
logmtl_lo = logmtl-alog10(mtl-mtlerr)
logmtl_hi = alog10(mtl+mtlerr)-logmtl

bformat='(F5.3,2x,F5.3,2x,F5.3,2x,F5.3,2x,F5.3)'
forprint, bin, binl, binr, mtl, mtlerr, textout='cosmetal.txt', /silent, /nocomment, format=bformat

;;logmtl_lo = logmtl-alog10(mtl-mtlerr)
;;logmtl_hi = alog10(mtl+mtlerr)-logmtl
;;;;;; logmtl_err = (2*mtlerr/mtl) ;;; linear approximation



; not right way to do error, but test it out
;baderr = ( (logmtl-alog10(mtl-mtlerr))+(alog10(mtl+mtlerr)-logmtl) )/2d
;fit = linfit(bin, alog10(mtl), measure_errors=baderr, sigma=sigma)

yerr = sqrt((mtlerr/mtl)^2.)/alog(10)

fit = linfit(bin, alog10(mtl), measure_errors=yerr, sigma=sigma)
;print, fit, sigma

; Lets do this properly, and see if we get similar answer. 
; Do a logorithmic fit on the non log data. 

fit2 = mpfitfun('logfunc', bin, mtl, mtlerr, [-0.6, -0.2], perr=perr, /quiet)

print, fit2, perr

;;; Now lets do the same, but this time, lets NOT include the lowest
;;; redshift bin
num = n_elements(bin)
ebin = bin[1:num-1]
emtl = mtl[1:num-1]
emtlerr = mtlerr[1:num-1]

fit3 = mpfitfun('logfunc', ebin, emtl, emtlerr, [-0.6, -0.2], perr=perr3, /quiet)

print, fit3, perr3

; lets fit to the metallicity weighted by the NHi rather than to the bins
fit4 = mpfitfun('logfunc', adla.zabs, 10^adla.mtl, 10^adla.sigmtl, [-0.6, -0.2], perr=perr4, weights=((10^adla.nhi)/2d20), /quiet)
;fit4 = mpfitfun('logfunc', adla.zabs, adla.mtl, [-0.6, -0.2], perr=perr4, /quiet, weights=adla.nhi)

print, fit4, perr4


; Lets bootstrap to get the error
;;;;;;;;;; ERROR ARRAY IS NOT RIGHT - LOG FUNCTION. However, ignored
;;;;;;;;;;                            with weights!!
;result = bootstrap(adla.zabs, 10^adla.mtl, weights=((10^adla.nhi)/2d20), check=check)

;print, result


;corr = r_correlate(adla.zabs, adla.mtl, /kendall)
;corr2 = r_correlate(bin,mtl, /kendall)


; lets take a look at the highest metallicities
zmax=4.7
hiz = where(adla.zabs ge zmax, num2)

mtlhiz = total( (10.^adla[hiz].mtl)*10.^adla[hiz].NHI)/total(10.^adla[hiz].NHI)
mtlhizerr = bootstrap(10.^adla[hiz].mtl, weights=10.^adla[hiz].NHI, /umean, check=check)
hizbin = mean(adla[hiz].zabs)
hizbin_lo = hizbin-min(adla[hiz].zabs)
hizbin_hi = max(adla[hiz].zabs)-hizbin
logmtlhiz = alog10(mtlhiz)
logmtlhiz_lo = logmtlhiz-alog10(mtlhiz-mtlhizerr)
logmtlhiz_hi = alog10(mtlhiz+mtlhizerr)-logmtlhiz

print, logmtlhiz, logmtlhiz_lo, logmtlhiz_hi
;forprint, adla[hiz].qso, adla[hiz].zabs, adla[hiz].mtl, adla[hiz].flgmtl


plotsym, 0, 0.8, /fill

  x_psopen, 'fig_dla_mtlvsz.ps', /maxs
  colors = GetColor(/Load)

  csz = 2.
  plot,[1],[1],/nodata, xrange=[0, 5.7], yrange=[-3.2,0.5], $
       xtitle='Redshift (z)', ytitle='Metallicity ([M/H])', thick=5, $
       /xstyle, /ystyle, charsiz=csz, xmarg=[9,2], ymarg=[5,2]
;     /xstyle, /ystyle, charsize=1.3

  xvals = [0,1,2,3,4,5]
  yvals = [0.6,0.6,0.6,0.6,0.6,0.6]
  xyouts, xvals-0.2, yvals - 0.05, $
          [' 0Gyr', '7.7Gyr', '10.3Gyr', '11.4Gyr', '12Gyr', '12.3Gyr'], charthick=2, charsiz=csz



  plotsym, 8, 1.2, /fill
  oplot, sdla.zabs, sdla.mtl, psym=8, color=colors.red, thick=4



  oplot, hdla.zabs, hdla.mtl, psym=8, color=colors.red, thick=4

  oplot, rdla.zabs, rdla.mtl, psym=8, color=colors.red, thick=4
  
  plotsym, 0, 1.2, /fill
  oploterror, bin, logmtl, lobin, logmtl_lo, /lobar, psym=8, color=colors.blue, errcolor=colors.blue, thick=5, errthick=5
  oploterror, bin, logmtl, hibin, logmtl_hi, /hibar, psym=8, color=colors.blue, errcolor=colors.blue, thick=5, errthick=5
  
  
  x = findgen(8001)/1000d
  
  oplot, x, fit[1]*x+fit[0], thick=5, linestyle=2, color=colors.blue
  x_psclose
  !p.multi = [0,1,1]

end
