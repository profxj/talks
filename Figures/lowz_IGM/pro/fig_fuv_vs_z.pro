;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro fig_fuv_vs_z, summ_fil, ALL_GAL=all_gal

  if not keyword_set(psfile) then psfile = 'fig_fuv_vs_z.ps'
  if not keyword_set(binsz) then binsz = 0.01
  if not keyword_set(lsz) then lsz = 2.2 
  if not keyword_set(lsz2) then lsz2 = 2.5 
  if not keyword_set(csz) then csz = 2.5 
  if not keyword_set(dlim) then dlim = 5. ;arcmin
  ;; Blanton lum function
  if not keyword_set(phi_str) then phi_str = 1.49  ; 10^-2 h^3 Mpc^-3
  if not keyword_set(alpha) then alpha = -1.05
  if not keyword_set(Mstar) then Mstar = -20.44  ; M* - 5 log h

  ;; Initialize
  restore, '~/Dropbox/COS-Halos/idl/Catalogs/dr7_qso_galex.dat', /verb


  c = x_constants()


  ;; Plot
  if keyword_set(psfile) then x_psopen,psfile,/maxs
  !p.multi=[0,1,1]
  clr = getcolor(/load)
  lclr = clr.lightgray
  ;lclr = clr.black
  pclr = clr.cyan

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; dN/dX
  yrng = [20., 14.5]
  xrng = [0.05, 2] 
  xmrg = [8,6.0]
  ylbl = -1


  for ss=0,4 do begin
  
     case ss of 
        0: FUV_lim = 99.
        1: begin ;; STIS 
           FUV_exp = [ [15., 3958.], $
                       [16., 7387.], $
                       [17., 25078.], $
                       [18., 104140.] ]
           FUV_orbit = FUV_exp
           FUV_orbit[1,*] = FUV_orbit[1,*] / 2400. 
           ysty = 9
           lbl = 'STIS'
        end
        2: FUV_lim = 17.2 ; STIS again but blocked out
        3: begin ;; COS
           FUV_exp = [ [15., 136.], $
                       [16., 344.], $
                       [17., 876.], $
                       [18., 2266.], $
                       [19., 6100.] ]
           FUV_orbit = FUV_exp
           FUV_orbit[1,*] = FUV_orbit[1,*] / 2400. 
           ysty = 9
           lbl = 'COS'
        end
        4: FUV_lim = 19.2
        else: 
     endcase
;  xspaces = replicate(' ',30) 
     plot,[0.],[0.],xrange=xrng,yrange=yrng,color=lclr,$
          background=clr.white, charsize=csz, $
          xmargin=xmrg, ymargin=[4.5,3],xtitle='z!dQSO!N', $
          ytitle='FUV (mag)',/nodata,xstyle=1,ystyle=1, /xlog
  
     ;; Quasars
     plotsym, 0, 0.4, /fill
     gd = where(qso.FUV_MAG LT FUV_lim)
     oplot, qso.z[gd], qso.FUV_MAG[gd], psym=8, color=pclr

     ;; Label orbits
     if ss GT 0 then begin
        sz = size(FUV_orbit, /dime)
        xpos = 1.19
        for kk=0L,sz[1]-1 do begin
           norb = round(FUV_orbit[1,kk]) 
           if norb EQ 0 then ast = '<1'  $
           else ast = string(round(FUV_orbit[1,kk]),format='(i3)')
           xyouts, xrng[1]*xpos, FUV_orbit[0,kk]+0.1, ast, charsiz=lsz, color=lclr, align=0.5
        endfor
        xyouts, 0.07, yrng[1]+0.5, lbl, charsiz=lsz2, color=lclr, align=0.
        xyouts, xrng[1]*xpos, yrng[1]-0.4, 'S/N=7', charsiz=lsz, color=lclr, align=0.5
        xyouts, xrng[1]*xpos, yrng[1]-0.1, '(15 km/s)', charsiz=lsz, color=lclr, align=0.5
     endif

     if ss GE 2 then begin
        x_curvefill, xrng, replicate(FUV_lim,2),  replicate(yrng[0],2), color=clr.gray
     endif
     if ss EQ 4 then begin
        xyouts, 0.35, 19.7, 'Game Changer', color=lclr, charsi=lsz2, align=0.5
     endif
     
              
  ;;  Label orbits on the RHS
  ;lblt = -5 + lindgen(5)
  ;nlbl = n_elements(lblt)
  ;for jj=0L,nlbl-1 do begin
  ;   xyouts, lblt[jj]-0.05, ymnx[1]+0.2, string(NH_lbl[jj],format='(f5.1)'), $
  ;           color=clr.black, charsize=csz, alignment=0.5
  ;endfor

  endfor

  x_psclose
  !p.multi=[0,1,1]

  return
end
