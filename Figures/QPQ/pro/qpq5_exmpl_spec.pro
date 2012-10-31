;+ 
; NAME:
; fuse_velplt
;  V1.1
;
;   lbg_trans_ew, 1548.195, WINDEW=1., VMNX=[-300, 400]
;------------------------------------------------------------------------------
pro qpq5_exmpl_spec, WINDEW=windEW, VMNX=vmnx

  compile_opt strictarr

  cd, getenv('QSO_DIR')+'tex/QPQ7/Analysis/AbsLines/pro/', curr=curr
  resolve_routine, 'qpq7_mkdata', /compile_full_file, /is_function
  resolve_routine, 'qpq7_grab_idlines', /compile_full_file, /is_function
  cd, curr
  cd, getenv('QSO_DIR')+'tex/QPQ6/Analysis/pro/', curr=curr
  resolve_routine, 'qpq6_get_instrument', /compile_full_file, /is_function
  cd, curr
;  
  if not keyword_set( LSZ ) then lsz = 1.4
  if not keyword_set( LTHICK ) then lthick = 4.
  if not keyword_set( PTHICK ) then pthick = 2.
  if not keyword_set( BLSZ ) then blsz = 1.7
  if not keyword_set( CSIZE ) then csize = 2.5

  if not keyword_set( CONTI ) then conti = 0.8

  if not keyword_set(VMNX) then vmnx = [-2000., 2000]
  ymnx = [0., 1.2]

  ;; Transitions
  wrest = [1215.6701, 1334.5323, 1548.195d]
  clbl = ['HI Ly!9a!X', 'CII 1334', 'CIV']
  nlin = 3L
  c = x_constants()
  spl = c.c / 1d5

  ;; Grab transition info

  ;; Read in pair structure
  qpq_fil = '~/Dropbox/QSOPairs/qpq5_pairs.fits'
  qpq5_strct = xmrdfits(qpq_fil, 1)
  npairs = n_elements(qpq5_strct)

  if not keyword_set(PSFILE) then psfile = 'qpq5_exmpl_spec.ps'

  ;; Read in QPQ
  fg_qso = ['APOJ1553+1921', 'APOJ1413+2715', 'APOJ1112+6611']
  lfg_qso = ['J1553+1921', 'J1413+2715', 'J1112+6611']
  nobj = 3

  ;;  Stack files
  npx = nobj
  npy = nlin
  if keyword_set( PSFILE ) then x_psopen, psfile, /maxs

  ;;starting at p.multi[0], plot p.multi[1] rows and p.multi[2] columns
  ;;with p.multi[3] z dimensions and going from top to bottom (p.multi[4]
  !p.multi=[npx*(npy+1),npx,npy+1,0,1]
  clr = getcolor(/load)
  fclr = clr.white

  ;; LOOP
  for qq=0L,nobj-1 do begin
     ;!p.multi=[(npx-qq)*(npy+1),npx,npy+1,0,1]

     ;; Find the QSO
     mtq = where(strmatch(strtrim(qpq5_strct.qso,2), fg_qso[qq]), nmtq)
     if nmtq NE 1 then stop
     id_strct = qpq7_grab_idlines(qpq5_strct[mtq]) 

     ;; Loop on Transitions
     for ii = 0L,nlin-1 do begin
        !p.multi=[(npx-qq)*(npy+1) - ii,npx,npy+1,0,1]
        
        ;; Y-axis
        if ymnx[1]-ymnx[0] GT 0.8 then ytint = 0.5 else begin
           if ymnx[1]-ymnx[0] GT 0.4 then ytint = 0.2 else ytint = 0.1
        endelse
        
        ;; Plot
        ysty=1
        if ii LT (nlin-1) then begin
           xspaces = replicate(' ',30) 
           xtitle = ''
        endif else begin 
           if keyword_set(XSPACES) then delvarx, xspaces
           xtitle='Relative Velocity (km s!u-1!N)'
        endelse
        if ii EQ 1 then ytit='Normalized Flux' else ytit = ''
        if (ii GT npy-1) then yspaces = replicate(' ',30) 
        
        plot, [0],  [0],  xrange=vmnx, $
              yrange=ymnx, xtickn=xspaces, xmargin=[7,0], $
              ymargin=[0,0], NODATA=nblnd, $ ;ytickn=yspaces, $
              charsize=csize, psym=10, background=clr.white, color=fclr, $
              xstyle=1, ystyle=ysty, thick=pthick, ytickinterval=ytint,$
              xtitle=xtitle, ytitle=ytit

        ;; Read in data
        if ii GT 0 then begin
           idLine = where(abs(id_strct.wrest-wrest[ii]) LT 1e-3, nid)
           if nid NE 1 then stop
           datfil = getenv('QSO_DIR')+id_strct[idLine].datfil
        endif else datfil = getenv('QSO_DIR')+qpq5_strct[mtq].lya_fil
           data = qpq7_mkdata(qpq5_strct[mtq], FLG=flg_data, datfil=datfil)

        ;; Gray Redshift error
        x_curvefill, [-1,1]*qpq5_strct[mtq].z_fsig, replicate(ymnx[0],2), replicate(ymnx[1],2), $
                     color=clr.darkgray

        ;; Plot the Data
        velo = (data.wave - wrest[ii] * (1+qpq5_strct[mtq].z_fg))*spl / $
               (wrest[ii] * (1+qpq5_strct[mtq].z_fg))
        oplot, velo, data.fx, color=fclr, psym=10, thick=4

        ;; Lines
        ;oplot, [0., 0.], ymnx, color=clr.cyan, linestyle=2, thick=3
        oplot, [-10000., 10000.], [1.,1.], color=clr.green, linestyle=3, thick=4

        ;; Label
        case ii of
           0: xyouts, vmnx[0]+0.5*(vmnx[1]-vmnx[0]), 1.07, lfg_qso[qq], color=fclr, charsi=lsz, align=0.5
           1: xyouts, vmnx[0]+0.03*(vmnx[1]-vmnx[0]), 0.30, $
                      'R='+strtrim(round(qpq5_strct[mtq].R_phys),2)+' kpc', $
                      color=fclr, charsi=lsz, align=0.
           else:
        endcase

        xyouts, vmnx[0]+0.03*(vmnx[1]-vmnx[0]), 0.06, clbl[ii], color=clr.yellow, charsi=lsz, align=0.
        instr = qpq6_get_instrument( datfil, FULLSPEC=fullspec )
        xyouts, vmnx[1]-0.25*(vmnx[1]-vmnx[0]), 0.06, instr, color=clr.cyan, charsi=lsz, align=0.

        !p.multi=[(npx-qq)*(npy+1) - ii,npx,npy+1,0,1]
        plot, [0],  [0],  xrange=vmnx, $
              yrange=ymnx, xtickn=xspaces, xmargin=[7,0], $
              ymargin=[0,0], NODATA=nblnd, $ ;ytickn=yspaces, $
              charsize=csize, psym=10, background=clr.white, color=fclr, $
              xstyle=1, ystyle=ysty, thick=pthick, ytickinterval=ytint,$
              xtitle=xtitle, ytitle=ytit, /noeras
     endfor
  endfor 
  
  ;xyouts, 0.03, 0.55, 'Normalized Flux', $
  ;        alignment=0.5, ORIENTATION=90., /normal, charsize=blsz

  if keyword_set( PSFILE ) then x_psclose

  !p.multi=[0,1,1]
  
  return
end
