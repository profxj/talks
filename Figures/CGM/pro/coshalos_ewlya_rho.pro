;; Compare COS-Halos MgII EW vs. Chen
pro coshalos_ewlya_rho, lambda, PSFILE=psfile, MEGASTRUCT=megastruct, $
                     EW_range=ew_range, N_range=N_range, $
                     YMRG=ymrg, CSZ=csz, NOPS=nops, XLOG=xlog

  compile_opt strictarr

  if not keyword_set(CSZ) then csz = 2.2
  if not keyword_set(LSZ) then lsz = 2.0
  if not keyword_set(ASZ) then asz = 1.8  ;; arrow size

  if not keyword_set(EW_range) then EW_range = [0.02, 6.]
  if not keyword_set(N_range) then N_range = [12.0, 15.0]

  ;; Ion specific stuff
  lambda = 1215.6701d
  getion, lambda, ioni, elm, Z=ionz, NM=strioni

  case strioni of 
     'I': striplus = '0'
     'II': striplus = '+'
     'III': striplus = '++'
     'IV': striplus = '+3'
     'VI': striplus = '+5'
     else: stop
  endcase

  ;; PSFILE
  if not keyword_set(PSFILE) then  psfile = 'coshalos_ewlya_rho.ps'

  ;; Initialize
  cd, '~/paper/COS_HALO/LowIons/paper/Analysis/pro/', curr=curr
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
     print, 'Analyzing ', n_elements(megastruct), ' systems'
  endif
  lgsSFR = alog10(megastruct.galaxy.SFR) - megastruct.logmfinal
  printcol, lindgen(43), megastruct.galaxy.field, lgsSFR


  ;; N_HI Scaling
  NHI_min = 13.
  NHI_max = 17.
  ;sym_min = 0.4
  ;sym_max = 2.0
  sym_min = 1.5
  sym_max = 1.5
  m_sym = (sym_max-sym_min)/(NHI_max-NHI_min)
  b_sym = sym_max - m_sym * NHI_max
   
  ;; impact parameter
  if not keyword_set(XLOG) then xrange = [0.0, 170.0] else xrange = [10., 200]
  xlabel = '!8R!X [kpc]'

  ;; Plot
  x_psopen, psfile, /maxs
  !p.multi=[0,1,1]
  clr = getcolor(/load)
  lclr = clr.lightgray

  xmrg = [8,15.7]
  if not keyword_set(YMRG) then ymrg = [7,1]


     ylabel = '!8EW!X!d'+elm+strioni+' '+ $
              strtrim(string(fix(lambda), format = '(i4)'), 2)+'!n [A]'
     yrange = EW_range
     ylg = 1

  for ss=0,5 do begin

     plot, [0], [0], /nodata, yr = yrange, ys = 1, xr = xrange, xs = 1, $
           ytit = ylabel, xtit = xlabel, col = lclr, $
           chars = csz, xmar=xmrg, ymar=ymrg, ylog=ylg, XLOG=xlog
        
     if ss EQ 0 then continue

     case ss of 
        1: strct = megastruct[13]
        2: strct = megastruct[[13,30]]
        3: strct = megastruct[[13,30,38]]
        4: strct = megastruct[[13,30,38,0]]
        5: strct = megastruct
        else: stop
     endcase
     ewstrct = mega_parse_trans(lambda, mega=strct)

     cos_halos_redblue, strct, red, blue, NRED=nred, NBLUE=nblue, FLAG=0
  
     for qq=0,1 do begin ;; Color
        case qq of
           0: begin
              idx = blue
              pclr = clr.cyan
              symc = 5 ;; Filled square
              symo = 6
              nobj = nblue
           end
           1: begin
              pclr = clr.tomato
              idx = red
              symc = 4
              symo = 4
              nobj = nred
           end
           else: stop
        endcase

        if nobj EQ 0 then continue
     
        for kk=1,3 do begin ;; Flag
           flgs = ewstrct[idx].flgEW
           vals = ewstrct[idx].EW/1e3
           sigs = ewstrct[idx].sigEW/1e3
           pidx = where(flgs EQ kk, npid)
           for jj=0L,npid-1 do begin
              syms = strct[idx[pidx[jj]]].ion[1,1].lognion*m_sym + b_sym
              syms = sym_min > syms < sym_max
              
              ;; Filled
              if kk NE 3 then begin
                 ;; Value
                 oplot, [strct[idx[pidx[jj]]].rhofinal], $
                        [vals[pidx[jj]]], $
                        psym=sym(symc), color=pclr, symsi=syms ;, thick=5
                 ;; Outline
                 oplot, [strct[idx[pidx[jj]]].rhofinal], $
                        [vals[pidx[jj]]], $
                        psym=symo, color=lclr, symsi=syms*1.0, thick=6
              endif else begin
                 ;; Outline
                 oplot, [strct[idx[pidx[jj]]].rhofinal], $
                        [vals[pidx[jj]]], $
                        psym=symo, color=pclr, symsi=syms*1.0, thick=6
              endelse
              if kk NE 1 then begin
                 ;; Arrow
                 if kk EQ 2 then plotsym, 2, asz, thick=5 else $
                    plotsym, 1, asz, thick=5
                 oplot, [strct[idx[pidx[jj]]].rhofinal], $
                        [vals[pidx[jj]]], psym=8, color=lclr
              endif
           endfor
        endfor
     endfor
  endfor

  if not keyword_set(NOPS) then begin
     x_psclose
     !p.multi=[0,1,1]
  endif

  return
end
