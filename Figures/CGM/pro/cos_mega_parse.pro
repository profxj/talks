;; Parses the mega structure and passes back a structure of EW values
;; and Column Densities
;; flags, etc.
;; v1.1
function cos_mega_parse, lambda, MEGA=mega, NOTWOSIG=notwosig

  if not keyword_set(MEGA) then begin
     ;; Read in structure
     ldir = getenv('DROPBOX_DIR')+'/COS-Halos/lowions/'
     tdir = getenv('DROPBOX_DIR')+'/COS-Halos/Targets/'
     
     restore, ldir+'/cosmetals_megastructure.sav' 
     mega = megastruct
  endif
  ngal = n_elements(mega)

  ;; Output structure
  tmp = { $
        field: '', $
        galid: '', $
        oflg: 0, $ ; Original flag
        EW: 0., $
        sigEW: 0., $
        trueEW: 0., $ ; No monkeying with it for limits
        truesigEW: 0., $ ; No monkeying with it for limits
        flgEW: 0, $ ; Flag (0=nothing, 1=detection, 2=lower limit, 3=upper limit, -1=No good)
        colm: 0., $
        sigcolm: 0., $
        flgcolm: 0 $ ; Flag (0=nothing, 1=detection, 2=lower limit, 3=upper limit, -1=No good)
        }

  strct = replicate(tmp, ngal)
  strct.field = mega.galaxy.field
  strct.galid = mega.galaxy.galid
        
  ;; Grab EW values from the structure and handle limits
  FOR ii = 0, ngal - 1 do begin
     ;; Identify the line
     idx = where(abs(mega[ii].ion.trans.lambda - lambda) lt 1e-3, ntrans)
     case ntrans of 
        0:  ;; No such measurement
        1: begin ;; Fill in structure
           strct[ii].oflg = long( (mega[ii].ion.trans.flg)[idx])
           if (strct[ii].oflg MOD 2) EQ 0 then begin
              print, 'mega_parse: Even flag.  Ignoring', strct[ii].oflg
           endif else begin
              strct[ii].trueEW = (mega[ii].ion.trans.wrest)[idx]
              strct[ii].truesigEW = (mega[ii].ion.trans.sigwrest)[idx]
              case strct[ii].oflg of
                 1: begin ;; Standard detection
                    strct[ii].flgEW = 1
                    strct[ii].EW = (mega[ii].ion.trans.wrest)[idx]
                    strct[ii].sigEW = (mega[ii].ion.trans.sigwrest)[idx]
                    ;; Column
                    strct[ii].flgcolm= 1
                    strct[ii].colm= (mega[ii].ion.trans.logn)[idx]
                    strct[ii].sigcolm= (mega[ii].ion.trans.siglogn)[idx]
                 end
                 3: begin ;; Blended, but treated as standard
                    strct[ii].flgEW = 1
                    strct[ii].EW = (mega[ii].ion.trans.wrest)[idx]
                    strct[ii].sigEW = (mega[ii].ion.trans.sigwrest)[idx]
                    strct[ii].trueEW = strct[ii].EW 
                    strct[ii].truesigEW = strct[ii].sigEW 
                    ;; Column
                    strct[ii].flgcolm= 1
                    strct[ii].colm= (mega[ii].ion.trans.logn)[idx]
                    strct[ii].sigcolm= (mega[ii].ion.trans.siglogn)[idx]
                 end
                 5: begin ;; Undetected. Take greater of 2sigma upper limit and the value
                    strct[ii].flgEW = 3
                    strct[ii].EW = 2.*(mega[ii].ion.trans.sigwrest)[idx] > $ 
                                   (mega[ii].ion.trans.wrest)[idx]
                    strct[ii].sigEW = 99.
                    ;; Column
                    strct[ii].flgcolm= 3
                    strct[ii].colm= (mega[ii].ion.trans.logn2sig)[idx]
                    strct[ii].sigcolm= 99.
                 end
                 7: begin ;; Blended.  Take greater of 2 sigma and the value
                    strct[ii].flgEW = 3
                    strct[ii].EW = 2.*(mega[ii].ion.trans.sigwrest)[idx] > $ 
                                   (mega[ii].ion.trans.wrest)[idx]
                    strct[ii].sigEW = 99.
                    ;; Column
                    strct[ii].flgcolm= 3
                    strct[ii].colm= (mega[ii].ion.trans.logn2sig)[idx]
                    strct[ii].sigcolm= 99.
                 end
                 9: begin ;; Saturated (column)
                    strct[ii].flgEW = 1
                    strct[ii].EW = (mega[ii].ion.trans.wrest)[idx]
                    strct[ii].sigEW = (mega[ii].ion.trans.sigwrest)[idx]
                    ;; Column
                    strct[ii].flgcolm= 2
                    strct[ii].colm= (mega[ii].ion.trans.logn)[idx]
                    strct[ii].sigcolm= 99.
                 end
                 11: begin ;; Saturated and blended
                    strct[ii].flgEW = 2
                    strct[ii].EW = (mega[ii].ion.trans.wrest)[idx]
                    strct[ii].sigEW = (mega[ii].ion.trans.sigwrest)[idx]
                    ;; Column
                    strct[ii].flgcolm= 2
                    strct[ii].colm= (mega[ii].ion.trans.logn)[idx]
                    strct[ii].sigcolm= 99.
                 end
                 13: begin
                    strct[ii].flgEW = -1 ;; Saturated and blended and upper limit!!
                    strct[ii].flgcolm = -1 ;; Saturated and blended and upper limit!!
                 end
                 else: stop             ;; Flag not included yet
              endcase
           endelse
        end
        else: stop  ;; Multiple.  Likely a problem
     endcase
  endfor

  return, strct

end
