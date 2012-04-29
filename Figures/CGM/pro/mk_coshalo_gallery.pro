pro mk_coshalo_gallery

  ;;
  coshalo_gallery, 'J0042-1037', '358_9', IMG=img, STRCT=strct, PSFILE=psfile
  spawn, 'ps2pdf '+psfile

  coshalo_gallery, 'J1009+0713',  '204_17', PSFILE=psfile, SMAX=20.
  spawn, 'ps2pdf '+psfile

  coshalo_gallery, 'J1555+3628', '88_11', PSFILE=psfile, SMAX=20.
  spawn, 'ps2pdf '+psfile

  coshalo_gallery, 'J1419+4207', '132_30', PSFILE=psfile
  spawn, 'ps2pdf '+psfile

  ;coshalo_gallery, 'J0935+0204', '15_28' -- MagE!
  ;coshalo_gallery, 'J1342-0053', '304_29' -- MagE!

  return
end
