! *********************************************************************
module module_ufile

contains

  subroutine rewnml(                                                            &
  &                       ifile,  jfile)
  use module_zcfile
  implicit none
  integer, intent(out) :: ifile,  jfile

  ifile = nfparm
  jfile = nfstdo
  rewind(ifile, err=999)
  return

  999  write(jfile, *) '*** ERROR IN REWINDING THE NAMELIST FILE ***'
  return
  end subroutine rewnml

end module module_ufile
! *********************************************************************
