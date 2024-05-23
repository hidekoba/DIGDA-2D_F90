!######################################################################

! --- information -----------------------------------------------------
!
!  Sediment model : Clay
!
!  HISTORY
!     '15.01.09  H.Kobayashi :
!
! ---------------------------------------------------------------------

module globals_digda_clay
  use module_zcfile
  use module_digda_grid
  implicit none

end module globals_digda_clay

! =====================================================================

module module_digda_clay
  use globals_digda_clay
  implicit none
contains

  subroutine set_digda_clay

    WRITE(NFSTDO, *) '*** Sediment : set_digda_clay ***'

  end subroutine set_digda_clay

  subroutine calc_clay_dissolution( &
 &  SCLA, &
 &  RCLA )

    double precision, intent(in) :: SCLA(NXY, SDZDIM)
    double precision, intent(out) :: RCLA(NXY, SDZDIM)

    integer :: IJ, K

    DO K = SKST, SKEN
      DO IJ = IJTSTR, IJTEND
        RCLA(IJ, K)  = 0.0D0
      END DO
    END DO

  end subroutine calc_clay_dissolution

end module module_digda_clay
