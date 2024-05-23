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

      subroutine SDCLAYSET

      WRITE(NFSTDO, *) '*** Carbon : SDCLAYSET ***'

      end subroutine SDCLAYSET

      subroutine DISCLA( &
     &      SCLA, &
     &      RCLA )

      double precision, intent(in) :: SCLA(NXY, SDZDIM)
      double precision, intent(out) :: RCLA(NXY, SDZDIM)

      integer :: IJ, K

      DO K = SKST, SKEN
         DO IJ = IJTSTR, IJTEND
            RCLA(IJ, K)  = 0.0D0
         END DO
      END DO

      end subroutine DISCLA

end module module_digda_clay
