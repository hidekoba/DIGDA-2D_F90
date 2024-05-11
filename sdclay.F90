!######################################################################

! --- information -----------------------------------------------------
!
!  Sediment model : Clay
!
!  HISTORY
!     '15.01.09  H.Kobayashi :
!
! ---------------------------------------------------------------------

module globals_sdclay
      use zcfile_module
      use sddim_module
      implicit none

      REAL*8 SCLA(NXY, SDZDIM)
      REAL*8 RCLA(NXY, SDZDIM)

      INTEGER IJ, K

end module globals_sdclay

! =====================================================================

module sdclay_module
      use globals_sdclay
      implicit none
contains

      SUBROUTINE SDCLAYSET

      WRITE(NFSTDO, *) '*** Carbon : SDCLAYSET ***'

      END SUBROUTINE SDCLAYSET
      SUBROUTINE DISCLA( &
     &      SCLA, &
     &      RCLA )

      double precision, intent(in) :: SCLA(NXY, SDZDIM)
      double precision, intent(out) :: RCLA(NXY, SDZDIM)

      DO K = SKST, SKEN
         DO IJ = IJTSTR, IJTEND
            RCLA(IJ, K)  = 0.0D0
         END DO
      END DO

      END SUBROUTINE DISCLA

end module sdclay_module
