!######################################################################
module sdclay_module
      implicit none
contains

      SUBROUTINE SDCLAYSET

      IMPLICIT NONE

! --- information -----------------------------------------------------
!
!  Sediment model : Clay
!
!  HISTORY
!     '15.01.09  H.Kobayashi :
!
! ---------------------------------------------------------------------

#include "zcfile.F"

#include "sddim.F"


      WRITE(NFSTDO, *) '*** Carbon : SDCLAYSET ***'

      END SUBROUTINE SDCLAYSET

! =====================================================================
      SUBROUTINE DISCLA(
     I      SCLA,
     O      RCLA )

      DO K = SKST, SKEN
         DO IJ = IJTSTR, IJTEND
            RCLA(IJ, K)  = 0.0D0
         END DO
      END DO

      END SUBROUTINE DISCLA

end module sdclay_module
