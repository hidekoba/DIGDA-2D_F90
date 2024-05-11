!######################################################################
module sdcaco3_module
      implicit none
contains

      SUBROUTINE SDCACO3SET

      IMPLICIT NONE

! --- information -----------------------------------------------------
!
!  Sediment model : CaCO3
!
!  HISTORY
!     '15.01.09  H.Kobayashi :
!
! ---------------------------------------------------------------------

#include "zcfile.F"
#include "zcphys.F"

#include "sddim.F"


      WRITE(NFSTDO, *) '*** Carbon : SDCACO3SET ***'
      REWIND(NFPARM)
      READ (NFPARM, NMSDCA, END=99)
 99   WRITE(NFSTDO, NMSDCA)

      RCCAR  = RCCAR / 3.1536D+7  !! [/yr] --> [/s]

      END SUBROUTINE SDCACO3SET

! =====================================================================
      SUBROUTINE DISCAR(
     I      SCAR, CCO32, SATCO32, PHI,
     O      RCAR )

      DO K = SKST, SKEN
      DO IJ = IJTSTR, IJTEND
         IF ( SCAR(IJ, K) /= 0.0D0 .and.
     &        CCO32(IJ, K) < SATCO32(IJ) ) THEN
            RCAR(IJ, K) = RCCAR
     &          * ( ( 1.0D0 - ( CCO32(IJ, K) / SATCO32(IJ) ) )
     &          ** NCAR )
         ELSE
            RCAR(IJ, K) = 0.0D0
         ENDIF
      ENDDO
      ENDDO

      END SUBROUTINE DISCAR

! =====================================================================
      SUBROUTINE SATCAR(
     I      TEM, SAL, BDEP,
     O      CSAT )

      DO IJ = IJTSTR, IJTEND
         PRES  =  BDEP(IJ) / 1.0D+3  !! [cm] --> [bar]
         KPRES = LOG(KPRIME)
     &      - DELV / ( RR * ( TEM(IJ) + kelvin ) )        * (PRES)     !! Effect of Pressure
     &      + 0.50D0 * DK / ( RR * ( TEM(IJ) + kelvin ) ) * (PRES)**2  !! Effect of Pressure
         KPRES = EXP(KPRES)
         CSAT(IJ)  = ( KPRES / 1.028D-2 ) * 1.0D+6                     !! Millero [1982]
      ENDDO

      END SUBROUTINE SATCAR(

end module sdcaco3_module
