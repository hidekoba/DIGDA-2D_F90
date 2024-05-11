!######################################################################
module sdopal_module
      implicit none
contains

      SUBROUTINE SDOPALSET

      IMPLICIT NONE

! --- information -----------------------------------------------------
!
!  Sediment model : Opal
!
!  HISTORY
!     '15.01.09  H.Kobayashi :
!
! ---------------------------------------------------------------------

#include "zcfile.F"

#include "sddim.F"

      REAL*8 SOPA(NXY, SDZDIM)
      REAL*8 PCLA(NXY, SDZDIM)
      REAL*8 POPA(NXY, SDZDIM)
      REAL*8 PCO(NXY, SDZDIM)
      REAL*8 GAMMAL(NXY, SDZDIM)
      REAL*8 CSIL(NXY, SDZDIM)
      REAL*8 SATSIL(NXY, SDZDIM)
      REAL*8 ROPA(NXY, SDZDIM)

      REAL*8 UO(NXY, SDZDIM)
      REAL*8 TEM(NXY)
      REAL*8 ETA1(NXY, SDZDIM), ETA2(NXY, SDZDIM)
      REAL*8 KO(NXY, SDZDIM)
      REAL*8 RCOPA
      SAVE   RCOPA
      NAMELIST /NMSDOP/ RCOPA, ROPRID, ROPARC
      DATA   RCOPA
     &    / 5.00D-3 /
          !! 5.0*10-3 [cm3/mol/s]

      LOGICAL ROPRID, ROPARC
      SAVE ROPRID, ROPARC
      DATA ROPRID / .FALSE. /  !! Ridgwell's DEI formulation ?
      DATA ROPARC / .TRUE. /  !! Archer's DEI formulation ?

      INTEGER IJ, K

      WRITE(NFSTDO, *) '*** Carbon : SDOPALSET ***'
      REWIND(NFPARM)
      READ (NFPARM, NMSDOP, END=99)
 99   WRITE(NFSTDO, NMSDOP)

      RCOPA  = RCOPA * 1.0D-9  !! [cm3/mol/s] --> [m3/mMol/s]

      END SUBROUTINE SDOPALSET

! =====================================================================
      SUBROUTINE DISOPA(
     I      SOPA, CSIL, PCLA, POPA, TEM,
     O      ROPA, SATSIL )

      DO K = SKST, SKEN
      DO IJ = IJTSTR, IJTEND
         IF ( SOPA(IJ, K) == 0.0D0 ) THEN
            GAMMAL(IJ, K) = 0.2D0
         ELSE
            PCO(IJ, K) = PCLA(IJ, K) / POPA(IJ, K)
            IF ( PCO(IJ, K) > 1.5D+1 ) THEN
               GAMMAL(IJ, K) = 0.2D0
            ELSE
               GAMMAL(IJ, K) = 1.0D0 - ( ( 4.5D-2 * PCO(IJ, K) )
     &                       ** 5.8D-1 )
            ENDIF
         ENDIF
         SATSIL(IJ, K) = ( GAMMAL(IJ, K) * 1.0D+1
     &      ** ( 6.44D0 - ( 9.68D+2 / ( 2.7315D+2 + TEM(IJ) ) ) ) )
      ENDDO
      ENDDO

      DO K = SKST, SKEN
      DO IJ = IJTSTR, IJTEND
         IF ( SOPA(IJ, K) /= 0.0D0 .and.
     &        CSIL(IJ, K) < SATSIL(IJ, K) ) THEN
            IF (ROPRID) THEN
               UO(IJ, K)   = 1.0D0 - ( CSIL(IJ, K) / SATSIL(IJ, K) )
               ETA1(IJ, K) = 2.25D-1 * ( 1.0D0 + ( TEM(IJ) / 1.5D+1 ) )
     &                     * UO(IJ, K) +  7.75D-1
     &                     * ( ( ( ( 1.0D0 + ( TEM(IJ) / 4.0D+2 ) )
     &                     ** 4.0D0 ) * UO(IJ, K) ) ** 9.25D0 )
               ETA2(IJ, K) = 2.6D-1 + 7.4D-1
     &                     * EXP( - LZ(K) / 7.0D-2 )
               KO(IJ, K)   = 5.0D-2 + 5.5D-2
     &                     * ( ( 1.64D-2 + PCO(IJ, K) )
     &                     ** ( - 7.5D-1 ) )
               ROPA(IJ, K) = ETA1(IJ, K) * ETA2(IJ, K) * KO(IJ, K)
     &                     / 3.1536D+7
            ELSEIF (ROPARC) THEN
               ROPA(IJ, K) = RCOPA
     &             * ( SATSIL(IJ, K) - CSIL(IJ, K) )
            ENDIF
         ELSE
            ROPA(IJ, K) = 0.0D0
         ENDIF
      ENDDO
      ENDDO

      END SUBROUTINE DISOPA

end module sdopal_module
