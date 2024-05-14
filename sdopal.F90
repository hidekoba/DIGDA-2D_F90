!######################################################################

! --- information -----------------------------------------------------
!
!  Sediment model : Opal
!
!  HISTORY
!     '15.01.09  H.Kobayashi :
!
! ---------------------------------------------------------------------

module globals_sdopal
      use zcfile_module
      use sddim_module
      implicit none

!     REAL*8 SOPA(NXY, SDZDIM)
!     REAL*8 PCLA(NXY, SDZDIM)
!     REAL*8 POPA(NXY, SDZDIM)
!     REAL*8 PCO(NXY, SDZDIM)
!     REAL*8 CSIL(NXY, SDZDIM)
!     REAL*8 SATSIL(NXY, SDZDIM)
!     REAL*8 ROPA(NXY, SDZDIM)

!     REAL*8 UO(NXY, SDZDIM)
!     REAL*8 TEM(NXY)
!     REAL*8 ETA1(NXY, SDZDIM), ETA2(NXY, SDZDIM)
!     REAL*8 KO(NXY, SDZDIM)

      REAL*8 RCOPA
      SAVE   RCOPA
      LOGICAL ROPRID, ROPARC
      SAVE ROPRID, ROPARC
      NAMELIST /NMSDOP/ RCOPA, ROPRID, ROPARC
      DATA   RCOPA &
     &    / 5.00D-3 /
          !! 5.0*10-3 [cm3/mol/s]
      DATA ROPRID / .FALSE. /  !! Ridgwell's DEI formulation ?
      DATA ROPARC / .TRUE. /  !! Archer's DEI formulation ?

end module globals_sdopal

! =====================================================================

module sdopal_module
      use globals_sdopal
      implicit none
contains

      SUBROUTINE SDOPALSET

      WRITE(NFSTDO, *) '*** Carbon : SDOPALSET ***'
      REWIND(NFPARM)
      READ (NFPARM, NMSDOP, END=99)
 99   WRITE(NFSTDO, NMSDOP)

      RCOPA  = RCOPA * 1.0D-9  !! [cm3/mol/s] --> [m3/mMol/s]

      END SUBROUTINE SDOPALSET

      SUBROUTINE DISOPA( &
     &      SOPA, CSIL, PCLA, POPA, TEM, &
     &      ROPA, SATSIL )

      double precision, intent(in) :: SOPA(NXY, SDZDIM)
      double precision, intent(in) :: CSIL(NXY, SDZDIM)
      double precision, intent(in) :: PCLA(NXY, SDZDIM)
      double precision, intent(in) :: POPA(NXY, SDZDIM)
      double precision, intent(in) :: TEM(NXY)
      double precision, intent(out) :: ROPA(NXY, SDZDIM)
      double precision, intent(out) :: SATSIL(NXY, SDZDIM)

      double precision :: GAMMAL(NXY, SDZDIM)
      double precision :: PCO(NXY, SDZDIM)
      double precision :: UO(NXY, SDZDIM)
      double precision :: ETA1(NXY, SDZDIM), ETA2(NXY, SDZDIM)
      double precision :: KO(NXY, SDZDIM)

      integer :: IJ, K

      DO K = SKST, SKEN
      DO IJ = IJTSTR, IJTEND
         IF ( SOPA(IJ, K) == 0.0D0 ) THEN
            GAMMAL(IJ, K) = 0.2D0
         ELSE
            PCO(IJ, K) = PCLA(IJ, K) / POPA(IJ, K)
            IF ( PCO(IJ, K) > 1.5D+1 ) THEN
               GAMMAL(IJ, K) = 0.2D0
            ELSE
               GAMMAL(IJ, K) = 1.0D0 - ( ( 4.5D-2 * PCO(IJ, K) ) &
     &                       ** 5.8D-1 )
            END IF
         END IF
         SATSIL(IJ, K) = ( GAMMAL(IJ, K) * 1.0D+1 &
     &      ** ( 6.44D0 - ( 9.68D+2 / ( 2.7315D+2 + TEM(IJ) ) ) ) )
      END DO
      END DO

      DO K = SKST, SKEN
      DO IJ = IJTSTR, IJTEND
         IF ( SOPA(IJ, K) /= 0.0D0 .and. &
     &        CSIL(IJ, K) < SATSIL(IJ, K) ) THEN
            IF (ROPRID) THEN
               UO(IJ, K)   = 1.0D0 - ( CSIL(IJ, K) / SATSIL(IJ, K) )
               ETA1(IJ, K) = 2.25D-1 * ( 1.0D0 + ( TEM(IJ) / 1.5D+1 ) ) &
     &                     * UO(IJ, K) +  7.75D-1 &
     &                     * ( ( ( ( 1.0D0 + ( TEM(IJ) / 4.0D+2 ) ) &
     &                     ** 4.0D0 ) * UO(IJ, K) ) ** 9.25D0 )
               ETA2(IJ, K) = 2.6D-1 + 7.4D-1 &
     &                     * EXP( - LZ(K) / 7.0D-2 )
               KO(IJ, K)   = 5.0D-2 + 5.5D-2 &
     &                     * ( ( 1.64D-2 + PCO(IJ, K) ) &
     &                     ** ( - 7.5D-1 ) )
               ROPA(IJ, K) = ETA1(IJ, K) * ETA2(IJ, K) * KO(IJ, K) &
     &                     / 3.1536D+7
            ELSEIF (ROPARC) THEN
               ROPA(IJ, K) = RCOPA &
     &             * ( SATSIL(IJ, K) - CSIL(IJ, K) )
            END IF
         ELSE
            ROPA(IJ, K) = 0.0D0
         END IF
      END DO
      END DO

      END SUBROUTINE DISOPA

end module sdopal_module
