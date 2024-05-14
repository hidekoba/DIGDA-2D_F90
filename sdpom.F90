!######################################################################

! --- information -----------------------------------------------------
!
!  Sediment model : POM
!
!  HISTORY
!     '15.01.09  H.Kobayashi :
!
! ---------------------------------------------------------------------

module globals_sdpom
      use zcfile_module
      use sddim_module
      implicit none

!     REAL*8 SPOM(NXY, SDZDIM)
!     REAL*8 CDO2(NXY, SDZDIM)
!     REAL*8 RPOM(NXY, SDZDIM)
!     REAL*8 RPOMMAX

      REAL*8 KOX, KAX
      SAVE   KOX, KAX
      NAMELIST /NMSDPM/ KOX, KAX
      DATA   KOX, KAX &
     &    / 2.00D-9, 7.90D-11 /
          !! 2.0*10-9, 7.9*10-11 [/s]

!     REAL*8 CP, OP

!     REAL*8 PHI(SDZDIM)

end module globals_sdpom

! =====================================================================

module sdpom_module
      use globals_sdpom
      implicit none
contains

      SUBROUTINE SDPOMSET

      WRITE(NFSTDO, *) '*** Carbon : SDPOMSET ***'
      REWIND(NFPARM)
      READ (NFPARM, NMSDPM, END=99)
 99   WRITE(NFSTDO, NMSDPM)

      END SUBROUTINE SDPOMSET

      SUBROUTINE DISPOM( &
     &      SPOM, CDO2, CP, OP, PHI, &
     &      RPOM )

      double precision, intent(in) :: SPOM(NXY, SDZDIM)
      double precision, intent(in) :: CDO2(NXY, SDZDIM)
      double precision, intent(in) :: CP, OP
      double precision, intent(in) :: PHI(SDZDIM)
      double precision, intent(out) :: RPOM(NXY, SDZDIM)

      double precision :: RPOMMAX

      integer :: IJ, K

      DO K = SKST, SKEN
      DO IJ = IJTSTR, IJTEND
         IF ( SPOM(IJ, K) /= 0.0D0 ) THEN
            IF ( CDO2(IJ, K) > 0.0D0 ) THEN
               RPOM(IJ, K) = KOX
            ELSEIF ( CDO2(IJ, K) == 0.0D0 ) THEN
               RPOM(IJ, K) = KAX
            END IF
            RPOMMAX = ( PHI(K) / ( 1.0D0 - PHI(K) ) ) &
     &              * ( CP / OP ) * CDO2(IJ, K) / SPOM(IJ, K) / TS
            RPOM(IJ, K) = MAX( MIN( RPOM(IJ, K), RPOMMAX ), KAX )
         ELSE
            RPOM(IJ, K) = 0.0D0
         END IF
      END DO
      END DO

      END SUBROUTINE DISPOM

end module sdpom_module
