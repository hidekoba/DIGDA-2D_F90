!######################################################################

! --- information -----------------------------------------------------
!
!  Sediment model : POM
!
!  HISTORY
!     '15.01.09  H.Kobayashi :
!
! ---------------------------------------------------------------------

module globals_digda_pom
  use module_zcfile
  use module_digda_grid
  use module_sdtime
  implicit none

      REAL*8 KOX, KAX
      SAVE   KOX, KAX
      NAMELIST /NMSDPM/ KOX, KAX
      DATA   KOX, KAX &
     &    / 2.00D-9, 7.90D-11 /
          !! 2.0*10-9, 7.9*10-11 [/s]

end module globals_digda_pom

! =====================================================================

module module_digda_pom
  use globals_digda_pom
  implicit none
contains

      subroutine SDPOMSET

      WRITE(NFSTDO, *) '*** Carbon : SDPOMSET ***'
      REWIND(NFPARM)
      READ (NFPARM, NMSDPM, END=99)
99 continue
      WRITE(NFSTDO, NMSDPM)

      end subroutine SDPOMSET

      subroutine DISPOM( &
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

      end subroutine DISPOM

end module module_digda_pom
