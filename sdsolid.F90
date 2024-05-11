!######################################################################
module sdsolid_module
      implicit none
contains

      SUBROUTINE SDSOLIDSET

      IMPLICIT NONE

! --- information -----------------------------------------------------
!
!  Sediment model : Particles
!
!  HISTORY
!     '15.01.09  H.Kobayashi :
!
! ---------------------------------------------------------------------

#include "zcfile.F"

#include "sddim.F"

      REAL*8 FF(NXY)
      REAL*8 BFZ(NXY)
      REAL*8 S(NXY,SDZDIM)
      REAL*8 SI(NXY,SDZDIM+1)
      REAL*8 SN(NXY,SDZDIM)
      REAL*8 C0(NXY,SDZDIM+1)
      REAL*8 C1(NXY,SDZDIM+1)
      REAL*8 C2(NXY,SDZDIM+1)
      REAL*8 FLX(NXY,SDZDIM+1), DFLX(NXY,SDZDIM+1)
      REAL*8 W(NXY,SDZDIM+1)
      REAL*8 RST(NXY,SDZDIM)

      REAL*8 PHI(SDZDIM)
      REAL*8 PHII(SDZDIM+1)
      REAL*8 DBI(SDZDIM+1)
      REAL*8 RHOCLA

      REAL*8 EA(SDZDIM), EB(SDZDIM), EC(SDZDIM)
      REAL*8 EAL(SDZDIM), EBE(SDZDIM), EGA(SDZDIM), EEP(NXY,SDZDIM)
      REAL*8 ED(NXY,SDZDIM),  EX(NXY,SDZDIM)
      REAL*8 FLGAM(SDZDIM+1),  FLDEL(SDZDIM+1), FLTMP

      INTEGER IUNIT
      DATA    IUNIT / 1 /
      REAL*8  AUNIT
      SAVE    AUNIT

      INTEGER IJ, K

      LOGICAL SADU, SADQ, SADQF, SDIC, SDII
      SAVE SADU, SADQ, SADQF, SDIC, SDII
      NAMELIST /NMSDSO/ SADU, SADQ, SADQF, SDIC, SDII
!     -- Advection choice --
      DATA SADU  / .FALSE. / !! Upwind ?
      DATA SADQ  / .TRUE.  / !! Quickest ?
      DATA SADQF / .TRUE.  / !!   + Flux limiter ?
!     -- Diffusion choice --
      DATA SDIC  / .FALSE. / !! Crank nicolson ?
      DATA SDII  / .TRUE.  / !! Implicit ?

      WRITE(NFSTDO, *) '*** Carbon : SDSOLIDSET ***'
      REWIND(NFPARM)
      READ (NFPARM, NMSDSO, END=99)
 99   WRITE(NFSTDO, NMSDSO)

      IF (IUNIT.EQ.1) THEN
         AUNIT = 1.0D0 / 8.64D4 / 365.0D0
      ELSEIF (IUNIT.EQ.2) THEN
         AUNIT = 1.0D0 / 8.64D4 / 30.0D0
      ELSEIF (IUNIT.EQ.3) THEN
         AUNIT = 1.0D0 / 8.64D4
      ELSEIF (IUNIT.EQ.4) THEN
         AUNIT = 1.0D0 / 3600.0D0
      ELSEIF (IUNIT.EQ.5) THEN
         AUNIT = 1.0D0 / 60.0D0
      ELSEIF (IUNIT.EQ.6) THEN
         AUNIT = 1.0D0
      ELSE
         AUNIT = 0.0D0
      ENDIF

!     -- Convert unit --
      DO K = SKST, SKEN
         LZ(K)   = LZ(K)   * 1.0D+2            !! [m] --> [cm]
         DLZ(K)  = DLZ(K)  * 1.0D+2            !! [m] --> [cm]
      ENDDO
      DO K = SKST, SKEN+1
         LZI(K)  = LZI(K)  * 1.0D+2            !! [m] --> [cm]
         DLZI(K) = DLZI(K) * 1.0D+2            !! [m] --> [cm]
      ENDDO

      RETURN

! =====================================================================
      SUBROUTINE SOLID(
     I      FF, RST, W, PHI, PHII, DBI, S,
     O      SN, BFZ )

      DO IJ = IJTSTR, IJTEND
         FLX(IJ, SKST) = FF(IJ)
      ENDDO

!     -- Advection : Upwind --

      IF (SADU) THEN

      DO K = SKST+1, SKEN
      DO IJ = IJTSTR, IJTEND
         IF ( W(IJ, K) >= 0.0D0 ) THEN
            SI(IJ, K) = S(IJ, K-1)
         ELSE IF ( W(IJ, K) < 0.0D0 ) THEN
            SI(IJ, K) = S(IJ, K)
         ENDIF
         FLX(IJ, K)   = ( 1.0D0 - PHII(K) ) * W(IJ, K) * SI(IJ, K)
      ENDDO
      ENDDO

      DO IJ = IJTSTR, IJTEND
      IF ( W(IJ, SKEN+1) >= 0.0D0 ) THEN
         SI(IJ, SKEN+1) = S(IJ, SKEN)
      ELSE IF ( W(IJ, SKEN+1) < 0.0D0 ) THEN
         SI(IJ, SKEN+1) = 0.0D0
      ENDIF
         FLX(IJ, SKEN+1) = ( 1.0D0 - PHII(SKEN+1) )
     &                   * W(IJ, SKEN+1) * SI(IJ, SKEN+1)
      ENDDO

      ENDIF

!     -- Advection : Quickest --

      IF (SADQ) THEN

      DO K = SKST+2, SKEN-1
      DO IJ = IJTSTR, IJTEND
         IF ( W(IJ, K) >= 0.0D0 ) THEN
            C2(IJ, K) = ( 1.0D0 / ( DLZI(K-1) + DLZI(K) ) )
     &         * ( ( ( S(IJ, K-2) - S(IJ, K-1) ) / DLZI(K-1) )
     &           - ( ( S(IJ, K-1) - S(IJ, K)   ) / DLZI(K)   ) )
         ELSEIF ( W(IJ, K) < 0.0D0 ) THEN
            C2(IJ, K) = ( 1.0D0 / ( DLZI(K) + DLZI(K+1) ) )
     &         * ( ( ( S(IJ, K-1) - S(IJ, K) ) / DLZI(K)   )
     &           - ( ( S(IJ, K) - S(IJ, K+1) ) / DLZI(K+1) ) )
         ENDIF
            C1(IJ, K) = ( ( S(IJ, K) - S(IJ, K-1) ) / DLZI(K) )
     &         + ( ( ( DLZ(K-1) - DLZ(K) ) / 2.0D0 ) * C2(IJ, K) )
            C0(IJ, K) = ( ( DLZ(K) * S(IJ, K-1) + DLZ(K-1) * S(IJ, K) )
     &                / ( 2.0D0 * DLZI(K) ) )
     &         - ( ( ( DLZ(K-1) * DLZ(K) ) / 4.0D0 ) * C2(IJ, K) )
      ENDDO
      ENDDO

      DO IJ = IJTSTR, IJTEND
      IF ( W(IJ, SKST+1) >= 0.0D0 ) THEN
         C2(IJ, SKST+1) = ( 1.0D0 / ( DLZI(SKST) + DLZI(SKST+1) ) )
     &      * ( ( ( 0.0D0       - S(IJ, SKST)   ) / DLZI(SKST) )
     &        - ( ( S(IJ, SKST) - S(IJ, SKST+1) ) / DLZI(SKST+1) ) )
      ELSEIF ( W(IJ, SKST+1) < 0.0D0 ) THEN
         C2(IJ, SKST+1) = ( 1.0D0 / ( DLZI(SKST+1) + DLZI(SKST+2) ) )
     &      * ( ( ( S(IJ, SKST)   - S(IJ, SKST+1) ) / DLZI(SKST+1) )
     &        - ( ( S(IJ, SKST+1) - S(IJ, SKST+2) ) / DLZI(SKST+2) ) )
      ENDIF
         C1(IJ, SKST+1) = ( ( S(IJ, SKST+1) - S(IJ, SKST) )
     &                      / DLZI(SKST+1) )
     &      + ( ( ( DLZ(SKST) - DLZ(SKST+1) ) / 2.0D0 )
     &          * C2(IJ, SKST+1) )
         C0(IJ, SKST+1) = ( ( DLZ(SKST+1) * S(IJ, SKST)
     &                      + DLZ(SKST)   * S(IJ, SKST+1) )
     &                  / ( 2.0D0 * DLZI(SKST+1) ) )
     &      - ( ( ( DLZ(SKST) * DLZ(SKST+1) ) / 4.0D0 )
     &          * C2(IJ, SKST+1) )
      ENDDO

      DO IJ = IJTSTR, IJTEND
      IF ( W(IJ, SKEN) >= 0.0D0 ) THEN
         C2(IJ, SKEN) = ( 1.0D0 / ( DLZI(SKEN-1) + DLZI(SKEN) ) )
     &      * ( ( ( S(IJ, SKEN-2) - S(IJ, SKEN-1) ) / DLZI(SKEN-1) )
     &        - ( ( S(IJ, SKEN-1) - S(IJ, SKEN)   ) / DLZI(SKEN)   ) )
      ELSEIF ( W(IJ, SKEN) < 0.0D0 ) THEN
         C2(IJ, SKEN) = ( 1.0D0 / ( DLZI(SKEN) + DLZI(SKEN+1) ) )
     &      * ( ( ( S(IJ, SKEN-1) - S(IJ, SKEN) ) / DLZI(SKEN)   )
     &        - ( ( S(IJ, SKEN)   - 0.0D0       ) / DLZI(SKEN+1) ) )
      ENDIF
         C1(IJ, SKEN) = ( ( S(IJ, SKEN) - S(IJ, SKEN-1) ) / DLZI(SKEN) )
     &      + ( ( ( DLZ(SKEN-1) - DLZ(SKEN) ) / 2.0D0 )
     &          * C2(IJ, SKEN) )
         C0(IJ, SKEN) = ( ( DLZ(SKEN)   * S(IJ, SKEN-1)
     &                    + DLZ(SKEN-1) * S(IJ, SKEN)   )
     &                / ( 2.0D0 * DLZI(SKEN) ) )
     &      - ( ( ( DLZ(SKEN-1) * DLZ(SKEN) ) / 4.0D0 )
     &          * C2(IJ, SKEN) )
      ENDDO

      DO IJ = IJTSTR, IJTEND
      IF ( W(IJ, SKEN+1) >= 0.0D0 ) THEN
         C2(IJ, SKEN+1) = ( 1.0D0 / ( DLZI(SKEN) + DLZI(SKEN+1) ) )
     &      * ( ( ( S(IJ, SKEN-1) - S(IJ, SKEN) ) / DLZI(SKEN)   )
     &        - ( ( S(IJ, SKEN)   - S(IJ, SKEN) ) / DLZI(SKEN+1) ) )
      ELSEIF ( W(IJ, SKEN+1) < 0.0D0 ) THEN
         C2(IJ, SKEN) = ( 1.0D0 / ( DLZI(SKEN+1) + DLZI(SKEN+1) ) )
     &      * ( ( ( S(IJ, SKEN) - 0.0D0 ) / DLZI(SKEN+1) )
     &        - ( ( 0.0D0       - 0.0D0 ) / DLZI(SKEN+1) ) )
      ENDIF
      IF ( W(IJ, SKEN+1) >= 0.0D0 ) THEN
         C1(IJ, SKEN+1) = ( ( S(IJ, SKEN) - S(IJ, SKEN) )
     &                    / DLZI(SKEN+1) )
         C0(IJ, SKEN+1) = ( ( DLZ(SKEN) * S(IJ, SKEN)
     &                      + DLZ(SKEN) * S(IJ, SKEN) )
     &                  / ( 2.0D0 * DLZI(SKEN+1) ) )
     &      - ( ( ( DLZ(SKEN) * DLZ(SKEN) ) / 4.0D0 )
     &          * C2(IJ, SKEN+1) )
      ELSEIF ( W(IJ, SKEN+1) < 0.0D0 ) THEN
         C1(IJ, SKEN+1) = ( ( 0.0D0 - S(IJ, SKEN) )
     &                    / DLZI(SKEN+1) )
         C0(IJ, SKEN+1) = ( ( DLZ(SKEN) * S(IJ, SKEN)
     &                      + DLZ(SKEN) * 0.0D0 )
     &                  / ( 2.0D0 * DLZI(SKEN+1) ) )
     &      - ( ( ( DLZ(SKEN) * DLZ(SKEN) ) / 4.0D0 )
     &          * C2(IJ, SKEN+1) )
      ENDIF
      ENDDO

      DO K = SKST+1, SKEN+1
      DO IJ = IJTSTR, IJTEND
         SI(IJ, K) = C0(IJ, K)
     &   - ( C1(IJ, K) * ( ( W(IJ, K) * TS ) / 2.0D0 ) )
     &   + ( C2(IJ, K) * ( ( ( ( W(IJ, K) * TS ) ** 2.0D0 ) / 3.0D0 )
     &             - ( ( ( DLZI(K) ** 2.0D0 ) / 12.0D0 ) ) ) )
      ENDDO
      ENDDO

!     -- Advection : Quickest Flux limiter --

      IF (SADQF) THEN

      DO IJ = IJTSTR, IJTEND
      IF ( W(IJ, SKST+1) >= 0.0D0 ) THEN
      FLDEL(SKST+1) = S(IJ, SKST+1) - 0.0D0
      FLGAM(SKST+1) = ( ( S(IJ, SKST+1) - S(IJ, SKST) )
     &                  / DLZI(SKST+1) )
     &              - ( ( S(IJ, SKST)   - 0.0D0       )
     &                  / DLZI(SKST)   )
      IF ( ABS( FLGAM(SKST+1)*DLZ(SKST) ) >= ABS( FLDEL(SKST+1) ) ) THEN
         SI(IJ, SKST+1) = S(IJ, SKST)
      ELSE IF ( FLDEL(SKST+1) >= 0 ) THEN
         FLTMP = 0.0D0
     &         + ( S(IJ, SKST) - 0.0D0 ) / W(IJ, SKST+1) * DLZI(SKST+1)
         SI(IJ, SKST+1) =
     &      MIN( SI(IJ, SKST+1), MIN( FLTMP, S(IJ, SKST+1) ) )
         SI(IJ, SKST+1) = MAX( SI(IJ, SKST+1), S(IJ, SKST) )
      ELSE IF ( FLDEL(SKST+1) <  0 ) THEN
         FLTMP = 0.0D0
     &         + ( S(IJ, SKST) - 0.0D0 ) / W(IJ, SKST+1) * DLZI(SKST+1)
         SI(IJ, SKST+1) =
     &      MAX( SI(IJ, SKST+1), MAX( FLTMP, S(IJ, SKST+1) ) )
         SI(IJ, SKST+1) = MIN( SI(IJ, SKST+1), S(IJ, SKST) )
      ENDIF
      ELSEIF ( W(IJ, SKST+1) < 0.0D0 ) THEN
      FLDEL(SKST+1) = S(IJ, SKST) - S(IJ, SKST+2)
      FLGAM(SKST+1) = ( ( S(IJ, SKST)   - S(IJ, SKST+1) )
     &                  / DLZI(SKST+1) )
     &              - ( ( S(IJ, SKST+1) - S(IJ, SKST+2) )
     &                  / DLZI(SKST+2) )
      IF ( ABS( FLGAM(SKST+1)*DLZ(SKST+1) )
     &  >= ABS( FLDEL(SKST+1) ) ) THEN
         SI(IJ, SKST+1) = S(IJ, SKST+1)
      ELSE IF ( FLDEL(SKST+1) >= 0 ) THEN
         FLTMP = S(IJ, SKST+2)
     &- ( S(IJ, SKST+1) - S(IJ, SKST+2) ) / W(IJ, SKST+1) * DLZI(SKST+1)
         SI(IJ, SKST+1) =
     &      MIN( SI(IJ, SKST+1), MIN( FLTMP, S(IJ, SKST) ) )
         SI(IJ, SKST+1) = MAX( SI(IJ, SKST+1), S(IJ, SKST+1) )
      ELSE IF ( FLDEL(SKST+1) <  0 ) THEN
         FLTMP = S(IJ, SKST+2)
     &- ( S(IJ, SKST+1) - S(IJ, SKST+2) ) / W(IJ, SKST+1) * DLZI(SKST+1)
         SI(IJ, SKST+1) =
     &      MAX( SI(IJ, SKST+1), MAX( FLTMP, S(IJ, SKST) ) )
         SI(IJ, SKST+1) = MIN( SI(IJ, SKST+1), S(IJ, SKST+1) )
      ENDIF
      ENDIF
      ENDDO

      DO K = SKST+2, SKEN-1
      DO IJ = IJTSTR, IJTEND
      IF ( W(IJ, K) >= 0.0D0 ) THEN
         FLDEL(K) = S(IJ, K) - S(IJ, K-2)
         FLGAM(K) = ( ( S(IJ, K)   - S(IJ, K-1) ) / DLZI(K)   )
     &            - ( ( S(IJ, K-1) - S(IJ, K-2) ) / DLZI(K-1) )
         IF ( ABS( FLGAM(K)*DLZ(K-1) ) >= ABS( FLDEL(K) ) ) THEN
            SI(IJ, K) = S(IJ, K-1)
         ELSE IF ( FLDEL(K) >= 0 ) THEN
            FLTMP = S(IJ, K-2)
     &            + ( S(IJ, K-1) - S(IJ, K-2) ) / W(IJ, K) * DLZI(K)
            SI(IJ, K) =
     &         MIN( SI(IJ, K), MIN( FLTMP, S(IJ, K) ) )
            SI(IJ, K) = MAX( SI(IJ, K), S(IJ, K-1) )
         ELSE IF ( FLDEL(K) <  0 ) THEN
            FLTMP = S(IJ, K-2)
     &            + ( S(IJ, K-1) - S(IJ, K-2) ) / W(IJ, K) * DLZI(K)
            SI(IJ, K) =
     &         MAX( SI(IJ, K), MAX( FLTMP, S(IJ, K) ) )
            SI(IJ, K) = MIN( SI(IJ, K), S(IJ, K-1) )
         ENDIF
      ELSEIF ( W(IJ, K) < 0.0D0 ) THEN
         FLDEL(K) = S(IJ, K-1) - S(IJ, K+1)
         FLGAM(K) = ( ( S(IJ, K-1) - S(IJ, K)   ) / DLZI(K)   )
     &            - ( ( S(IJ, K)   - S(IJ, K+1) ) / DLZI(K+1) )
         IF ( ABS( FLGAM(K)*DLZ(K) ) >= ABS( FLDEL(K) ) ) THEN
            SI(IJ, K) = S(IJ, K)
         ELSE IF ( FLDEL(K) >= 0 ) THEN
            FLTMP = S(IJ, K+1)
     &            - ( S(IJ, K) - S(IJ, K+1) ) / W(IJ, K) * DLZI(K)
            SI(IJ, K) =
     &         MIN( SI(IJ, K), MIN( FLTMP, S(IJ, K-1) ) )
            SI(IJ, K) = MAX( SI(IJ, K), S(IJ, K) )
         ELSE IF ( FLDEL(K) <  0 ) THEN
            FLTMP = S(IJ, K+1)
     &            - ( S(IJ, K) - S(IJ, K+1) ) / W(IJ, K) * DLZI(K)
            SI(IJ, K) =
     &         MAX( SI(IJ, K), MAX( FLTMP, S(IJ, K-1) ) )
            SI(IJ, K) = MIN( SI(IJ, K), S(IJ, K) )
         ENDIF
      ENDIF
      ENDDO
      ENDDO

      DO IJ = IJTSTR, IJTEND
      IF ( W(IJ, SKEN) >= 0.0D0 ) THEN
      FLDEL(SKEN) = S(IJ, SKEN) - S(IJ, SKEN-2)
      FLGAM(SKEN) = ( ( S(IJ, SKEN)   - S(IJ, SKEN-1) ) / DLZI(SKEN)   )
     &            - ( ( S(IJ, SKEN-1) - S(IJ, SKEN-2) ) / DLZI(SKEN-1) )
      IF ( ABS( FLGAM(SKEN)*DLZ(SKEN-1) ) >= ABS( FLDEL(SKEN) ) ) THEN
         SI(IJ, SKEN) = S(IJ, SKEN-1)
      ELSE IF ( FLDEL(SKEN) >= 0 ) THEN
         FLTMP = S(IJ, SKEN-2)
     &+ ( S(IJ, SKEN-1) - S(IJ, SKEN-2) ) / W(IJ, SKEN) * DLZI(SKEN)
         SI(IJ, SKEN) =
     &      MIN( SI(IJ, SKEN), MIN( FLTMP, S(IJ, SKEN) ) )
         SI(IJ, SKEN) = MAX( SI(IJ, SKEN), S(IJ, SKEN-1) )
      ELSE IF ( FLDEL(SKEN) <  0 ) THEN
         FLTMP = S(IJ, SKEN-2)
     &+ ( S(IJ, SKEN-1) - S(IJ, SKEN-2) ) / W(IJ, SKEN) * DLZI(SKEN)
         SI(IJ, SKEN) =
     &      MAX( SI(IJ, SKEN), MAX( FLTMP, S(IJ, SKEN) ) )
         SI(IJ, SKEN) = MIN( SI(IJ, SKEN), S(IJ, SKEN-1) )
      ENDIF
      ELSEIF ( W(IJ, SKEN) < 0.0D0 ) THEN
      FLDEL(SKEN) = S(IJ, SKEN-1) - 0.0D0
      FLGAM(SKEN) = ( ( S(IJ, SKEN-1) - S(IJ, SKEN) ) / DLZI(SKEN)   )
     &            - ( ( S(IJ, SKEN)   - 0.0D0       ) / DLZI(SKEN+1) )
      IF ( ABS( FLGAM(SKEN)*DLZ(SKEN) ) >= ABS( FLDEL(SKEN) ) ) THEN
         SI(IJ, SKEN) = S(IJ, SKEN)
      ELSE IF ( FLDEL(SKEN) >= 0 ) THEN
         FLTMP = 0.0D0
     &- ( S(IJ, SKEN) - 0.0D0 ) / W(IJ, SKEN) * DLZI(SKEN)
         SI(IJ, SKEN) =
     &      MIN( SI(IJ, SKEN), MIN( FLTMP, S(IJ, SKEN-1) ) )
         SI(IJ, SKEN) = MAX( SI(IJ, SKEN), S(IJ, SKEN) )
      ELSE IF ( FLDEL(SKEN) <  0 ) THEN
         FLTMP = 0.0D0
     &- ( S(IJ, SKEN) - 0.0D0 ) / W(IJ, SKEN) * DLZI(SKEN)
         SI(IJ, SKEN) =
     &      MAX( SI(IJ, SKEN), MAX( FLTMP, S(IJ, SKEN-1) ) )
         SI(IJ, SKEN) = MIN( SI(IJ, SKEN), S(IJ, SKEN) )
      ENDIF
      ENDIF
      ENDDO

      DO IJ = IJTSTR, IJTEND
      IF ( W(IJ, SKEN+1) >= 0.0D0 ) THEN
      FLDEL(SKEN+1) = S(IJ, SKEN) - S(IJ, SKEN-1)
      FLGAM(SKEN+1) = ( ( S(IJ, SKEN) - S(IJ, SKEN)   ) / DLZI(SKEN+1) )
     &              - ( ( S(IJ, SKEN) - S(IJ, SKEN-1) ) / DLZI(SKEN)   )
      IF ( ABS( FLGAM(SKEN+1)*DLZ(SKEN) ) >= ABS( FLDEL(SKEN+1) ) ) THEN
         SI(IJ, SKEN+1) = S(IJ, SKEN)
      ELSE IF ( FLDEL(SKEN+1) >= 0 ) THEN
         FLTMP = S(IJ, SKEN-1)
     &+ ( S(IJ, SKEN) - S(IJ, SKEN-1) ) / W(IJ, SKEN+1) * DLZI(SKEN+1)
         SI(IJ, SKEN+1) =
     &      MIN( SI(IJ, SKEN+1), MIN( FLTMP, S(IJ, SKEN) ) )
         SI(IJ, SKEN+1) = MAX( SI(IJ, SKEN+1), S(IJ, SKEN) )
      ELSE IF ( FLDEL(SKEN+1) <  0 ) THEN
         FLTMP = S(IJ, SKEN-1)
     &+ ( S(IJ, SKEN) - S(IJ, SKEN-1) ) / W(IJ, SKEN+1) * DLZI(SKEN+1)
         SI(IJ, SKEN+1) =
     &      MAX( SI(IJ, SKEN+1), MAX( FLTMP, S(IJ, SKEN) ) )
         SI(IJ, SKEN+1) = MIN( SI(IJ, SKEN+1), S(IJ, SKEN) )
      ENDIF
      ELSEIF ( W(IJ, SKEN+1) < 0.0D0 ) THEN
      FLDEL(SKEN+1) = S(IJ, SKEN) - 0.0D0
      FLGAM(SKEN+1) = ( ( S(IJ, SKEN) - 0.0D0 ) / DLZI(SKEN+1) )
     &              - ( ( 0.0D0       - 0.0D0 ) / DLZI(SKEN+1) )
      IF ( ABS( FLGAM(SKEN+1)*DLZ(SKEN) ) >= ABS( FLDEL(SKEN+1) ) ) THEN
         SI(IJ, SKEN+1) = 0.0D0
      ELSE IF ( FLDEL(SKEN+1) >= 0 ) THEN
         FLTMP = 0.0D0
     &- ( 0.0D0 - 0.0D0 ) / W(IJ, SKEN+1) * DLZI(SKEN+1)
         SI(IJ, SKEN+1) =
     &      MIN( SI(IJ, SKEN+1), MIN( FLTMP, S(IJ, SKEN) ) )
         SI(IJ, SKEN+1) = MAX( SI(IJ, SKEN+1), 0.0D0 )
      ELSE IF ( FLDEL(SKEN+1) <  0 ) THEN
         FLTMP = 0.0D0
     &- ( 0.0D0 - 0.0D0 ) / W(IJ, SKEN+1) * DLZI(SKEN+1)
         SI(IJ, SKEN+1) =
     &      MAX( SI(IJ, SKEN+1), MAX( FLTMP, S(IJ, SKEN) ) )
         SI(IJ, SKEN+1) = MIN( SI(IJ, SKEN+1), 0.0D0 )
      ENDIF
      ENDIF
      ENDDO

      ENDIF

      DO K = SKST+1, SKEN+1
      DO IJ = IJTSTR, IJTEND
         FLX(IJ, K) = ( 1.0D0 - PHII(K) ) * W(IJ, K) * SI(IJ, K)
      ENDDO
      ENDDO

      ENDIF

!     -- Diffusion : Crank nicolson

      IF (SDIC) THEN

      DO K = SKST, SKEN
         EAL(K) = ( ( 1.0D0 - PHI(K) ) / TS )
         EBE(K) = ( ( ( 1.0D0 - PHII(K) )   * DBI(K) )
     &            / ( 2.0D0 * DLZ(K) * DLZI(K)   ) )
         EGA(K) = ( ( ( 1.0D0 - PHII(K+1) ) * DBI(K+1) )
     &            / ( 2.0D0 * DLZ(K) * DLZI(K+1) ) )
      ENDDO

      DO K = SKST+1, SKEN-1
         EA(K)     = - EBE(K)
         EB(K)     = EAL(K) + EBE(K) + EGA(K)
         EC(K)     = - EGA(K)
      ENDDO

         EA(SKST)  = - EBE(SKST)
         EB(SKST)  = EAL(SKST) + EGA(SKST)
         EC(SKST)  = - EGA(SKST)
         EA(SKEN)  = - EBE(SKEN)
         EB(SKEN)  = EAL(SKEN) + EBE(SKEN)
         EC(SKEN)  = - EGA(SKEN)

      DO IJ = IJTSTR, IJTEND
         ED(IJ, SKST) =
     &      + ( EAL(SKST) - EGA(SKST) ) * S(IJ, SKST)
     &      + EGA(SKST) * S(IJ, SKST+1)
     &      + ( ( FLX(IJ, SKST) - FLX(IJ, SKST+1) ) / DLZ(SKST) )
     &      - RST(IJ, SKST)
      ENDDO

      DO K = SKST+1, SKEN-1
      DO IJ = IJTSTR, IJTEND
         ED(IJ, K) =
     &      + EBE(K) * S(IJ, K-1)
     &      + ( EAL(K) - EBE(K) - EGA(K) ) * S(IJ, K)
     &      + EGA(K) * S(IJ, K+1)
     &      + ( ( FLX(IJ, K) - FLX(IJ, K+1) ) / DLZ(K) )
     &      - RST(IJ, K)
      ENDDO
      ENDDO

      DO IJ = IJTSTR, IJTEND
         ED(IJ, SKEN) =
     &      + EBE(SKEN) * S(IJ, SKEN-1)
     &      + ( EAL(SKEN) - EBE(SKEN) ) * S(IJ, SKEN)
     &      + ( ( FLX(IJ, SKEN) - FLX(IJ, SKEN+1) ) / DLZ(SKEN) )
     &      - RST(IJ, SKEN)
      ENDDO

      CALL STRIDAG( EA, EB, EC, ED, EX )

      ENDIF

!     -- Diffusion : Implicit

      IF (SDII) THEN

      DO K = SKST, SKEN
         EAL(K) = ( ( 1.0D0 - PHI(K) ) / TS )
         EBE(K) = ( ( ( 1.0D0 - PHII(K) )   * DBI(K) )
     &            / ( DLZ(K) * DLZI(K)   ) )
         EGA(K) = ( ( ( 1.0D0 - PHII(K+1) ) * DBI(K+1) )
     &            / ( DLZ(K) * DLZI(K+1) ) )
      ENDDO

      DO K = SKST+1, SKEN-1
         EA(K)     = - EBE(K)
         EB(K)     = EAL(K) + EBE(K) + EGA(K)
         EC(K)     = - EGA(K)
      ENDDO

         EA(SKST)  = - EBE(SKST)
         EB(SKST)  = EAL(SKST) + EGA(SKST)
         EC(SKST)  = - EGA(SKST)
         EA(SKEN)  = - EBE(SKEN)
         EB(SKEN)  = EAL(SKEN) + EBE(SKEN)
         EC(SKEN)  = - EGA(SKEN)

      DO K = SKST, SKEN
      DO IJ = IJTSTR, IJTEND
         ED(IJ, K) =
     &      + EAL(K) * S(IJ, K)
     &      + ( ( FLX(IJ, K) - FLX(IJ, K+1) ) / DLZ(K) )
     &      - RST(IJ, K)
      ENDDO
      ENDDO

      CALL STRIDAG( EA, EB, EC, ED, EX )

      ENDIF

!     -- t = t + dt --

      DO K = SKST, SKEN
      DO IJ = IJTSTR, IJTEND
         SN(IJ, K) = MAX( EX(IJ, K), 0.0D0 )
         EX(IJ, K) = 0.0D0
      ENDDO
      ENDDO

      DO IJ = IJTSTR, IJTEND
         BFZ(IJ) = FLX(IJ, SKEN+1)
      ENDDO

      RETURN

! =====================================================================
      ENTRY SOLID_CLAY(
     I      FF, RST, W, PHI, PHII, DBI, S, RHOCLA,
     O      SN, BFZ )

      DO IJ = IJTSTR, IJTEND
         FLX(IJ, SKST) = FF(IJ)
      ENDDO

!     -- Advection : Upwind --

      IF (SADU) THEN

      DO K = SKST+1, SKEN
      DO IJ = IJTSTR, IJTEND
         IF ( W(IJ, K) >= 0.0D0 ) THEN
            SI(IJ, K) = S(IJ, K-1)
         ELSE IF ( W(IJ, K) < 0.0D0 ) THEN
            SI(IJ, K) = S(IJ, K)
         ENDIF
         FLX(IJ, K)   = ( 1.0D0 - PHII(K) ) * W(IJ, K) * SI(IJ, K)
      ENDDO
      ENDDO

      DO IJ = IJTSTR, IJTEND
      IF ( W(IJ, SKEN+1) >= 0.0D0 ) THEN
         SI(IJ, SKEN+1) = S(IJ, SKEN)
      ELSE IF ( W(IJ, SKEN+1) < 0.0D0 ) THEN
         SI(IJ, SKEN+1) = 1.0D+6 * RHOCLA
      ENDIF
         FLX(IJ, SKEN+1) = ( 1.0D0 - PHII(SKEN+1) )
     &                   * W(IJ, SKEN+1) * SI(IJ, SKEN+1)
      ENDDO

      ENDIF

!     -- Advection : Quickest --

      IF (SADQ) THEN

      DO K = SKST+2, SKEN-1
      DO IJ = IJTSTR, IJTEND
         IF ( W(IJ, K) >= 0.0D0 ) THEN
            C2(IJ, K) = ( 1.0D0 / ( DLZI(K-1) + DLZI(K) ) )
     &         * ( ( ( S(IJ, K-2) - S(IJ, K-1) ) / DLZI(K-1) )
     &           - ( ( S(IJ, K-1) - S(IJ, K)   ) / DLZI(K)   ) )
         ELSEIF ( W(IJ, K) < 0.0D0 ) THEN
            C2(IJ, K) = ( 1.0D0 / ( DLZI(K) + DLZI(K+1) ) )
     &         * ( ( ( S(IJ, K-1) - S(IJ, K) ) / DLZI(K)   )
     &           - ( ( S(IJ, K) - S(IJ, K+1) ) / DLZI(K+1) ) )
         ENDIF
            C1(IJ, K) = ( ( S(IJ, K) - S(IJ, K-1) ) / DLZI(K) )
     &         + ( ( ( DLZ(K-1) - DLZ(K) ) / 2.0D0 ) * C2(IJ, K) )
            C0(IJ, K) = ( ( DLZ(K) * S(IJ, K-1) + DLZ(K-1) * S(IJ, K) )
     &                / ( 2.0D0 * DLZI(K) ) )
     &         - ( ( ( DLZ(K-1) * DLZ(K) ) / 4.0D0 ) * C2(IJ, K) )
      ENDDO
      ENDDO

      DO IJ = IJTSTR, IJTEND
      IF ( W(IJ, SKST+1) >= 0.0D0 ) THEN
         C2(IJ, SKST+1) = ( 1.0D0 / ( DLZI(SKST) + DLZI(SKST+1) ) )
     &      * ( ( ( 0.0D0       - S(IJ, SKST)   ) / DLZI(SKST) )
     &        - ( ( S(IJ, SKST) - S(IJ, SKST+1) ) / DLZI(SKST+1) ) )
      ELSEIF ( W(IJ, SKST+1) < 0.0D0 ) THEN
         C2(IJ, SKST+1) = ( 1.0D0 / ( DLZI(SKST+1) + DLZI(SKST+2) ) )
     &      * ( ( ( S(IJ, SKST)   - S(IJ, SKST+1) ) / DLZI(SKST+1) )
     &        - ( ( S(IJ, SKST+1) - S(IJ, SKST+2) ) / DLZI(SKST+2) ) )
      ENDIF
         C1(IJ, SKST+1) = ( ( S(IJ, SKST+1) - S(IJ, SKST) )
     &                      / DLZI(SKST+1) )
     &      + ( ( ( DLZ(SKST) - DLZ(SKST+1) ) / 2.0D0 )
     &          * C2(IJ, SKST+1) )
         C0(IJ, SKST+1) = ( ( DLZ(SKST+1) * S(IJ, SKST)
     &                      + DLZ(SKST)   * S(IJ, SKST+1) )
     &                  / ( 2.0D0 * DLZI(SKST+1) ) )
     &      - ( ( ( DLZ(SKST) * DLZ(SKST+1) ) / 4.0D0 )
     &          * C2(IJ, SKST+1) )
      ENDDO

      DO IJ = IJTSTR, IJTEND
      IF ( W(IJ, SKEN) >= 0.0D0 ) THEN
         C2(IJ, SKEN) = ( 1.0D0 / ( DLZI(SKEN-1) + DLZI(SKEN) ) )
     &      * ( ( ( S(IJ, SKEN-2) - S(IJ, SKEN-1) ) / DLZI(SKEN-1) )
     &        - ( ( S(IJ, SKEN-1) - S(IJ, SKEN)   ) / DLZI(SKEN)   ) )
      ELSEIF ( W(IJ, SKEN) < 0.0D0 ) THEN
         C2(IJ, SKEN) = ( 1.0D0 / ( DLZI(SKEN) + DLZI(SKEN+1) ) )
     &      * ( ( ( S(IJ, SKEN-1) - S(IJ, SKEN) ) / DLZI(SKEN)   )
     &        - ( ( S(IJ, SKEN)   - 1.0D+6 * RHOCLA ) / DLZI(SKEN+1) ) )
      ENDIF
         C1(IJ, SKEN) = ( ( S(IJ, SKEN) - S(IJ, SKEN-1) ) / DLZI(SKEN) )
     &      + ( ( ( DLZ(SKEN-1) - DLZ(SKEN) ) / 2.0D0 )
     &          * C2(IJ, SKEN) )
         C0(IJ, SKEN) = ( ( DLZ(SKEN)   * S(IJ, SKEN-1)
     &                    + DLZ(SKEN-1) * S(IJ, SKEN)   )
     &                / ( 2.0D0 * DLZI(SKEN) ) )
     &      - ( ( ( DLZ(SKEN-1) * DLZ(SKEN) ) / 4.0D0 )
     &          * C2(IJ, SKEN) )
      ENDDO

      DO IJ = IJTSTR, IJTEND
      IF ( W(IJ, SKEN+1) >= 0.0D0 ) THEN
         C2(IJ, SKEN+1) = ( 1.0D0 / ( DLZI(SKEN) + DLZI(SKEN+1) ) )
     &      * ( ( ( S(IJ, SKEN-1) - S(IJ, SKEN) ) / DLZI(SKEN)   )
     &        - ( ( S(IJ, SKEN)   - S(IJ, SKEN) ) / DLZI(SKEN+1) ) )
      ELSEIF ( W(IJ, SKEN+1) < 0.0D0 ) THEN
         C2(IJ, SKEN) = ( 1.0D0 / ( DLZI(SKEN+1) + DLZI(SKEN+1) ) )
     &      * ( ( ( S(IJ, SKEN) - 1.0D+6 * RHOCLA ) / DLZI(SKEN+1) )
     &    - ( ( 1.0D+6 * RHOCLA - 1.0D+6 * RHOCLA ) / DLZI(SKEN+1) ) )
      ENDIF
      IF ( W(IJ, SKEN+1) >= 0.0D0 ) THEN
         C1(IJ, SKEN+1) = ( ( S(IJ, SKEN) - S(IJ, SKEN) )
     &                    / DLZI(SKEN+1) )
         C0(IJ, SKEN+1) = ( ( DLZ(SKEN) * S(IJ, SKEN)
     &                      + DLZ(SKEN) * S(IJ, SKEN) )
     &                  / ( 2.0D0 * DLZI(SKEN+1) ) )
     &      - ( ( ( DLZ(SKEN) * DLZ(SKEN) ) / 4.0D0 )
     &          * C2(IJ, SKEN+1) )
      ELSEIF ( W(IJ, SKEN+1) < 0.0D0 ) THEN
         C1(IJ, SKEN+1) = ( ( 1.0D+6 * RHOCLA - S(IJ, SKEN) )
     &                    / DLZI(SKEN+1) )
         C0(IJ, SKEN+1) = ( ( DLZ(SKEN) * S(IJ, SKEN)
     &                      + DLZ(SKEN) * 1.0D+6 * RHOCLA )
     &                  / ( 2.0D0 * DLZI(SKEN+1) ) )
     &      - ( ( ( DLZ(SKEN) * DLZ(SKEN) ) / 4.0D0 )
     &          * C2(IJ, SKEN+1) )
      ENDIF
      ENDDO

      DO K = SKST+1, SKEN+1
      DO IJ = IJTSTR, IJTEND
         SI(IJ, K) = C0(IJ, K)
     &   - ( C1(IJ, K) * ( ( W(IJ, K) * TS ) / 2.0D0 ) )
     &   + ( C2(IJ, K) * ( ( ( ( W(IJ, K) * TS ) ** 2.0D0 ) / 3.0D0 )
     &             - ( ( ( DLZI(K) ** 2.0D0 ) / 12.0D0 ) ) ) )
      ENDDO
      ENDDO

!     -- Advection : Quickest Flux limiter --

      IF (SADQF) THEN

      DO IJ = IJTSTR, IJTEND
      IF ( W(IJ, SKST+1) >= 0.0D0 ) THEN
      FLDEL(SKST+1) = S(IJ, SKST+1) - 0.0D0
      FLGAM(SKST+1) = ( ( S(IJ, SKST+1) - S(IJ, SKST) )
     &                  / DLZI(SKST+1) )
     &              - ( ( S(IJ, SKST)   - 0.0D0       )
     &                  / DLZI(SKST)   )
      IF ( ABS( FLGAM(SKST+1)*DLZ(SKST) ) >= ABS( FLDEL(SKST+1) ) ) THEN
         SI(IJ, SKST+1) = S(IJ, SKST)
      ELSE IF ( FLDEL(SKST+1) >= 0 ) THEN
         FLTMP = 0.0D0
     &         + ( S(IJ, SKST) - 0.0D0 ) / W(IJ, SKST+1) * DLZI(SKST+1)
         SI(IJ, SKST+1) =
     &      MIN( SI(IJ, SKST+1), MIN( FLTMP, S(IJ, SKST+1) ) )
         SI(IJ, SKST+1) = MAX( SI(IJ, SKST+1), S(IJ, SKST) )
      ELSE IF ( FLDEL(SKST+1) <  0 ) THEN
         FLTMP = 0.0D0
     &         + ( S(IJ, SKST) - 0.0D0 ) / W(IJ, SKST+1) * DLZI(SKST+1)
         SI(IJ, SKST+1) =
     &      MAX( SI(IJ, SKST+1), MAX( FLTMP, S(IJ, SKST+1) ) )
         SI(IJ, SKST+1) = MIN( SI(IJ, SKST+1), S(IJ, SKST) )
      ENDIF
      ELSEIF ( W(IJ, SKST+1) < 0.0D0 ) THEN
      FLDEL(SKST+1) = S(IJ, SKST) - S(IJ, SKST+2)
      FLGAM(SKST+1) = ( ( S(IJ, SKST)   - S(IJ, SKST+1) )
     &                  / DLZI(SKST+1) )
     &              - ( ( S(IJ, SKST+1) - S(IJ, SKST+2) )
     &                  / DLZI(SKST+2) )
      IF ( ABS( FLGAM(SKST+1)*DLZ(SKST+1) )
     &  >= ABS( FLDEL(SKST+1) ) ) THEN
         SI(IJ, SKST+1) = S(IJ, SKST+1)
      ELSE IF ( FLDEL(SKST+1) >= 0 ) THEN
         FLTMP = S(IJ, SKST+2)
     &- ( S(IJ, SKST+1) - S(IJ, SKST+2) ) / W(IJ, SKST+1) * DLZI(SKST+1)
         SI(IJ, SKST+1) =
     &      MIN( SI(IJ, SKST+1), MIN( FLTMP, S(IJ, SKST) ) )
         SI(IJ, SKST+1) = MAX( SI(IJ, SKST+1), S(IJ, SKST+1) )
      ELSE IF ( FLDEL(SKST+1) <  0 ) THEN
         FLTMP = S(IJ, SKST+2)
     &- ( S(IJ, SKST+1) - S(IJ, SKST+2) ) / W(IJ, SKST+1) * DLZI(SKST+1)
         SI(IJ, SKST+1) =
     &      MAX( SI(IJ, SKST+1), MAX( FLTMP, S(IJ, SKST) ) )
         SI(IJ, SKST+1) = MIN( SI(IJ, SKST+1), S(IJ, SKST+1) )
      ENDIF
      ENDIF
      ENDDO

      DO K = SKST+2, SKEN-1
      DO IJ = IJTSTR, IJTEND
      IF ( W(IJ, K) >= 0.0D0 ) THEN
         FLDEL(K) = S(IJ, K) - S(IJ, K-2)
         FLGAM(K) = ( ( S(IJ, K)   - S(IJ, K-1) ) / DLZI(K)   )
     &            - ( ( S(IJ, K-1) - S(IJ, K-2) ) / DLZI(K-1) )
         IF ( ABS( FLGAM(K)*DLZ(K-1) ) >= ABS( FLDEL(K) ) ) THEN
            SI(IJ, K) = S(IJ, K-1)
         ELSE IF ( FLDEL(K) >= 0 ) THEN
            FLTMP = S(IJ, K-2)
     &            + ( S(IJ, K-1) - S(IJ, K-2) ) / W(IJ, K) * DLZI(K)
            SI(IJ, K) =
     &         MIN( SI(IJ, K), MIN( FLTMP, S(IJ, K) ) )
            SI(IJ, K) = MAX( SI(IJ, K), S(IJ, K-1) )
         ELSE IF ( FLDEL(K) <  0 ) THEN
            FLTMP = S(IJ, K-2)
     &            + ( S(IJ, K-1) - S(IJ, K-2) ) / W(IJ, K) * DLZI(K)
            SI(IJ, K) =
     &         MAX( SI(IJ, K), MAX( FLTMP, S(IJ, K) ) )
            SI(IJ, K) = MIN( SI(IJ, K), S(IJ, K-1) )
         ENDIF
      ELSEIF ( W(IJ, K) < 0.0D0 ) THEN
         FLDEL(K) = S(IJ, K-1) - S(IJ, K+1)
         FLGAM(K) = ( ( S(IJ, K-1) - S(IJ, K)   ) / DLZI(K)   )
     &            - ( ( S(IJ, K)   - S(IJ, K+1) ) / DLZI(K+1) )
         IF ( ABS( FLGAM(K)*DLZ(K) ) >= ABS( FLDEL(K) ) ) THEN
            SI(IJ, K) = S(IJ, K)
         ELSE IF ( FLDEL(K) >= 0 ) THEN
            FLTMP = S(IJ, K+1)
     &            - ( S(IJ, K) - S(IJ, K+1) ) / W(IJ, K) * DLZI(K)
            SI(IJ, K) =
     &         MIN( SI(IJ, K), MIN( FLTMP, S(IJ, K-1) ) )
            SI(IJ, K) = MAX( SI(IJ, K), S(IJ, K) )
         ELSE IF ( FLDEL(K) <  0 ) THEN
            FLTMP = S(IJ, K+1)
     &            - ( S(IJ, K) - S(IJ, K+1) ) / W(IJ, K) * DLZI(K)
            SI(IJ, K) =
     &         MAX( SI(IJ, K), MAX( FLTMP, S(IJ, K-1) ) )
            SI(IJ, K) = MIN( SI(IJ, K), S(IJ, K) )
         ENDIF
      ENDIF
      ENDDO
      ENDDO

      DO IJ = IJTSTR, IJTEND
      IF ( W(IJ, SKEN) >= 0.0D0 ) THEN
      FLDEL(SKEN) = S(IJ, SKEN) - S(IJ, SKEN-2)
      FLGAM(SKEN) = ( ( S(IJ, SKEN)   - S(IJ, SKEN-1) ) / DLZI(SKEN)   )
     &            - ( ( S(IJ, SKEN-1) - S(IJ, SKEN-2) ) / DLZI(SKEN-1) )
      IF ( ABS( FLGAM(SKEN)*DLZ(SKEN-1) ) >= ABS( FLDEL(SKEN) ) ) THEN
         SI(IJ, SKEN) = S(IJ, SKEN-1)
      ELSE IF ( FLDEL(SKEN) >= 0 ) THEN
         FLTMP = S(IJ, SKEN-2)
     &+ ( S(IJ, SKEN-1) - S(IJ, SKEN-2) ) / W(IJ, SKEN) * DLZI(SKEN)
         SI(IJ, SKEN) =
     &      MIN( SI(IJ, SKEN), MIN( FLTMP, S(IJ, SKEN) ) )
         SI(IJ, SKEN) = MAX( SI(IJ, SKEN), S(IJ, SKEN-1) )
      ELSE IF ( FLDEL(SKEN) <  0 ) THEN
         FLTMP = S(IJ, SKEN-2)
     &+ ( S(IJ, SKEN-1) - S(IJ, SKEN-2) ) / W(IJ, SKEN) * DLZI(SKEN)
         SI(IJ, SKEN) =
     &      MAX( SI(IJ, SKEN), MAX( FLTMP, S(IJ, SKEN) ) )
         SI(IJ, SKEN) = MIN( SI(IJ, SKEN), S(IJ, SKEN-1) )
      ENDIF
      ELSEIF ( W(IJ, SKEN) < 0.0D0 ) THEN
      FLDEL(SKEN) = S(IJ, SKEN-1) - 1.0D+6 * RHOCLA
      FLGAM(SKEN) = ( ( S(IJ, SKEN-1) - S(IJ, SKEN) ) / DLZI(SKEN)   )
     &            - ( ( S(IJ, SKEN) - 1.0D+6 * RHOCLA ) / DLZI(SKEN+1) )
      IF ( ABS( FLGAM(SKEN)*DLZ(SKEN) ) >= ABS( FLDEL(SKEN) ) ) THEN
         SI(IJ, SKEN) = S(IJ, SKEN)
      ELSE IF ( FLDEL(SKEN) >= 0 ) THEN
         FLTMP = 1.0D+6 * RHOCLA
     &- ( S(IJ, SKEN) - 1.0D+6 * RHOCLA ) / W(IJ, SKEN) * DLZI(SKEN)
         SI(IJ, SKEN) =
     &      MIN( SI(IJ, SKEN), MIN( FLTMP, S(IJ, SKEN-1) ) )
         SI(IJ, SKEN) = MAX( SI(IJ, SKEN), S(IJ, SKEN) )
      ELSE IF ( FLDEL(SKEN) <  0 ) THEN
         FLTMP = 1.0D+6 * RHOCLA
     &- ( S(IJ, SKEN) - 1.0D+6 * RHOCLA ) / W(IJ, SKEN) * DLZI(SKEN)
         SI(IJ, SKEN) =
     &      MAX( SI(IJ, SKEN), MAX( FLTMP, S(IJ, SKEN-1) ) )
         SI(IJ, SKEN) = MIN( SI(IJ, SKEN), S(IJ, SKEN) )
      ENDIF
      ENDIF
      ENDDO

      DO IJ = IJTSTR, IJTEND
      IF ( W(IJ, SKEN+1) >= 0.0D0 ) THEN
      FLDEL(SKEN+1) = S(IJ, SKEN) - S(IJ, SKEN-1)
      FLGAM(SKEN+1) = ( ( S(IJ, SKEN) - S(IJ, SKEN)   ) / DLZI(SKEN+1) )
     &              - ( ( S(IJ, SKEN) - S(IJ, SKEN-1) ) / DLZI(SKEN)   )
      IF ( ABS( FLGAM(SKEN+1)*DLZ(SKEN) ) >= ABS( FLDEL(SKEN+1) ) ) THEN
         SI(IJ, SKEN+1) = S(IJ, SKEN)
      ELSE IF ( FLDEL(SKEN+1) >= 0 ) THEN
         FLTMP = S(IJ, SKEN-1)
     &+ ( S(IJ, SKEN) - S(IJ, SKEN-1) ) / W(IJ, SKEN+1) * DLZI(SKEN+1)
         SI(IJ, SKEN+1) =
     &      MIN( SI(IJ, SKEN+1), MIN( FLTMP, S(IJ, SKEN) ) )
         SI(IJ, SKEN+1) = MAX( SI(IJ, SKEN+1), S(IJ, SKEN) )
      ELSE IF ( FLDEL(SKEN+1) <  0 ) THEN
         FLTMP = S(IJ, SKEN-1)
     &+ ( S(IJ, SKEN) - S(IJ, SKEN-1) ) / W(IJ, SKEN+1) * DLZI(SKEN+1)
         SI(IJ, SKEN+1) =
     &      MAX( SI(IJ, SKEN+1), MAX( FLTMP, S(IJ, SKEN) ) )
         SI(IJ, SKEN+1) = MIN( SI(IJ, SKEN+1), S(IJ, SKEN) )
      ENDIF
      ELSEIF ( W(IJ, SKEN+1) < 0.0D0 ) THEN
      FLDEL(SKEN+1) = S(IJ, SKEN) - 1.0D+6 * RHOCLA
      FLGAM(SKEN+1)
     &= ( ( S(IJ, SKEN)     - 1.0D+6 * RHOCLA ) / DLZI(SKEN+1) )
     &- ( ( 1.0D+6 * RHOCLA - 1.0D+6 * RHOCLA ) / DLZI(SKEN+1) )
      IF ( ABS( FLGAM(SKEN+1)*DLZ(SKEN) ) >= ABS( FLDEL(SKEN+1) ) ) THEN
         SI(IJ, SKEN+1) = 1.0D+6 * RHOCLA
      ELSE IF ( FLDEL(SKEN+1) >= 0 ) THEN
         FLTMP = 1.0D+6 * RHOCLA
     &         - ( 1.0D+6 * RHOCLA - 1.0D+6 * RHOCLA )
     &         / W(IJ, SKEN+1) * DLZI(SKEN+1)
         SI(IJ, SKEN+1) =
     &      MIN( SI(IJ, SKEN+1), MIN( FLTMP, S(IJ, SKEN) ) )
         SI(IJ, SKEN+1) = MAX( SI(IJ, SKEN+1), 1.0D+6 * RHOCLA )
      ELSE IF ( FLDEL(SKEN+1) <  0 ) THEN
         FLTMP = 1.0D+6 * RHOCLA
     &         - ( 1.0D+6 * RHOCLA - 1.0D+6 * RHOCLA )
     &         / W(IJ, SKEN+1) * DLZI(SKEN+1)
         SI(IJ, SKEN+1) =
     &      MAX( SI(IJ, SKEN+1), MAX( FLTMP, S(IJ, SKEN) ) )
         SI(IJ, SKEN+1) = MIN( SI(IJ, SKEN+1), 1.0D+6 * RHOCLA )
      ENDIF
      ENDIF
      ENDDO

      ENDIF

      DO K = SKST+1, SKEN+1
      DO IJ = IJTSTR, IJTEND
         FLX(IJ, K) = ( 1.0D0 - PHII(K) ) * W(IJ, K) * SI(IJ, K)
      ENDDO
      ENDDO

      ENDIF

!     -- Diffusion : Crank nicolson

      IF (SDIC) THEN

      DO K = SKST, SKEN
         EAL(K) = ( ( 1.0D0 - PHI(K) ) / TS )
         EBE(K) = ( ( ( 1.0D0 - PHII(K) )   * DBI(K) )
     &            / ( 2.0D0 * DLZ(K) * DLZI(K)   ) )
         EGA(K) = ( ( ( 1.0D0 - PHII(K+1) ) * DBI(K+1) )
     &            / ( 2.0D0 * DLZ(K) * DLZI(K+1) ) )
      ENDDO

      DO K = SKST+1, SKEN-1
         EA(K)     = - EBE(K)
         EB(K)     = EAL(K) + EBE(K) + EGA(K)
         EC(K)     = - EGA(K)
      ENDDO

         EA(SKST)  = - EBE(SKST)
         EB(SKST)  = EAL(SKST) + EGA(SKST)
         EC(SKST)  = - EGA(SKST)
         EA(SKEN)  = - EBE(SKEN)
         EB(SKEN)  = EAL(SKEN) + EBE(SKEN)
         EC(SKEN)  = - EGA(SKEN)

      DO IJ = IJTSTR, IJTEND
         ED(IJ, SKST) =
     &      + ( EAL(SKST) - EGA(SKST) ) * S(IJ, SKST)
     &      + EGA(SKST) * S(IJ, SKST+1)
     &      + ( ( FLX(IJ, SKST) - FLX(IJ, SKST+1) ) / DLZ(SKST) )
     &      - RST(IJ, SKST)
      ENDDO

      DO K = SKST+1, SKEN-1
      DO IJ = IJTSTR, IJTEND
         ED(IJ, K) =
     &      + EBE(K) * S(IJ, K-1)
     &      + ( EAL(K) - EBE(K) - EGA(K) ) * S(IJ, K)
     &      + EGA(K) * S(IJ, K+1)
     &      + ( ( FLX(IJ, K) - FLX(IJ, K+1) ) / DLZ(K) )
     &      - RST(IJ, K)
      ENDDO
      ENDDO

      DO IJ = IJTSTR, IJTEND
         ED(IJ, SKEN) =
     &      + EBE(SKEN) * S(IJ, SKEN-1)
     &      + ( EAL(SKEN) - EBE(SKEN) ) * S(IJ, SKEN)
     &      + ( ( FLX(IJ, SKEN) - FLX(IJ, SKEN+1) ) / DLZ(SKEN) )
     &      - RST(IJ, SKEN)
      ENDDO

      CALL STRIDAG( EA, EB, EC, ED, EX )

      ENDIF

!     -- Diffusion : Implicit

      IF (SDII) THEN

      DO K = SKST, SKEN
         EAL(K) = ( ( 1.0D0 - PHI(K) ) / TS )
         EBE(K) = ( ( ( 1.0D0 - PHII(K) )   * DBI(K) )
     &            / ( DLZ(K) * DLZI(K)   ) )
         EGA(K) = ( ( ( 1.0D0 - PHII(K+1) ) * DBI(K+1) )
     &            / ( DLZ(K) * DLZI(K+1) ) )
      ENDDO

      DO K = SKST+1, SKEN-1
         EA(K)     = - EBE(K)
         EB(K)     = EAL(K) + EBE(K) + EGA(K)
         EC(K)     = - EGA(K)
      ENDDO

         EA(SKST)  = - EBE(SKST)
         EB(SKST)  = EAL(SKST) + EGA(SKST)
         EC(SKST)  = - EGA(SKST)
         EA(SKEN)  = - EBE(SKEN)
         EB(SKEN)  = EAL(SKEN) + EBE(SKEN)
         EC(SKEN)  = - EGA(SKEN)

      DO K = SKST, SKEN
      DO IJ = IJTSTR, IJTEND
         ED(IJ, K) =
     &      + EAL(K) * S(IJ, K)
     &      + ( ( FLX(IJ, K) - FLX(IJ, K+1) ) / DLZ(K) )
     &      - RST(IJ, K)
      ENDDO
      ENDDO

      CALL STRIDAG( EA, EB, EC, ED, EX )

      ENDIF

!     -- t = t + dt --

      DO K = SKST, SKEN
      DO IJ = IJTSTR, IJTEND
         SN(IJ, K) = MAX( EX(IJ, K), 0.0D0 )
         EX(IJ, K) = 0.0D0
      ENDDO
      ENDDO

      DO IJ = IJTSTR, IJTEND
         BFZ(IJ) = FLX(IJ, SKEN+1)
      ENDDO

      END SUBROUTINE SOLID

! =====================================================================
      SUBROUTINE AGE(
     I      FF, W, PHI, PHII, DBI, S,
     O      SN )

      DO IJ = IJTSTR, IJTEND
         FLX(IJ, SKST) = FF(IJ)
      ENDDO

!     -- Advection : Upwind --

      IF (SADU) THEN

      DO K = SKST+1, SKEN
      DO IJ = IJTSTR, IJTEND
         IF ( W(IJ, K) >= 0.0D0 ) THEN
            SI(IJ, K) = S(IJ, K-1)
         ELSE IF ( W(IJ, K) < 0.0D0 ) THEN
            SI(IJ, K) = S(IJ, K)
         ENDIF
         FLX(IJ, K)   = ( 1.0D0 - PHII(K) ) * W(IJ, K) * SI(IJ, K)
      ENDDO
      ENDDO

      DO IJ = IJTSTR, IJTEND
      IF ( W(IJ, SKEN+1) >= 0.0D0 ) THEN
         SI(IJ, SKEN+1) = S(IJ, SKEN)
      ELSE IF ( W(IJ, SKEN+1) < 0.0D0 ) THEN
         SI(IJ, SKEN+1) = 0.0D0
      ENDIF
         FLX(IJ, SKEN+1) = ( 1.0D0 - PHII(SKEN+1) )
     &                   * W(IJ, SKEN+1) * SI(IJ, SKEN+1)
      ENDDO

      ENDIF

!     -- Advection : Quickest --

      IF (SADQ) THEN

      DO K = SKST+2, SKEN-1
      DO IJ = IJTSTR, IJTEND
         IF ( W(IJ, K) >= 0.0D0 ) THEN
            C2(IJ, K) = ( 1.0D0 / ( DLZI(K-1) + DLZI(K) ) )
     &         * ( ( ( S(IJ, K-2) - S(IJ, K-1) ) / DLZI(K-1) )
     &           - ( ( S(IJ, K-1) - S(IJ, K)   ) / DLZI(K)   ) )
         ELSEIF ( W(IJ, K) < 0.0D0 ) THEN
            C2(IJ, K) = ( 1.0D0 / ( DLZI(K) + DLZI(K+1) ) )
     &         * ( ( ( S(IJ, K-1) - S(IJ, K) ) / DLZI(K)   )
     &           - ( ( S(IJ, K) - S(IJ, K+1) ) / DLZI(K+1) ) )
         ENDIF
            C1(IJ, K) = ( ( S(IJ, K) - S(IJ, K-1) ) / DLZI(K) )
     &         + ( ( ( DLZ(K-1) - DLZ(K) ) / 2.0D0 ) * C2(IJ, K) )
            C0(IJ, K) = ( ( DLZ(K) * S(IJ, K-1) + DLZ(K-1) * S(IJ, K) )
     &                / ( 2.0D0 * DLZI(K) ) )
     &         - ( ( ( DLZ(K-1) * DLZ(K) ) / 4.0D0 ) * C2(IJ, K) )
      ENDDO
      ENDDO

      DO IJ = IJTSTR, IJTEND
      IF ( W(IJ, SKST+1) >= 0.0D0 ) THEN
         C2(IJ, SKST+1) = ( 1.0D0 / ( DLZI(SKST) + DLZI(SKST+1) ) )
     &      * ( ( ( 0.0D0       - S(IJ, SKST)   ) / DLZI(SKST) )
     &        - ( ( S(IJ, SKST) - S(IJ, SKST+1) ) / DLZI(SKST+1) ) )
      ELSEIF ( W(IJ, SKST+1) < 0.0D0 ) THEN
         C2(IJ, SKST+1) = ( 1.0D0 / ( DLZI(SKST+1) + DLZI(SKST+2) ) )
     &      * ( ( ( S(IJ, SKST)   - S(IJ, SKST+1) ) / DLZI(SKST+1) )
     &        - ( ( S(IJ, SKST+1) - S(IJ, SKST+2) ) / DLZI(SKST+2) ) )
      ENDIF
         C1(IJ, SKST+1) = ( ( S(IJ, SKST+1) - S(IJ, SKST) )
     &                      / DLZI(SKST+1) )
     &      + ( ( ( DLZ(SKST) - DLZ(SKST+1) ) / 2.0D0 )
     &          * C2(IJ, SKST+1) )
         C0(IJ, SKST+1) = ( ( DLZ(SKST+1) * S(IJ, SKST)
     &                      + DLZ(SKST)   * S(IJ, SKST+1) )
     &                  / ( 2.0D0 * DLZI(SKST+1) ) )
     &      - ( ( ( DLZ(SKST) * DLZ(SKST+1) ) / 4.0D0 )
     &          * C2(IJ, SKST+1) )
      ENDDO

      DO IJ = IJTSTR, IJTEND
      IF ( W(IJ, SKEN) >= 0.0D0 ) THEN
         C2(IJ, SKEN) = ( 1.0D0 / ( DLZI(SKEN-1) + DLZI(SKEN) ) )
     &      * ( ( ( S(IJ, SKEN-2) - S(IJ, SKEN-1) ) / DLZI(SKEN-1) )
     &        - ( ( S(IJ, SKEN-1) - S(IJ, SKEN)   ) / DLZI(SKEN)   ) )
      ELSEIF ( W(IJ, SKEN) < 0.0D0 ) THEN
         C2(IJ, SKEN) = ( 1.0D0 / ( DLZI(SKEN) + DLZI(SKEN+1) ) )
     &      * ( ( ( S(IJ, SKEN-1) - S(IJ, SKEN) ) / DLZI(SKEN)   )
     &        - ( ( S(IJ, SKEN)   - 0.0D0       ) / DLZI(SKEN+1) ) )
      ENDIF
         C1(IJ, SKEN) = ( ( S(IJ, SKEN) - S(IJ, SKEN-1) ) / DLZI(SKEN) )
     &      + ( ( ( DLZ(SKEN-1) - DLZ(SKEN) ) / 2.0D0 )
     &          * C2(IJ, SKEN) )
         C0(IJ, SKEN) = ( ( DLZ(SKEN)   * S(IJ, SKEN-1)
     &                    + DLZ(SKEN-1) * S(IJ, SKEN)   )
     &                / ( 2.0D0 * DLZI(SKEN) ) )
     &      - ( ( ( DLZ(SKEN-1) * DLZ(SKEN) ) / 4.0D0 )
     &          * C2(IJ, SKEN) )
      ENDDO

      DO IJ = IJTSTR, IJTEND
      IF ( W(IJ, SKEN+1) >= 0.0D0 ) THEN
         C2(IJ, SKEN+1) = ( 1.0D0 / ( DLZI(SKEN) + DLZI(SKEN+1) ) )
     &      * ( ( ( S(IJ, SKEN-1) - S(IJ, SKEN) ) / DLZI(SKEN)   )
     &        - ( ( S(IJ, SKEN)   - S(IJ, SKEN) ) / DLZI(SKEN+1) ) )
      ELSEIF ( W(IJ, SKEN+1) < 0.0D0 ) THEN
         C2(IJ, SKEN) = ( 1.0D0 / ( DLZI(SKEN+1) + DLZI(SKEN+1) ) )
     &      * ( ( ( S(IJ, SKEN) - 0.0D0 ) / DLZI(SKEN+1) )
     &        - ( ( 0.0D0       - 0.0D0 ) / DLZI(SKEN+1) ) )
      ENDIF
      IF ( W(IJ, SKEN+1) >= 0.0D0 ) THEN
         C1(IJ, SKEN+1) = ( ( S(IJ, SKEN) - S(IJ, SKEN) )
     &                    / DLZI(SKEN+1) )
         C0(IJ, SKEN+1) = ( ( DLZ(SKEN) * S(IJ, SKEN)
     &                      + DLZ(SKEN) * S(IJ, SKEN) )
     &                  / ( 2.0D0 * DLZI(SKEN+1) ) )
     &      - ( ( ( DLZ(SKEN) * DLZ(SKEN) ) / 4.0D0 )
     &          * C2(IJ, SKEN+1) )
      ELSEIF ( W(IJ, SKEN+1) < 0.0D0 ) THEN
         C1(IJ, SKEN+1) = ( ( 0.0D0 - S(IJ, SKEN) )
     &                    / DLZI(SKEN+1) )
         C0(IJ, SKEN+1) = ( ( DLZ(SKEN) * S(IJ, SKEN)
     &                      + DLZ(SKEN) * 0.0D0 )
     &                  / ( 2.0D0 * DLZI(SKEN+1) ) )
     &      - ( ( ( DLZ(SKEN) * DLZ(SKEN) ) / 4.0D0 )
     &          * C2(IJ, SKEN+1) )
      ENDIF
      ENDDO

      DO K = SKST+1, SKEN+1
      DO IJ = IJTSTR, IJTEND
         SI(IJ, K) = C0(IJ, K)
     &   - ( C1(IJ, K) * ( ( W(IJ, K) * TS ) / 2.0D0 ) )
     &   + ( C2(IJ, K) * ( ( ( ( W(IJ, K) * TS ) ** 2.0D0 ) / 3.0D0 )
     &             - ( ( ( DLZI(K) ** 2.0D0 ) / 12.0D0 ) ) ) )
      ENDDO
      ENDDO

!     -- Advection : Quickest Flux limiter --

      IF (SADQF) THEN

      DO IJ = IJTSTR, IJTEND
      IF ( W(IJ, SKST+1) >= 0.0D0 ) THEN
      FLDEL(SKST+1) = S(IJ, SKST+1) - 0.0D0
      FLGAM(SKST+1) = ( ( S(IJ, SKST+1) - S(IJ, SKST) )
     &                  / DLZI(SKST+1) )
     &              - ( ( S(IJ, SKST)   - 0.0D0       )
     &                  / DLZI(SKST)   )
      IF ( ABS( FLGAM(SKST+1)*DLZ(SKST) ) >= ABS( FLDEL(SKST+1) ) ) THEN
         SI(IJ, SKST+1) = S(IJ, SKST)
      ELSE IF ( FLDEL(SKST+1) >= 0 ) THEN
         FLTMP = 0.0D0
     &         + ( S(IJ, SKST) - 0.0D0 ) / W(IJ, SKST+1) * DLZI(SKST+1)
         SI(IJ, SKST+1) =
     &      MIN( SI(IJ, SKST+1), MIN( FLTMP, S(IJ, SKST+1) ) )
         SI(IJ, SKST+1) = MAX( SI(IJ, SKST+1), S(IJ, SKST) )
      ELSE IF ( FLDEL(SKST+1) <  0 ) THEN
         FLTMP = 0.0D0
     &         + ( S(IJ, SKST) - 0.0D0 ) / W(IJ, SKST+1) * DLZI(SKST+1)
         SI(IJ, SKST+1) =
     &      MAX( SI(IJ, SKST+1), MAX( FLTMP, S(IJ, SKST+1) ) )
         SI(IJ, SKST+1) = MIN( SI(IJ, SKST+1), S(IJ, SKST) )
      ENDIF
      ELSEIF ( W(IJ, SKST+1) < 0.0D0 ) THEN
      FLDEL(SKST+1) = S(IJ, SKST) - S(IJ, SKST+2)
      FLGAM(SKST+1) = ( ( S(IJ, SKST)   - S(IJ, SKST+1) )
     &                  / DLZI(SKST+1) )
     &              - ( ( S(IJ, SKST+1) - S(IJ, SKST+2) )
     &                  / DLZI(SKST+2) )
      IF ( ABS( FLGAM(SKST+1)*DLZ(SKST+1) )
     &  >= ABS( FLDEL(SKST+1) ) ) THEN
         SI(IJ, SKST+1) = S(IJ, SKST+1)
      ELSE IF ( FLDEL(SKST+1) >= 0 ) THEN
         FLTMP = S(IJ, SKST+2)
     &- ( S(IJ, SKST+1) - S(IJ, SKST+2) ) / W(IJ, SKST+1) * DLZI(SKST+1)
         SI(IJ, SKST+1) =
     &      MIN( SI(IJ, SKST+1), MIN( FLTMP, S(IJ, SKST) ) )
         SI(IJ, SKST+1) = MAX( SI(IJ, SKST+1), S(IJ, SKST+1) )
      ELSE IF ( FLDEL(SKST+1) <  0 ) THEN
         FLTMP = S(IJ, SKST+2)
     &- ( S(IJ, SKST+1) - S(IJ, SKST+2) ) / W(IJ, SKST+1) * DLZI(SKST+1)
         SI(IJ, SKST+1) =
     &      MAX( SI(IJ, SKST+1), MAX( FLTMP, S(IJ, SKST) ) )
         SI(IJ, SKST+1) = MIN( SI(IJ, SKST+1), S(IJ, SKST+1) )
      ENDIF
      ENDIF
      ENDDO

      DO K = SKST+2, SKEN-1
      DO IJ = IJTSTR, IJTEND
      IF ( W(IJ, K) >= 0.0D0 ) THEN
         FLDEL(K) = S(IJ, K) - S(IJ, K-2)
         FLGAM(K) = ( ( S(IJ, K)   - S(IJ, K-1) ) / DLZI(K)   )
     &            - ( ( S(IJ, K-1) - S(IJ, K-2) ) / DLZI(K-1) )
         IF ( ABS( FLGAM(K)*DLZ(K-1) ) >= ABS( FLDEL(K) ) ) THEN
            SI(IJ, K) = S(IJ, K-1)
         ELSE IF ( FLDEL(K) >= 0 ) THEN
            FLTMP = S(IJ, K-2)
     &            + ( S(IJ, K-1) - S(IJ, K-2) ) / W(IJ, K) * DLZI(K)
            SI(IJ, K) =
     &         MIN( SI(IJ, K), MIN( FLTMP, S(IJ, K) ) )
            SI(IJ, K) = MAX( SI(IJ, K), S(IJ, K-1) )
         ELSE IF ( FLDEL(K) <  0 ) THEN
            FLTMP = S(IJ, K-2)
     &            + ( S(IJ, K-1) - S(IJ, K-2) ) / W(IJ, K) * DLZI(K)
            SI(IJ, K) =
     &         MAX( SI(IJ, K), MAX( FLTMP, S(IJ, K) ) )
            SI(IJ, K) = MIN( SI(IJ, K), S(IJ, K-1) )
         ENDIF
      ELSEIF ( W(IJ, K) < 0.0D0 ) THEN
         FLDEL(K) = S(IJ, K-1) - S(IJ, K+1)
         FLGAM(K) = ( ( S(IJ, K-1) - S(IJ, K)   ) / DLZI(K)   )
     &            - ( ( S(IJ, K)   - S(IJ, K+1) ) / DLZI(K+1) )
         IF ( ABS( FLGAM(K)*DLZ(K) ) >= ABS( FLDEL(K) ) ) THEN
            SI(IJ, K) = S(IJ, K)
         ELSE IF ( FLDEL(K) >= 0 ) THEN
            FLTMP = S(IJ, K+1)
     &            - ( S(IJ, K) - S(IJ, K+1) ) / W(IJ, K) * DLZI(K)
            SI(IJ, K) =
     &         MIN( SI(IJ, K), MIN( FLTMP, S(IJ, K-1) ) )
            SI(IJ, K) = MAX( SI(IJ, K), S(IJ, K) )
         ELSE IF ( FLDEL(K) <  0 ) THEN
            FLTMP = S(IJ, K+1)
     &            - ( S(IJ, K) - S(IJ, K+1) ) / W(IJ, K) * DLZI(K)
            SI(IJ, K) =
     &         MAX( SI(IJ, K), MAX( FLTMP, S(IJ, K-1) ) )
            SI(IJ, K) = MIN( SI(IJ, K), S(IJ, K) )
         ENDIF
      ENDIF
      ENDDO
      ENDDO

      DO IJ = IJTSTR, IJTEND
      IF ( W(IJ, SKEN) >= 0.0D0 ) THEN
      FLDEL(SKEN) = S(IJ, SKEN) - S(IJ, SKEN-2)
      FLGAM(SKEN) = ( ( S(IJ, SKEN)   - S(IJ, SKEN-1) ) / DLZI(SKEN)   )
     &            - ( ( S(IJ, SKEN-1) - S(IJ, SKEN-2) ) / DLZI(SKEN-1) )
      IF ( ABS( FLGAM(SKEN)*DLZ(SKEN-1) ) >= ABS( FLDEL(SKEN) ) ) THEN
         SI(IJ, SKEN) = S(IJ, SKEN-1)
      ELSE IF ( FLDEL(SKEN) >= 0 ) THEN
         FLTMP = S(IJ, SKEN-2)
     &+ ( S(IJ, SKEN-1) - S(IJ, SKEN-2) ) / W(IJ, SKEN) * DLZI(SKEN)
         SI(IJ, SKEN) =
     &      MIN( SI(IJ, SKEN), MIN( FLTMP, S(IJ, SKEN) ) )
         SI(IJ, SKEN) = MAX( SI(IJ, SKEN), S(IJ, SKEN-1) )
      ELSE IF ( FLDEL(SKEN) <  0 ) THEN
         FLTMP = S(IJ, SKEN-2)
     &+ ( S(IJ, SKEN-1) - S(IJ, SKEN-2) ) / W(IJ, SKEN) * DLZI(SKEN)
         SI(IJ, SKEN) =
     &      MAX( SI(IJ, SKEN), MAX( FLTMP, S(IJ, SKEN) ) )
         SI(IJ, SKEN) = MIN( SI(IJ, SKEN), S(IJ, SKEN-1) )
      ENDIF
      ELSEIF ( W(IJ, SKEN) < 0.0D0 ) THEN
      FLDEL(SKEN) = S(IJ, SKEN-1) - 0.0D0
      FLGAM(SKEN) = ( ( S(IJ, SKEN-1) - S(IJ, SKEN) ) / DLZI(SKEN)   )
     &            - ( ( S(IJ, SKEN)   - 0.0D0       ) / DLZI(SKEN+1) )
      IF ( ABS( FLGAM(SKEN)*DLZ(SKEN) ) >= ABS( FLDEL(SKEN) ) ) THEN
         SI(IJ, SKEN) = S(IJ, SKEN)
      ELSE IF ( FLDEL(SKEN) >= 0 ) THEN
         FLTMP = 0.0D0
     &- ( S(IJ, SKEN) - 0.0D0 ) / W(IJ, SKEN) * DLZI(SKEN)
         SI(IJ, SKEN) =
     &      MIN( SI(IJ, SKEN), MIN( FLTMP, S(IJ, SKEN-1) ) )
         SI(IJ, SKEN) = MAX( SI(IJ, SKEN), S(IJ, SKEN) )
      ELSE IF ( FLDEL(SKEN) <  0 ) THEN
         FLTMP = 0.0D0
     &- ( S(IJ, SKEN) - 0.0D0 ) / W(IJ, SKEN) * DLZI(SKEN)
         SI(IJ, SKEN) =
     &      MAX( SI(IJ, SKEN), MAX( FLTMP, S(IJ, SKEN-1) ) )
         SI(IJ, SKEN) = MIN( SI(IJ, SKEN), S(IJ, SKEN) )
      ENDIF
      ENDIF
      ENDDO

      DO IJ = IJTSTR, IJTEND
      IF ( W(IJ, SKEN+1) >= 0.0D0 ) THEN
      FLDEL(SKEN+1) = S(IJ, SKEN) - S(IJ, SKEN-1)
      FLGAM(SKEN+1) = ( ( S(IJ, SKEN) - S(IJ, SKEN)   ) / DLZI(SKEN+1) )
     &              - ( ( S(IJ, SKEN) - S(IJ, SKEN-1) ) / DLZI(SKEN)   )
      IF ( ABS( FLGAM(SKEN+1)*DLZ(SKEN) ) >= ABS( FLDEL(SKEN+1) ) ) THEN
         SI(IJ, SKEN+1) = S(IJ, SKEN)
      ELSE IF ( FLDEL(SKEN+1) >= 0 ) THEN
         FLTMP = S(IJ, SKEN-1)
     &+ ( S(IJ, SKEN) - S(IJ, SKEN-1) ) / W(IJ, SKEN+1) * DLZI(SKEN+1)
         SI(IJ, SKEN+1) =
     &      MIN( SI(IJ, SKEN+1), MIN( FLTMP, S(IJ, SKEN) ) )
         SI(IJ, SKEN+1) = MAX( SI(IJ, SKEN+1), S(IJ, SKEN) )
      ELSE IF ( FLDEL(SKEN+1) <  0 ) THEN
         FLTMP = S(IJ, SKEN-1)
     &+ ( S(IJ, SKEN) - S(IJ, SKEN-1) ) / W(IJ, SKEN+1) * DLZI(SKEN+1)
         SI(IJ, SKEN+1) =
     &      MAX( SI(IJ, SKEN+1), MAX( FLTMP, S(IJ, SKEN) ) )
         SI(IJ, SKEN+1) = MIN( SI(IJ, SKEN+1), S(IJ, SKEN) )
      ENDIF
      ELSEIF ( W(IJ, SKEN+1) < 0.0D0 ) THEN
      FLDEL(SKEN+1) = S(IJ, SKEN) - 0.0D0
      FLGAM(SKEN+1) = ( ( S(IJ, SKEN) - 0.0D0 ) / DLZI(SKEN+1) )
     &              - ( ( 0.0D0       - 0.0D0 ) / DLZI(SKEN+1) )
      IF ( ABS( FLGAM(SKEN+1)*DLZ(SKEN) ) >= ABS( FLDEL(SKEN+1) ) ) THEN
         SI(IJ, SKEN+1) = 0.0D0
      ELSE IF ( FLDEL(SKEN+1) >= 0 ) THEN
         FLTMP = 0.0D0
     &- ( 0.0D0 - 0.0D0 ) / W(IJ, SKEN+1) * DLZI(SKEN+1)
         SI(IJ, SKEN+1) =
     &      MIN( SI(IJ, SKEN+1), MIN( FLTMP, S(IJ, SKEN) ) )
         SI(IJ, SKEN+1) = MAX( SI(IJ, SKEN+1), 0.0D0 )
      ELSE IF ( FLDEL(SKEN+1) <  0 ) THEN
         FLTMP = 0.0D0
     &- ( 0.0D0 - 0.0D0 ) / W(IJ, SKEN+1) * DLZI(SKEN+1)
         SI(IJ, SKEN+1) =
     &      MAX( SI(IJ, SKEN+1), MAX( FLTMP, S(IJ, SKEN) ) )
         SI(IJ, SKEN+1) = MIN( SI(IJ, SKEN+1), 0.0D0 )
      ENDIF
      ENDIF
      ENDDO

      ENDIF

      DO K = SKST+1, SKEN+1
      DO IJ = IJTSTR, IJTEND
         FLX(IJ, K) = ( 1.0D0 - PHII(K) ) * W(IJ, K) * SI(IJ, K)
      ENDDO
      ENDDO

      ENDIF

!     -- Diffusion : Crank nicolson

      IF (SDIC) THEN

      DO K = SKST, SKEN
         EAL(K) = ( ( 1.0D0 - PHI(K) ) / TS )
         EBE(K) = ( ( ( 1.0D0 - PHII(K) )   * DBI(K) )
     &            / ( 2.0D0 * DLZ(K) * DLZI(K)   ) )
         EGA(K) = ( ( ( 1.0D0 - PHII(K+1) ) * DBI(K+1) )
     &            / ( 2.0D0 * DLZ(K) * DLZI(K+1) ) )
      ENDDO

      DO K = SKST+1, SKEN-1
         EA(K)     = - EBE(K)
         EB(K)     = EAL(K) + EBE(K) + EGA(K)
         EC(K)     = - EGA(K)
      ENDDO

         EA(SKST)  = - EBE(SKST)
         EB(SKST)  = EAL(SKST) + EGA(SKST)
         EC(SKST)  = - EGA(SKST)
         EA(SKEN)  = - EBE(SKEN)
         EB(SKEN)  = EAL(SKEN) + EBE(SKEN)
         EC(SKEN)  = - EGA(SKEN)

      DO IJ = IJTSTR, IJTEND
         ED(IJ, SKST) =
     &      + ( EAL(SKST) - EGA(SKST) ) * S(IJ, SKST)
     &      + EGA(SKST) * S(IJ, SKST+1)
     &      + ( ( FLX(IJ, SKST) - FLX(IJ, SKST+1) ) / DLZ(SKST) )
     &      + ( ( 1.0D0 - PHI(SKST) ) / 3.1536D+7 )
      ENDDO

      DO K = SKST+1, SKEN-1
      DO IJ = IJTSTR, IJTEND
         ED(IJ, K) =
     &      + EBE(K) * S(IJ, K-1)
     &      + ( EAL(K) - EBE(K) - EGA(K) ) * S(IJ, K)
     &      + EGA(K) * S(IJ, K+1)
     &      + ( ( FLX(IJ, K) - FLX(IJ, K+1) ) / DLZ(K) )
     &      + ( ( 1.0D0 - PHI(K) ) / 3.1536D+7 )
      ENDDO
      ENDDO

      DO IJ = IJTSTR, IJTEND
         ED(IJ, SKEN) =
     &      + EBE(SKEN) * S(IJ, SKEN-1)
     &      + ( EAL(SKEN) - EBE(SKEN) ) * S(IJ, SKEN)
     &      + ( ( FLX(IJ, SKEN) - FLX(IJ, SKEN+1) ) / DLZ(SKEN) )
     &      + ( ( 1.0D0 - PHI(SKEN) ) / 3.1536D+7 )
      ENDDO

      CALL STRIDAG( EA, EB, EC, ED, EX )

      ENDIF

!     -- Diffusion : Implicit

      IF (SDII) THEN

      DO K = SKST, SKEN
         EAL(K) = ( ( 1.0D0 - PHI(K) ) / TS )
         EBE(K) = ( ( ( 1.0D0 - PHII(K) )   * DBI(K) )
     &            / ( DLZ(K) * DLZI(K)   ) )
         EGA(K) = ( ( ( 1.0D0 - PHII(K+1) ) * DBI(K+1) )
     &            / ( DLZ(K) * DLZI(K+1) ) )
      ENDDO

      DO K = SKST+1, SKEN-1
         EA(K)     = - EBE(K)
         EB(K)     = EAL(K) + EBE(K) + EGA(K)
         EC(K)     = - EGA(K)
      ENDDO

         EA(SKST)  = - EBE(SKST)
         EB(SKST)  = EAL(SKST) + EGA(SKST)
         EC(SKST)  = - EGA(SKST)
         EA(SKEN)  = - EBE(SKEN)
         EB(SKEN)  = EAL(SKEN) + EBE(SKEN)
         EC(SKEN)  = - EGA(SKEN)

      DO K = SKST, SKEN
      DO IJ = IJTSTR, IJTEND
         ED(IJ, K) =
     &      + EAL(K) * S(IJ, K)
     &      + ( ( FLX(IJ, K) - FLX(IJ, K+1) ) / DLZ(K) )
     &      + ( ( 1.0D0 - PHI(K) ) / 3.1536D+7 )
      ENDDO
      ENDDO

      CALL STRIDAG( EA, EB, EC, ED, EX )

      ENDIF

!     -- t = t + dt --

      DO K = SKST, SKEN
      DO IJ = IJTSTR, IJTEND
         SN(IJ, K) = MAX( EX(IJ, K), 0.0D0 )
         EX(IJ, K) = 0.0D0
      ENDDO
      ENDDO

      END SUBROUTINE AGE

! =====================================================================
      SUBROUTINE STRIDAG( EA, EB, EC, ED, EX )

      IMPLICIT NONE

#include "sddim.F"

      REAL*8  EA(SDZDIM), EB(SDZDIM), EC(SDZDIM)
      REAL*8  ED(NXY,SDZDIM)
      REAL*8  EX(NXY,SDZDIM)
      REAL*8  EW(NXY,SDZDIM)
      REAL*8  BET(NXY)
      INTEGER IJ, K

      DO IJ = IJTSTR, IJTEND
         BET(IJ)      = EB(SKST)
         EX(IJ, SKST) = ED(IJ, SKST) / BET(IJ)
      ENDDO

      DO K = SKST+1, SKEN
      DO IJ = IJTSTR, IJTEND
         EW(IJ, K) = EC(K-1) / BET(IJ)
         BET(IJ)   = EB(K) - EA(K) * EW(IJ, K)
         EX(IJ, K) = (ED(IJ, K) - EA(K) * EX(IJ, K-1) ) / BET(IJ)
      ENDDO
      ENDDO

      DO K = SKEN-1, SKST, -1
      DO IJ = IJTSTR, IJTEND
         EX(IJ, K) = EX(IJ, K) - EW(IJ, K+1) * EX(IJ, K+1)
      ENDDO
      ENDDO

      END SUBROUTINE STRIDAG

end module sdsolid_module
