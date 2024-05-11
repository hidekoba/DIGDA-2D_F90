!######################################################################

! --- information -----------------------------------------------------
!
!  Sediment model : Dissolved material
!
!  HISTORY
!     '15.01.09  H.Kobayashi :
!
! ---------------------------------------------------------------------

module globals_sdliquid
      use zcfile_module
      use sddim_module
      implicit none

      REAL*8 C(NXY,SDZDIM)
      REAL*8 CN(NXY,SDZDIM)
      REAL*8 BW(NXY)
      REAL*8 DFLX(NXY,SDZDIM+1)
      REAL*8 RLT(NXY,SDZDIM)

      REAL*8 PHI(SDZDIM)
      REAL*8 PHII(SDZDIM+1)
      REAL*8 DEI(SDZDIM+1)

      REAL*8 EA(SDZDIM), EB(NXY, SDZDIM), EC(SDZDIM)
      REAL*8 EAL(SDZDIM), EBE(SDZDIM), EGA(SDZDIM), EEP(NXY,SDZDIM)
      REAL*8 ED(NXY,SDZDIM),  EX(NXY,SDZDIM)

      INTEGER IJ, K

      LOGICAL CDIC, CDII
      SAVE CDIC, CDII
      NAMELIST /NMSDLI/ CDIC, CDII
!     -- Diffusion choice --
      DATA CDIC  / .FALSE. / !! Crank nicolson ?
      DATA CDII  / .TRUE.  / !! Implicit ?

end module globals_sdliquid

! =====================================================================

module sdliquid_module
      use globals_sdliquid
      implicit none
contains

      SUBROUTINE SDLIQUIDSET

      WRITE(NFSTDO, *) '*** Carbon : SDLIQUIDSET ***'
      REWIND(NFPARM)
      READ (NFPARM, NMSDLI, END=99)
 99   WRITE(NFSTDO, NMSDLI)

!     -- Convert unit --
      DO K = SKST, SKEN
         LZ(K)   = LZ(K)   * 1.0D+2            !! [m] --> [cm]
         DLZ(K)  = DLZ(K)  * 1.0D+2            !! [m] --> [cm]
      ENDDO
      DO K = SKST, SKEN+1
         LZI(K)  = LZI(K)  * 1.0D+2            !! [m] --> [cm]
         DLZI(K) = DLZI(K) * 1.0D+2            !! [m] --> [cm]
      ENDDO

      END SUBROUTINE SDLIQUIDSET

! =====================================================================

      SUBROUTINE LIQUID( &
     &      BW, RLT, PHI, PHII, DEI, C, &
     &      CN )

      double precision, intent(in)  :: BW(NXY)
      double precision, intent(in)  :: RLT(NXY,SDZDIM)
      double precision, intent(in)  :: PHI(SDZDIM)
      double precision, intent(in)  :: PHII(SDZDIM+1)
      double precision, intent(in)  :: DEI(SDZDIM+1)
      double precision, intent(in)  :: C(NXY,SDZDIM)
      double precision, intent(out) :: CN(NXY,SDZDIM)

!     -- Diffusion : Crank nicolson --

      IF (CDIC) THEN

      DO K = SKST, SKEN
         EAL(K) = ( PHI(K) / TS )
         EBE(K) = ( ( PHII(K)   * DEI(K)   ) &
     &            / ( 2.0D0 * DLZ(K) * DLZI(K)   ) )
         EGA(K) = ( ( PHII(K+1) * DEI(K+1) ) &
     &            / ( 2.0D0 * DLZ(K) * DLZI(K+1) ) )
      ENDDO

      DO K = SKST+1, SKEN-1
         EA(K)         = - EBE(K)
         EC(K)         = - EGA(K)
      DO IJ = IJTSTR, IJTEND
         EB(IJ, K)     = EAL(K) + EBE(K) + EGA(K)
      ENDDO
      ENDDO

         EA(SKST)      = - EBE(SKST)
         EC(SKST)      = - EGA(SKST)
         EA(SKEN)      = - EBE(SKEN)
         EC(SKEN)      = - EGA(SKEN)
      DO IJ = IJTSTR, IJTEND
         EB(IJ, SKST)  = EAL(SKST) + EGA(SKST) &
     &   + ( PHII(SKST) * DEI(SKST) / DLZI(SKST) / DLZI(SKST) )
         EB(IJ, SKEN)  = EAL(SKEN) + EBE(SKEN)
      ENDDO

      DO IJ = IJTSTR, IJTEND
         ED(IJ, SKST) = &
     &      + ( EAL(SKST) - EGA(SKST) ) * C(IJ, SKST) &
     &      + EGA(SKST) * C(IJ, SKST+1) &
     &      + RLT(IJ, SKST) &
     &   + ( PHII(SKST) * DEI(SKST) * BW(IJ) / DLZI(SKST) / DLZI(SKST) )
      ENDDO

      DO K = SKST+1, SKEN-1
      DO IJ = IJTSTR, IJTEND
         ED(IJ, K) = &
     &      + EBE(K) * C(IJ, K-1) &
     &      + ( EAL(K) - EBE(K) - EGA(K) ) * C(IJ, K) &
     &      + EGA(K) * C(IJ, K+1) &
     &      + RLT(IJ, K)
      ENDDO
      ENDDO

      DO IJ = IJTSTR, IJTEND
         ED(IJ, SKEN) = &
     &      + EBE(SKEN) * C(IJ, SKEN-1) &
     &      + ( EAL(SKEN) - EGA(SKEN) ) * C(IJ, SKEN) &
     &      + RLT(IJ, SKEN)
      ENDDO

      CALL LTRIDAG( EA, EB, EC, ED, EX )

      ENDIF

!     -- Diffusion : Implicit --

      IF (CDII) THEN

      DO K = SKST, SKEN
         EAL(K) = ( PHI(K) / TS )
         EBE(K) = ( ( PHII(K)   * DEI(K)   ) &
     &            / ( DLZ(K) * DLZI(K)   ) )
         EGA(K) = ( ( PHII(K+1) * DEI(K+1) ) &
     &            / ( DLZ(K) * DLZI(K+1) ) )
      ENDDO

      DO K = SKST+1, SKEN-1
         EA(K)         = - EBE(K)
         EC(K)         = - EGA(K)
      DO IJ = IJTSTR, IJTEND
         EB(IJ, K)     = EAL(K) + EBE(K) + EGA(K)
      ENDDO
      ENDDO

         EA(SKST)      = - EBE(SKST)
         EC(SKST)      = - EGA(SKST)
         EA(SKEN)      = - EBE(SKEN)
         EC(SKEN)      = - EGA(SKEN)
      DO IJ = IJTSTR, IJTEND
         EB(IJ, SKST)  = EAL(SKST) + EGA(SKST) &
     &   + ( PHII(SKST) * DEI(SKST) / DLZI(SKST) / DLZI(SKST) )
         EB(IJ, SKEN)  = EAL(SKEN) + EBE(SKEN)
      ENDDO

      DO K = SKST, SKEN
      DO IJ = IJTSTR, IJTEND
         ED(IJ, K) = &
     &      + EAL(K) * C(IJ, K) &
     &      + RLT(IJ, K)
      ENDDO
      ENDDO

      DO IJ = IJTSTR, IJTEND
         ED(IJ, SKST) = ED(IJ, SKST) &
     &   + ( PHII(SKST) * DEI(SKST) * BW(IJ) / DLZI(SKST) / DLZI(SKST) )
      ENDDO

      CALL LTRIDAG( EA, EB, EC, ED, EX )

      ENDIF

!     -- t = t + dt --

      DO K = SKST, SKEN
      DO IJ = IJTSTR, IJTEND
         CN(IJ, K) = MAX( EX(IJ, K), 0.0D0 )
         EX(IJ, K) = 0.0D0
      ENDDO
      ENDDO

      END SUBROUTINE LIQUID

! =====================================================================

      SUBROUTINE LIQUID_O2( &
     &      BW, RLT, PHI, PHII, DEI, C, &
     &      CN )

      double precision, intent(in)  :: BW(NXY)
      double precision, intent(inout) :: RLT(NXY,SDZDIM)
      double precision, intent(in)  :: PHI(SDZDIM)
      double precision, intent(in)  :: PHII(SDZDIM+1)
      double precision, intent(in)  :: DEI(SDZDIM+1)
      double precision, intent(in)  :: C(NXY,SDZDIM)
      double precision, intent(out) :: CN(NXY,SDZDIM)


      DO K = SKST, SKEN
      DO IJ = IJTSTR, IJTEND
         IF ( C(IJ, K) /= 0.0D0 ) THEN
            RLT(IJ, K) = RLT(IJ, K) / C(IJ, K)
         ELSE
            RLT(IJ, K) = 0.0D0
         ENDIF
      ENDDO
      ENDDO

!     -- Diffusion : Crank nicolson --

      IF (CDIC) THEN

      DO K = SKST, SKEN
         EAL(K) = ( PHI(K) / TS )
         EBE(K) = ( ( PHII(K)   * DEI(K)   ) &
     &            / ( 2.0D0 * DLZ(K) * DLZI(K)   ) )
         EGA(K) = ( ( PHII(K+1) * DEI(K+1) ) &
     &            / ( 2.0D0 * DLZ(K) * DLZI(K+1) ) )
      DO IJ = IJTSTR, IJTEND
         EEP(IJ, K) = - 0.50D0 * RLT(IJ, K)
      ENDDO
      ENDDO

      DO K = SKST+1, SKEN-1
         EA(K)         = - EBE(K)
         EC(K)         = - EGA(K)
      DO IJ = IJTSTR, IJTEND
         EB(IJ, K)     = EAL(K) + EBE(K) + EGA(K) + EEP(IJ, K)
      ENDDO
      ENDDO

         EA(SKST)      = - EBE(SKST)
         EC(SKST)      = - EGA(SKST)
         EA(SKEN)      = - EBE(SKEN)
         EC(SKEN)      = - EGA(SKEN)
      DO IJ = IJTSTR, IJTEND
         EB(IJ, SKST)  = EAL(SKST) + EGA(SKST) + EEP(IJ, SKST) &
     &   + ( PHII(SKST) * DEI(SKST) / DLZI(SKST) / DLZI(SKST) )
         EB(IJ, SKEN)  = EAL(SKEN) + EBE(SKEN) + EEP(IJ, SKEN)
      ENDDO

      DO IJ = IJTSTR, IJTEND
         ED(IJ, SKST) = &
     &      + ( EAL(SKST) - EGA(SKST) ) * C(IJ, SKST) &
     &      + EGA(SKST) * C(IJ, SKST+1) &
     &      + 0.50D0 * RLT(IJ, SKST) * C(IJ, SKST) &
     &   + ( PHII(SKST) * DEI(SKST) * BW(IJ) / DLZI(SKST) / DLZI(SKST) )
      ENDDO

      DO K = SKST+1, SKEN-1
      DO IJ = IJTSTR, IJTEND
         ED(IJ, K) = &
     &      + EBE(K) * C(IJ, K-1) &
     &      + ( EAL(K) - EBE(K) - EGA(K) ) * C(IJ, K) &
     &      + EGA(K) * C(IJ, K+1) &
     &      + 0.50D0 * RLT(IJ, K) * C(IJ, K)
      ENDDO
      ENDDO

      DO IJ = IJTSTR, IJTEND
         ED(IJ, SKEN) = &
     &      + EBE(SKEN) * C(IJ, SKEN-1) &
     &      + ( EAL(SKEN) - EGA(SKEN) ) * C(IJ, SKEN) &
     &      + 0.50D0 * RLT(IJ, SKEN) * C(IJ, SKEN)
      ENDDO

      CALL LTRIDAG( EA, EB, EC, ED, EX )

      ENDIF

!     -- Diffusion : Implicit --

      IF (CDII) THEN

      DO K = SKST, SKEN
         EAL(K) = ( PHI(K) / TS )
         EBE(K) = ( ( PHII(K)   * DEI(K)   ) &
     &            / ( DLZ(K) * DLZI(K)   ) )
         EGA(K) = ( ( PHII(K+1) * DEI(K+1) ) &
     &            / ( DLZ(K) * DLZI(K+1) ) )
      DO IJ = IJTSTR, IJTEND
         EEP(IJ, K) = - RLT(IJ, K)
      ENDDO
      ENDDO

      DO K = SKST+1, SKEN-1
         EA(K)         = - EBE(K)
         EC(K)         = - EGA(K)
      DO IJ = IJTSTR, IJTEND
         EB(IJ, K)     = EAL(K) + EBE(K) + EGA(K) + EEP(IJ, K)
      ENDDO
      ENDDO

         EA(SKST)      = - EBE(SKST)
         EC(SKST)      = - EGA(SKST)
         EA(SKEN)      = - EBE(SKEN)
         EC(SKEN)      = - EGA(SKEN)
      DO IJ = IJTSTR, IJTEND
         EB(IJ, SKST)  = EAL(SKST) + EGA(SKST) + EEP(IJ, SKST) &
     &   + ( PHII(SKST) * DEI(SKST) / DLZI(SKST) / DLZI(SKST) )
         EB(IJ, SKEN)  = EAL(SKEN) + EBE(SKEN) + EEP(IJ, SKEN)
      ENDDO

      DO K = SKST, SKEN
      DO IJ = IJTSTR, IJTEND
         ED(IJ, K) = &
     &      + EAL(K) * C(IJ, K)
      ENDDO
      ENDDO

      DO IJ = IJTSTR, IJTEND
         ED(IJ, SKST) = ED(IJ, SKST) &
     &   + ( PHII(SKST) * DEI(SKST) * BW(IJ) / DLZI(SKST) / DLZI(SKST) )
      ENDDO

      CALL LTRIDAG( EA, EB, EC, ED, EX )

      ENDIF

!     -- t = t + dt --

      DO K = SKST, SKEN
      DO IJ = IJTSTR, IJTEND
         CN(IJ, K) = MAX( EX(IJ, K), 0.0D0 )
         EX(IJ, K) = 0.0D0
      ENDDO
      ENDDO

      END SUBROUTINE LIQUID_O2

! =====================================================================

      SUBROUTINE LTRIDAG( EA, EB, EC, ED, EX )

      double precision, intent(in)  :: EA(SDZDIM), EB(NXY,SDZDIM), EC(SDZDIM)
      double precision, intent(in)  :: ED(NXY,SDZDIM)
      double precision, intent(out) :: EX(NXY,SDZDIM)

      double precision :: EW(NXY,SDZDIM)
      double precision :: BET(NXY)

      integer :: IJ, K

      DO IJ = IJTSTR, IJTEND
         BET(IJ)      = EB(IJ, SKST)
         EX(IJ, SKST) = ED(IJ, SKST) / BET(IJ)
      ENDDO

      DO K = SKST+1, SKEN
      DO IJ = IJTSTR, IJTEND
         EW(IJ, K) = EC(K-1) / BET(IJ)
         BET(IJ)   = EB(IJ, K) - EA(K) * EW(IJ, K)
         EX(IJ, K) = (ED(IJ, K) - EA(K) * EX(IJ, K-1) ) / BET(IJ)
      ENDDO
      ENDDO

      DO K = SKEN-1, SKST, -1
      DO IJ = IJTSTR, IJTEND
         EX(IJ, K) = EX(IJ, K) - EW(IJ, K+1) * EX(IJ, K+1)
      ENDDO
      ENDDO

      END SUBROUTINE LTRIDAG

end module sdliquid_module
