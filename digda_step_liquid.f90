!######################################################################

! --- information -----------------------------------------------------
!
!  Sediment model : Dissolved material
!
!  HISTORY
!     '15.01.09  H.Kobayashi :
!
! ---------------------------------------------------------------------

module globals_digda_step_liquid
  use module_zcfile
  use module_digda_grid
  use module_sdtime
  implicit none

      double precision :: DFLX(NXY, SDZDIM+1)

      double precision :: EA(SDZDIM), EB(NXY, SDZDIM), EC(SDZDIM)
      double precision :: EAL(SDZDIM), EBE(SDZDIM), EGA(SDZDIM), EEP(NXY, SDZDIM)
      double precision :: ED(NXY, SDZDIM),  EX(NXY, SDZDIM)

      logical CDCN, CDII
      save CDCN, CDII
      namelist /NMSDLI/ CDCN, CDII
!     -- Diffusion choice --
      data CDCN  / .FALSE. / !! Crank nicolson ?
      data CDII  / .TRUE.  / !! Implicit ?

end module globals_digda_step_liquid

! =====================================================================

module module_digda_step_liquid
  use globals_digda_step_liquid
  implicit none
contains

      subroutine SDLIQUIDSET

      integer :: IJ, K

      WRITE(NFSTDO, *) '*** Carbon : SDLIQUIDSET ***'
      REWIND(NFPARM)
      READ (NFPARM, NMSDLI, END=99)
99 continue
      WRITE(NFSTDO, NMSDLI)

      end subroutine SDLIQUIDSET

! =====================================================================

      subroutine LIQUID( &
     &      BW, RLT, PHI, PHII, DEI, C, &
     &      CN )

      double precision, intent(in)  :: BW(NXY)
      double precision, intent(in)  :: RLT(NXY,SDZDIM)
      double precision, intent(in)  :: PHI(SDZDIM)
      double precision, intent(in)  :: PHII(SDZDIM+1)
      double precision, intent(in)  :: DEI(SDZDIM+1)
      double precision, intent(in)  :: C(NXY,SDZDIM)
      double precision, intent(out) :: CN(NXY,SDZDIM)

      integer :: IJ, K

!     -- Diffusion : Crank nicolson --

      IF (CDCN) THEN

      DO K = SKST, SKEN
         EAL(K) = ( PHI(K) / TS )
         EBE(K) = ( ( PHII(K)   * DEI(K)   ) &
     &            / ( 2.0D0 * DLZ(K) * DLZI(K)   ) )
         EGA(K) = ( ( PHII(K+1) * DEI(K+1) ) &
     &            / ( 2.0D0 * DLZ(K) * DLZI(K+1) ) )
      END DO

      DO K = SKST+1, SKEN-1
         EA(K)         = - EBE(K)
         EC(K)         = - EGA(K)
      DO IJ = IJTSTR, IJTEND
         EB(IJ, K)     = EAL(K) + EBE(K) + EGA(K)
      END DO
      END DO

         EA(SKST)      = - EBE(SKST)
         EC(SKST)      = - EGA(SKST)
         EA(SKEN)      = - EBE(SKEN)
         EC(SKEN)      = - EGA(SKEN)
      DO IJ = IJTSTR, IJTEND
         EB(IJ, SKST)  = EAL(SKST) + EGA(SKST) &
     &   + ( PHII(SKST) * DEI(SKST) / DLZI(SKST) / DLZI(SKST) )
         EB(IJ, SKEN)  = EAL(SKEN) + EBE(SKEN)
      END DO

      DO IJ = IJTSTR, IJTEND
         ED(IJ, SKST) = &
     &      + ( EAL(SKST) - EGA(SKST) ) * C(IJ, SKST) &
     &      + EGA(SKST) * C(IJ, SKST+1) &
     &      + RLT(IJ, SKST) &
     &   + ( PHII(SKST) * DEI(SKST) * BW(IJ) / DLZI(SKST) / DLZI(SKST) )
      END DO

      DO K = SKST+1, SKEN-1
      DO IJ = IJTSTR, IJTEND
         ED(IJ, K) = &
     &      + EBE(K) * C(IJ, K-1) &
     &      + ( EAL(K) - EBE(K) - EGA(K) ) * C(IJ, K) &
     &      + EGA(K) * C(IJ, K+1) &
     &      + RLT(IJ, K)
      END DO
      END DO

      DO IJ = IJTSTR, IJTEND
         ED(IJ, SKEN) = &
     &      + EBE(SKEN) * C(IJ, SKEN-1) &
     &      + ( EAL(SKEN) - EGA(SKEN) ) * C(IJ, SKEN) &
     &      + RLT(IJ, SKEN)
      END DO

      CALL LTRIDAG( EA, EB, EC, ED, EX )

      END IF

!     -- Diffusion : Implicit --

      IF (CDII) THEN

      DO K = SKST, SKEN
         EAL(K) = ( PHI(K) / TS )
         EBE(K) = ( ( PHII(K)   * DEI(K)   ) &
     &            / ( DLZ(K) * DLZI(K)   ) )
         EGA(K) = ( ( PHII(K+1) * DEI(K+1) ) &
     &            / ( DLZ(K) * DLZI(K+1) ) )
      END DO

      DO K = SKST+1, SKEN-1
         EA(K)         = - EBE(K)
         EC(K)         = - EGA(K)
      DO IJ = IJTSTR, IJTEND
         EB(IJ, K)     = EAL(K) + EBE(K) + EGA(K)
      END DO
      END DO

         EA(SKST)      = - EBE(SKST)
         EC(SKST)      = - EGA(SKST)
         EA(SKEN)      = - EBE(SKEN)
         EC(SKEN)      = - EGA(SKEN)
      DO IJ = IJTSTR, IJTEND
         EB(IJ, SKST)  = EAL(SKST) + EGA(SKST) &
     &   + ( PHII(SKST) * DEI(SKST) / DLZI(SKST) / DLZI(SKST) )
         EB(IJ, SKEN)  = EAL(SKEN) + EBE(SKEN)
      END DO

      DO K = SKST, SKEN
      DO IJ = IJTSTR, IJTEND
         ED(IJ, K) = &
     &      + EAL(K) * C(IJ, K) &
     &      + RLT(IJ, K)
      END DO
      END DO

      DO IJ = IJTSTR, IJTEND
         ED(IJ, SKST) = ED(IJ, SKST) &
     &   + ( PHII(SKST) * DEI(SKST) * BW(IJ) / DLZI(SKST) / DLZI(SKST) )
      END DO

      CALL LTRIDAG( EA, EB, EC, ED, EX )

      END IF

!     -- t = t + dt --

      DO K = SKST, SKEN
      DO IJ = IJTSTR, IJTEND
         CN(IJ, K) = MAX( EX(IJ, K), 0.0D0 )
         EX(IJ, K) = 0.0D0
      END DO
      END DO

      end subroutine LIQUID

! =====================================================================

      subroutine LIQUID_O2( &
     &      BW, RLT, PHI, PHII, DEI, C, &
     &      CN )

      double precision, intent(in)  :: BW(NXY)
      double precision, intent(inout) :: RLT(NXY,SDZDIM)
      double precision, intent(in)  :: PHI(SDZDIM)
      double precision, intent(in)  :: PHII(SDZDIM+1)
      double precision, intent(in)  :: DEI(SDZDIM+1)
      double precision, intent(in)  :: C(NXY,SDZDIM)
      double precision, intent(out) :: CN(NXY,SDZDIM)

      integer :: IJ, K

      DO K = SKST, SKEN
      DO IJ = IJTSTR, IJTEND
         IF ( C(IJ, K) /= 0.0D0 ) THEN
            RLT(IJ, K) = RLT(IJ, K) / C(IJ, K)
         ELSE
            RLT(IJ, K) = 0.0D0
         END IF
      END DO
      END DO

!     -- Diffusion : Crank nicolson --

      IF (CDCN) THEN

      DO K = SKST, SKEN
         EAL(K) = ( PHI(K) / TS )
         EBE(K) = ( ( PHII(K)   * DEI(K)   ) &
     &            / ( 2.0D0 * DLZ(K) * DLZI(K)   ) )
         EGA(K) = ( ( PHII(K+1) * DEI(K+1) ) &
     &            / ( 2.0D0 * DLZ(K) * DLZI(K+1) ) )
      DO IJ = IJTSTR, IJTEND
         EEP(IJ, K) = - 0.50D0 * RLT(IJ, K)
      END DO
      END DO

      DO K = SKST+1, SKEN-1
         EA(K)         = - EBE(K)
         EC(K)         = - EGA(K)
      DO IJ = IJTSTR, IJTEND
         EB(IJ, K)     = EAL(K) + EBE(K) + EGA(K) + EEP(IJ, K)
      END DO
      END DO

         EA(SKST)      = - EBE(SKST)
         EC(SKST)      = - EGA(SKST)
         EA(SKEN)      = - EBE(SKEN)
         EC(SKEN)      = - EGA(SKEN)
      DO IJ = IJTSTR, IJTEND
         EB(IJ, SKST)  = EAL(SKST) + EGA(SKST) + EEP(IJ, SKST) &
     &   + ( PHII(SKST) * DEI(SKST) / DLZI(SKST) / DLZI(SKST) )
         EB(IJ, SKEN)  = EAL(SKEN) + EBE(SKEN) + EEP(IJ, SKEN)
      END DO

      DO IJ = IJTSTR, IJTEND
         ED(IJ, SKST) = &
     &      + ( EAL(SKST) - EGA(SKST) ) * C(IJ, SKST) &
     &      + EGA(SKST) * C(IJ, SKST+1) &
     &      + 0.50D0 * RLT(IJ, SKST) * C(IJ, SKST) &
     &   + ( PHII(SKST) * DEI(SKST) * BW(IJ) / DLZI(SKST) / DLZI(SKST) )
      END DO

      DO K = SKST+1, SKEN-1
      DO IJ = IJTSTR, IJTEND
         ED(IJ, K) = &
     &      + EBE(K) * C(IJ, K-1) &
     &      + ( EAL(K) - EBE(K) - EGA(K) ) * C(IJ, K) &
     &      + EGA(K) * C(IJ, K+1) &
     &      + 0.50D0 * RLT(IJ, K) * C(IJ, K)
      END DO
      END DO

      DO IJ = IJTSTR, IJTEND
         ED(IJ, SKEN) = &
     &      + EBE(SKEN) * C(IJ, SKEN-1) &
     &      + ( EAL(SKEN) - EGA(SKEN) ) * C(IJ, SKEN) &
     &      + 0.50D0 * RLT(IJ, SKEN) * C(IJ, SKEN)
      END DO

      CALL LTRIDAG( EA, EB, EC, ED, EX )

      END IF

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
      END DO
      END DO

      DO K = SKST+1, SKEN-1
         EA(K)         = - EBE(K)
         EC(K)         = - EGA(K)
      DO IJ = IJTSTR, IJTEND
         EB(IJ, K)     = EAL(K) + EBE(K) + EGA(K) + EEP(IJ, K)
      END DO
      END DO

         EA(SKST)      = - EBE(SKST)
         EC(SKST)      = - EGA(SKST)
         EA(SKEN)      = - EBE(SKEN)
         EC(SKEN)      = - EGA(SKEN)
      DO IJ = IJTSTR, IJTEND
         EB(IJ, SKST)  = EAL(SKST) + EGA(SKST) + EEP(IJ, SKST) &
     &   + ( PHII(SKST) * DEI(SKST) / DLZI(SKST) / DLZI(SKST) )
         EB(IJ, SKEN)  = EAL(SKEN) + EBE(SKEN) + EEP(IJ, SKEN)
      END DO

      DO K = SKST, SKEN
      DO IJ = IJTSTR, IJTEND
         ED(IJ, K) = &
     &      + EAL(K) * C(IJ, K)
      END DO
      END DO

      DO IJ = IJTSTR, IJTEND
         ED(IJ, SKST) = ED(IJ, SKST) &
     &   + ( PHII(SKST) * DEI(SKST) * BW(IJ) / DLZI(SKST) / DLZI(SKST) )
      END DO

      CALL LTRIDAG( EA, EB, EC, ED, EX )

      END IF

!     -- t = t + dt --

      DO K = SKST, SKEN
      DO IJ = IJTSTR, IJTEND
         CN(IJ, K) = MAX( EX(IJ, K), 0.0D0 )
         EX(IJ, K) = 0.0D0
      END DO
      END DO

      end subroutine LIQUID_O2

! =====================================================================

      subroutine LTRIDAG( EA, EB, EC, ED, EX )

      double precision, intent(in)  :: EA(SDZDIM), EB(NXY,SDZDIM), EC(SDZDIM)
      double precision, intent(in)  :: ED(NXY,SDZDIM)
      double precision, intent(out) :: EX(NXY,SDZDIM)

      double precision :: EW(NXY,SDZDIM)
      double precision :: BET(NXY)

      integer :: IJ, K

      DO IJ = IJTSTR, IJTEND
         BET(IJ)      = EB(IJ, SKST)
         EX(IJ, SKST) = ED(IJ, SKST) / BET(IJ)
      END DO

      DO K = SKST+1, SKEN
      DO IJ = IJTSTR, IJTEND
         EW(IJ, K) = EC(K-1) / BET(IJ)
         BET(IJ)   = EB(IJ, K) - EA(K) * EW(IJ, K)
         EX(IJ, K) = (ED(IJ, K) - EA(K) * EX(IJ, K-1) ) / BET(IJ)
      END DO
      END DO

      DO K = SKEN-1, SKST, -1
      DO IJ = IJTSTR, IJTEND
         EX(IJ, K) = EX(IJ, K) - EW(IJ, K+1) * EX(IJ, K+1)
      END DO
      END DO

      end subroutine LTRIDAG

end module module_digda_step_liquid
