!######################################################################

! --- information -----------------------------------------------------
!
!  Calculate [CO2(g)] from chemical equilibrium in ocean
!
!  HISTORY
!     '04.04.23  A.Oka :
!     '06.06.03  A.Oka : ocean density [mMol/m^3] --> [Mol/kg]
!
! ---------------------------------------------------------------------

module globals_sdchemeq
      use zcfile_module
      use zcphys_module
      use sddim_module
      use zchem_module
      implicit none

      REAL*8  CTIN(NXY, SDZDIM) !! Total carbon [mMol/m^3] (input)
      REAL*8  ATIN(NXY, SDZDIM) !! Alkalinity [mMol/m^3] (input)
      REAL*8  SITIN(NXY)        !! Total silicon [mMol/m^3] (input)
      REAL*8  PTIN(NXY)         !! Total phosphorous [mMol/m^3] (input)

      REAL*8  TEM(NXY)          !! Temperature
      REAL*8  SAL(NXY)          !! Salinity
      REAL*8  PCO2O(NXY)        !! [CO2(g)] [ppmv]
      REAL*8  ALP(NXY)          !! Solubility of CO2 [mMol/m^3/ppmv]
      REAL*8  CT(NXY, SDZDIM)   !! Total carbon [Mol/kg]
      REAL*8  AT(NXY, SDZDIM)   !! Alkalinity [Mol/kg]
      REAL*8  SIT(NXY)          !! Total silicon [Mol/kg]
      REAL*8  PT(NXY)           !! Total phosphorous [Mol/kg]

      REAL*8  CO2(NXY, SDZDIM)  !! surface [CO2]   [mMol/m^3]
      REAL*8  HCO3(NXY, SDZDIM) !! surface [HCO3-] [mMol/m^3]
      REAL*8  CO32(NXY, SDZDIM) !! surface [CO32-] [mMol/m^3]

      REAL*8  K0(NXY),  K1(NXY),   K2(NXY)
      REAL*8  KB(NXY),  KW(NXY)
      REAL*8  K1P(NXY), K2P(NXY),  K3P(NXY)
      REAL*8  KSI(NXY), KSO4(NXY), KF(NXY)

      REAL*8  BT(NXY)           !! Total boric acid [Mol/kg]
      REAL*8  SO4T(NXY)         !! Total SO4 [Mol/kg]
      REAL*8  FLT(NXY)          !! Total fluorosis [Mol/kg]

      REAL*8  PH(NXY, SDZDIM), CU, RCU, PM
      SAVE    PH, CU, RCU, PM

      REAL*8  H, X, XLW, XHI, F, DELT
      INTEGER IJ, K, ITR
      LOGICAL OFIRST
      SAVE OFIRST

      REAL*8  PHINI, DPHINI, DPH, ALKERR, ROCN
      INTEGER ITRMAX
      SAVE ITRMAX, DPH, DPHINI, ALKERR

      NAMELIST /NMSDCHEM/ PHINI, DPHINI, DPH, ALKERR, ITRMAX, ROCN

      DATA PHINI   / 8.0 /
      DATA DPHINI  / 8.0 /
      DATA DPH     / 3.0 /
      DATA ALKERR  / 0.1 /
      DATA ITRMAX  / 10000 /

      DATA ROCN / 1.0245 /  !! ocean density [kg/L]

end module globals_sdchemeq

! =====================================================================

module sdchemeq_module
    use globals_sdchemeq
    implicit none
contains

      SUBROUTINE SDCHEMEQSET

      WRITE(NFSTDO, *) '*** Carbon : SDCHEMEQSET ***'
      REWIND(NFPARM)
      READ (NFPARM, NMSDCHEM, END=99)
 99   WRITE(NFSTDO, NMSDCHEM)

      DO IJ = 1, NXY
         DO K = 1, SDZDIM
            PH(IJ, K) = PHINI
         ENDDO
      ENDDO

      PM  = 1.0D+6          !! [atm] --> [ppmv]

      CU  = 1.0D-6 / ROCN   !! [mMol/m^3] --> [Mol/kg]
      RCU = 1.0D+6 * ROCN   !! [Mol/kg] --> [mMol/m^3]
      ALKERR = ALKERR * CU

      OFIRST = .TRUE.

      END SUBROUTINE SDCHEMEQSET

! =====================================================================

      SUBROUTINE SDCHEMEQ( &
     &     TEM, SAL, CTIN, ATIN, &
     &     SITIN, PTIN, &
     &     CO32 )

      double precision, intent(in)  :: TEM(NXY)
      double precision, intent(in)  :: SAL(NXY)
      double precision, intent(in)  :: CTIN(NXY, SDZDIM)
      double precision, intent(in)  :: ATIN(NXY, SDZDIM)
      double precision, intent(in)  :: SITIN(NXY)
      double precision, intent(in)  :: PTIN(NXY)
      double precision, intent(out) :: CO32(NXY, SDZDIM)

      DO IJ = IJTSTR, IJTEND
!        -- [mMol/m^3]-->[Mol/kg] --
         DO K = SKST, SKEN
         CT(IJ, K)  = CTIN(IJ, K)  * CU
         AT(IJ, K)  = ATIN(IJ, K)  * CU
         ENDDO
         SIT(IJ) = SITIN(IJ) * CU
         PT(IJ)  = PTIN(IJ)  * CU
      ENDDO

      DO IJ = IJTSTR, IJTEND
      IF ( TEM(IJ) /= 0.0D0 .and. SAL(IJ) /= 0.0D0 ) THEN
!        -- Chemical equilibrium constant --
         K0(IJ)   = FK0( TEM(IJ)+kelvin, SAL(IJ) )
         K1(IJ)   = FK1( TEM(IJ)+kelvin, SAL(IJ) )
         K2(IJ)   = FK2( TEM(IJ)+kelvin, SAL(IJ) )
         KB(IJ)   = FKB( TEM(IJ)+kelvin, SAL(IJ) )
         KW(IJ)   = FKW( TEM(IJ)+kelvin, SAL(IJ) )
         K1P(IJ)  = FK1P( TEM(IJ)+kelvin, SAL(IJ) )
         K2P(IJ)  = FK2P( TEM(IJ)+kelvin, SAL(IJ) )
         K3P(IJ)  = FK3P( TEM(IJ)+kelvin, SAL(IJ) )
         KSI(IJ)  = FKSI( TEM(IJ)+kelvin, SAL(IJ) )
         KSO4(IJ) = FKSO4( TEM(IJ)+kelvin, SAL(IJ) )
         KF(IJ)   = FKF( TEM(IJ)+kelvin, SAL(IJ), KSO4(IJ) )
!        -- Total [Mol/kg] --
         BT(IJ)   = FBT( SAL(IJ) )
         SO4T(IJ) = FSO4T( SAL(IJ) )
         FLT(IJ)  = FFLT( SAL(IJ) )
      ENDIF
      ENDDO

      IF (OFIRST) THEN
         DELT = DPHINI
         OFIRST = .FALSE.
      ELSE
         DELT = DPH
      ENDIF

      DO K = SKST, SKEN
      DO IJ = IJTSTR, IJTEND
      IF ( TEM(IJ) /= 0.0D0 .and. SAL(IJ) /= 0.0D0 ) THEN

         XLW = PH(IJ, K) - DELT
         XHI = PH(IJ, K) + DELT

!        -- Iteration : bisection method --
         DO ITR = 1, ITRMAX

            X = 0.5D0 * ( XLW + XHI )
            H = 10.0D0 ** ( - X )

!           IF ( AMSKT(IJ, KSTR). EQ. 0.0D0 ) GOTO 999
            IF ( TEM(IJ) .EQ. 0.0D0 .and. SAL(IJ) .EQ. 0.0D0 ) GOTO 999

            F  = FALK( &
     &           H, &
     &           K1(IJ),  K2(IJ), &
     &           KB(IJ),  KW(IJ), &
     &           K1P(IJ), K2P(IJ),  K3P(IJ), &
     &           KSI(IJ), KSO4(IJ), KF(IJ), &
     &           CT(IJ, K),  BT(IJ),   PT(IJ), &
     &           SIT(IJ), SO4T(IJ), FLT(IJ) )

            IF ( ABS(F-AT(IJ, K)).LT.ALKERR ) GOTO 999

            IF ( ( F - AT(IJ, K) ).GT.0.D0 ) THEN
               XHI = X
            ELSE
               XLW = X
            ENDIF

         ENDDO

         WRITE(NFSTDO,*) ' ITERATION COUNT OVER (SEDIMENT) :', IJ
         WRITE(NFSTDO,*) '    PH  = ', X
         WRITE(NFSTDO,*) '    F   = ', F * RCU
         WRITE(NFSTDO,*) '    ALK = ', ATIN(IJ, K)
         WRITE(NFSTDO,*) '    F - ALK = ', F * RCU - ATIN(IJ, K)
         WRITE(NFSTDO,*) '            < ALKERR = ', ALKERR * RCU
         WRITE(NFSTDO,*) '   [Info] '
!        WRITE(NFSTDO,*) '    at  I = ', MOD(IJ,NXDIM) - ISTR + 1
!        WRITE(NFSTDO,*) '        J = ', IJ/NXDIM - JSTR + 2
!        WRITE(NFSTDO,*) '        K = ', K - SKST + 1
         WRITE(NFSTDO,*) '    TEM = ', TEM(IJ)
         WRITE(NFSTDO,*) '    SAL = ', SAL(IJ)
         WRITE(NFSTDO,*) '    DIC = ', CTIN(IJ, K)

         XLW = PH(IJ, K) - DELT
         XHI = PH(IJ, K) + DELT
         DO ITR = 1, ITRMAX
            X = 0.5D0 * ( XLW + XHI )
            H = 10.0D0 ** ( - X )
            F  = FALK( &
     &           H, &
     &           K1(IJ),  K2(IJ), &
     &           KB(IJ),  KW(IJ), &
     &           K1P(IJ), K2P(IJ),  K3P(IJ), &
     &           KSI(IJ), KSO4(IJ), KF(IJ), &
     &           CT(IJ, K),  BT(IJ),   PT(IJ), &
     &           SIT(IJ), SO4T(IJ), FLT(IJ) )
            IF ( ABS(F-AT(IJ, K)).LT.ALKERR ) GOTO 999
            IF ( ( F - AT(IJ, K) ).GT.0.D0 ) THEN
               XHI = X
            ELSE
               XLW = X
            ENDIF
         ENDDO

!        STOP

 999     CONTINUE
         PH(IJ, K) = X

      ENDIF
      ENDDO
      ENDDO

      DO K = SKST, SKEN
      DO IJ = IJTSTR, IJTEND
      IF ( TEM(IJ) /= 0.0D0 .and. SAL(IJ) /= 0.0D0 ) THEN
         X = 10.D0 ** ( - PH(IJ, K) )
!        PCO2O(IJ) = CT(IJ) *X*X &
!    &               / ( X*X + X*K1(IJ) + K1(IJ)*K2(IJ) ) &
!    &               / K0(IJ) &
!    &               * PM      !! [atm]-->[ppmv]
!        ALP(IJ) = K0(IJ)      !! [Mol/kg/atm] = [mMol/m^3/ppmv]
!        -- CO2 --
         CO2(IJ, K)  = CT(IJ, K) *X*X &
     &               / ( X*X + X*K1(IJ) + K1(IJ)*K2(IJ) ) / CU
!        -- BI-CARBONATE ION --
         HCO3(IJ, K) = K1(IJ) * CO2(IJ, K) / X
!        -- CARBONATE ION --
         CO32(IJ, K) = K2(IJ) * HCO3(IJ, K) / X
      ENDIF
      ENDDO
      ENDDO

      END SUBROUTINE SDCHEMEQ

end module sdchemeq_module
