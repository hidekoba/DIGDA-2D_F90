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

module globals_calc_digda_chemeq
  use module_zcfile
  use module_zcphys
  use module_digda_grid
  use module_zchem
  implicit none

  double precision :: PHINI = 8.0
  double precision :: DPHINI = 8.0
  double precision :: DPH = 3.0
  double precision :: ALKERR = 0.1
  double precision :: ROCN = 1.0245 !! ocean density [kg/L]
  integer :: ITRMAX = 10000

  namelist /NMSDCHEM/ PHINI, DPHINI, DPH, ALKERR, ITRMAX, ROCN

end module globals_calc_digda_chemeq

! =====================================================================

module module_digda_chemeq
  use globals_calc_digda_chemeq
  implicit none

  double precision, allocatable :: PH(:,:)
  double precision :: CU, RCU, PM
  logical :: OFIRST

contains

  subroutine set_digda_chemeq

    integer :: IJ, K

    WRITE(NFSTDO, *) '*** Sediment : set_digda_chemeq ***'
    REWIND(NFPARM)
    READ (NFPARM, NMSDCHEM, END=99)
99 CONTINUE
    WRITE(NFSTDO, NMSDCHEM)

    DO IJ = 1, NXY
      DO K = 1, SDZDIM
        PH(IJ, K) = PHINI
      END DO
    END DO

    PM  = 1.0D+6          !! [atm] --> [ppmv]

    CU  = 1.0D-6 / ROCN   !! [mMol/m^3] --> [Mol/kg]
    RCU = 1.0D+6 * ROCN   !! [Mol/kg] --> [mMol/m^3]
    ALKERR = ALKERR * CU

    OFIRST = .true.

  end subroutine set_digda_chemeq

  subroutine calc_digda_chemeq( &
&   TEM, SAL, CTIN, ATIN, &
&   SITIN, PTIN, &
&   CO32 )

    double precision, intent(in)  :: TEM(NXY)          !! Temperature
    double precision, intent(in)  :: SAL(NXY)          !! Salinity
    double precision, intent(in)  :: CTIN(NXY, SDZDIM) !! Total carbon [mMol/m^3]
    double precision, intent(in)  :: ATIN(NXY, SDZDIM) !! Alkalinity [mMol/m^3]
    double precision, intent(in)  :: SITIN(NXY)        !! Total silicon [mMol/m^3]
    double precision, intent(in)  :: PTIN(NXY)         !! Total phosphorous [mMol/m^3]
    double precision, intent(out) :: CO32(NXY, SDZDIM) !! [CO32-] [mMol/m^3]

    double precision :: PCO2O(NXY)        !! [CO2(g)] [ppmv]
    double precision :: ALP(NXY)          !! Solubility of CO2 [mMol/m^3/ppmv]
    double precision :: CT(NXY, SDZDIM)   !! Total carbon [Mol/kg]
    double precision :: AT(NXY, SDZDIM)   !! Alkalinity [Mol/kg]
    double precision :: SIT(NXY)          !! Total silicon [Mol/kg]
    double precision :: PT(NXY)           !! Total phosphorous [Mol/kg]
    double precision :: BT(NXY)           !! Total boric acid [Mol/kg]
    double precision :: SO4T(NXY)         !! Total SO4 [Mol/kg]
    double precision :: FLT(NXY)          !! Total fluorosis [Mol/kg]
    double precision :: CO2(NXY, SDZDIM)  !! [CO2]   [mMol/m^3]
    double precision :: HCO3(NXY, SDZDIM) !! [HCO3-] [mMol/m^3]

    double precision :: K0(NXY),  K1(NXY),   K2(NXY)
    double precision :: KB(NXY),  KW(NXY)
    double precision :: K1P(NXY), K2P(NXY),  K3P(NXY)
    double precision :: KSI(NXY), KSO4(NXY), KF(NXY)

    double precision :: H, X, XLW, XHI, F, DELT

    integer :: IJ, K, ITR

    DO IJ = IJTSTR, IJTEND
    ! -- [mMol/m^3]-->[Mol/kg] --
      DO K = SKST, SKEN
        CT(IJ, K)  = CTIN(IJ, K)  * CU
        AT(IJ, K)  = ATIN(IJ, K)  * CU
      END DO
      SIT(IJ) = SITIN(IJ) * CU
      PT(IJ)  = PTIN(IJ)  * CU
    END DO

    DO IJ = IJTSTR, IJTEND
   !IF ( TEM(IJ) /= 0.0D0 .and. SAL(IJ) /= 0.0D0 ) THEN
    IF ( TEM(IJ) /= rmiss .and. SAL(IJ) /= rmiss ) THEN
    ! -- Chemical equilibrium constant --
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
!     -- Total [Mol/kg] --
      BT(IJ)   = FBT( SAL(IJ) )
      SO4T(IJ) = FSO4T( SAL(IJ) )
      FLT(IJ)  = FFLT( SAL(IJ) )
    END IF
    END DO

    IF (OFIRST) THEN
      DELT = DPHINI
      OFIRST = .false.
    ELSE
      DELT = DPH
    END IF

    DO K = SKST, SKEN
    DO IJ = IJTSTR, IJTEND
   !IF ( TEM(IJ) /= 0.0D0 .and. SAL(IJ) /= 0.0D0 ) THEN
    IF ( TEM(IJ) /= rmiss .and. SAL(IJ) /= rmiss ) THEN

      XLW = PH(IJ, K) - DELT
      XHI = PH(IJ, K) + DELT

      !-- Iteration : bisection method --
      DO ITR = 1, ITRMAX
        X = 0.5D0 * ( XLW + XHI )
        H = 10.0D0 ** ( - X )
       !IF ( AMSKT(IJ, KSTR). EQ. 0.0D0 ) GOTO 999
       !IF ( TEM(IJ) .EQ. 0.0D0 .and. SAL(IJ) .EQ. 0.0D0 ) GOTO 999
        IF ( TEM(IJ) .EQ. rmiss .and. SAL(IJ) .EQ. rmiss ) GOTO 999
        F = FALK( &
     &      H, &
     &      K1(IJ),  K2(IJ), &
     &      KB(IJ),  KW(IJ), &
     &      K1P(IJ), K2P(IJ),  K3P(IJ), &
     &      KSI(IJ), KSO4(IJ), KF(IJ), &
     &      CT(IJ, K),  BT(IJ),   PT(IJ), &
     &      SIT(IJ), SO4T(IJ), FLT(IJ) )
        IF ( ABS(F-AT(IJ, K)).LT.ALKERR ) GOTO 999
        IF ( ( F - AT(IJ, K) ).GT.0.D0 ) THEN
          XHI = X
        ELSE
          XLW = X
        END IF
      END DO

      WRITE(NFSTDO,*) ' ITERATION COUNT OVER (SEDIMENT) :', IJ
      WRITE(NFSTDO,*) '    PH  = ', X
      WRITE(NFSTDO,*) '    F   = ', F * RCU
      WRITE(NFSTDO,*) '    ALK = ', ATIN(IJ, K)
      WRITE(NFSTDO,*) '    F - ALK = ', F * RCU - ATIN(IJ, K)
      WRITE(NFSTDO,*) '            < ALKERR = ', ALKERR * RCU
      WRITE(NFSTDO,*) '   [Info] '
     !WRITE(NFSTDO,*) '    at  I = ', MOD(IJ,NXDIM) - ISTR + 1
     !WRITE(NFSTDO,*) '        J = ', IJ/NXDIM - JSTR + 2
     !WRITE(NFSTDO,*) '        K = ', K - SKST + 1
      WRITE(NFSTDO,*) '    TEM = ', TEM(IJ)
      WRITE(NFSTDO,*) '    SAL = ', SAL(IJ)
      WRITE(NFSTDO,*) '    DIC = ', CTIN(IJ, K)

      XLW = PH(IJ, K) - DELT
      XHI = PH(IJ, K) + DELT
      DO ITR = 1, ITRMAX
        X = 0.5D0 * ( XLW + XHI )
        H = 10.0D0 ** ( - X )
        F = FALK( &
     &      H, &
     &      K1(IJ),  K2(IJ), &
     &      KB(IJ),  KW(IJ), &
     &      K1P(IJ), K2P(IJ),  K3P(IJ), &
     &      KSI(IJ), KSO4(IJ), KF(IJ), &
     &      CT(IJ, K),  BT(IJ),   PT(IJ), &
     &      SIT(IJ), SO4T(IJ), FLT(IJ) )
        IF ( ABS(F-AT(IJ, K)).LT.ALKERR ) GOTO 999
        IF ( ( F - AT(IJ, K) ).GT.0.D0 ) THEN
          XHI = X
        ELSE
          XLW = X
        END IF
      END DO

!     STOP

999 CONTINUE
      PH(IJ, K) = X

    END IF
    END DO
    END DO

    DO K = SKST, SKEN
    DO IJ = IJTSTR, IJTEND
   !IF ( TEM(IJ) /= 0.0D0 .and. SAL(IJ) /= 0.0D0 ) THEN
    IF ( TEM(IJ) /= rmiss .and. SAL(IJ) /= rmiss ) THEN
      X = 10.D0 ** ( - PH(IJ, K) )
!     PCO2O(IJ) = CT(IJ) *X*X &
!  &            / ( X*X + X*K1(IJ) + K1(IJ)*K2(IJ) ) &
!  &            / K0(IJ) &
!  &            * PM      !! [atm]-->[ppmv]
!     ALP(IJ) = K0(IJ)      !! [Mol/kg/atm] = [mMol/m^3/ppmv]
      ! -- CO2 --
      CO2(IJ, K)  = CT(IJ, K) *X*X &
   &              / ( X*X + X*K1(IJ) + K1(IJ)*K2(IJ) ) / CU
      ! -- BI-CARBONATE ION --
      HCO3(IJ, K) = K1(IJ) * CO2(IJ, K) / X
      ! -- CARBONATE ION --
      CO32(IJ, K) = K2(IJ) * HCO3(IJ, K) / X
    END IF
    END DO
    END DO

  end subroutine calc_digda_chemeq

end module module_digda_chemeq
