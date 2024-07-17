!######################################################################

! --- information -----------------------------------------------------
!
!  Sediment model : CaCO3
!
!  HISTORY
!     '15.01.09  H.Kobayashi :
!
! ---------------------------------------------------------------------

module globals_digda_carbonate
  implicit none

  double precision :: RCCAR = 365.D0 !! 1.0 [1/day]
  double precision :: NCAR = 4.5D0
  namelist /NMSDCA/ RCCAR, NCAR

 ! Millero [1995]
  double precision :: KPRIME = 4.750D-7 !! 4.750D-7 [mol2 / kg2]
  double precision :: DELV = -4.400D+1  !! -4.400D+1 [cm3 / mol]
  double precision :: RR = 8.314D+1     !! 8.314D+1 [cm3 bar / K mol]
  double precision :: DK = -1.330D-2    !! -1.330D-2 [cm3 / bar mol]
  double precision :: PRES, KPRES

end module globals_digda_carbonate

! =====================================================================

module module_digda_carbonate
  use globals_digda_carbonate
  use module_zcfile
  use module_zcphys
  use module_digda_grid
  implicit none
contains

  subroutine set_digda_carbonate

    WRITE(NFSTDO, *) '*** Sediment : set_digda_carbonate ***'
    REWIND(NFPARM)
    READ (NFPARM, NMSDCA, END=99)
99 continue
    WRITE(NFSTDO, NMSDCA)

    RCCAR  = RCCAR / 3.1536D+7  !! [/yr] --> [/s]

  end subroutine set_digda_carbonate

  subroutine calc_carbonate_dissolution( &
 &  SCAR, CCO32, SATCO32, PHI, &
 &  RCAR )

    double precision, intent(in) :: SCAR(NXY, SDZDIM)
    double precision, intent(in) :: CCO32(NXY, SDZDIM)
    double precision, intent(in) :: SATCO32(NXY)
    double precision, intent(in) :: PHI(SDZDIM)
    double precision, intent(out) :: RCAR(NXY, SDZDIM)

    integer :: IJ, K

    DO K = SKST, SKEN
      DO IJ = IJTSTR, IJTEND
        IF ( SCAR(IJ, K) /= 0.0D0 .and. &
   &         CCO32(IJ, K) < SATCO32(IJ) ) THEN
          RCAR(IJ, K) = RCCAR &
   &        * ( ( 1.0D0 - ( CCO32(IJ, K) / SATCO32(IJ) ) ) &
   &        ** NCAR )
        ELSE
          RCAR(IJ, K) = 0.0D0
        END IF
      END DO
    END DO

  end subroutine calc_carbonate_dissolution

  subroutine calc_carbonate_saturation( &
&   TEM, SAL, BDEP,&
&   CSAT )

    double precision, intent(in) :: TEM(NXY), SAL(NXY)
    double precision, intent(in) :: BDEP(NXY)
    double precision, intent(out) :: CSAT(NXY)

    integer :: IJ, K

    DO IJ = IJTSTR, IJTEND
      PRES  =  BDEP(IJ) / 1.0D+3  !! [cm] --> [bar]
      KPRES = LOG(KPRIME) &
   &        - DELV / ( RR * ( TEM(IJ) + kelvin ) )        * (PRES)    & !! Effect of Pressure
   &        + 0.50D0 * DK / ( RR * ( TEM(IJ) + kelvin ) ) * (PRES)**2   !! Effect of Pressure
      KPRES = EXP(KPRES)
      CSAT(IJ)  = ( KPRES / 1.028D-2 ) * 1.0D+6                         !! Millero [1982]
    END DO

  end subroutine calc_carbonate_saturation

end module module_digda_carbonate
