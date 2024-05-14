!######################################################################

! --- information -----------------------------------------------------
!
!  Sediment model : CaCO3
!
!  HISTORY
!     '15.01.09  H.Kobayashi :
!
! ---------------------------------------------------------------------

module globals_sdcaco3
      use zcfile_module
      use zcphys_module
      use sddim_module
      implicit none

!     REAL*8 SCAR(NXY, SDZDIM)
!     REAL*8 CCO32(NXY, SDZDIM)
!     REAL*8 SATCO32(NXY)
!     REAL*8 RCAR(NXY, SDZDIM)

      REAL*8 RCCAR, NCAR
      SAVE   RCCAR, NCAR
      NAMELIST /NMSDCA/ RCCAR, NCAR
      DATA   RCCAR, NCAR &
     &    / 365.D0, 4.50D0 /
          !! 1.0 [/day], 4.5

!     REAL*8 TEM(NXY), SAL(NXY), BDEP(NXY), CSAT(NXY)
      REAL*8 KPRIME, DELV, RR, DK, PRES, KPRES !! Millero [1995]
      SAVE   KPRIME, DELV, RR, DK
      DATA KPRIME &
     &    /  4.750D-7 /
          !! [mol2 / kg2]
      DATA DELV &
     &    / -4.400D+1 /
          !! [cm3 / mol]
      DATA RR &
     &    /  8.314D+1 /
          !! [cm3 bar / K mol]
      DATA DK &
     &    / -1.330D-2 /
          !! [cm3 / bar mol]

!     REAL*8 PHI(SDZDIM)

end module globals_sdcaco3

! =====================================================================

module sdcaco3_module
      use globals_sdcaco3
      implicit none
contains

      SUBROUTINE SDCACO3SET

      WRITE(NFSTDO, *) '*** Carbon : SDCACO3SET ***'
      REWIND(NFPARM)
      READ (NFPARM, NMSDCA, END=99)
 99   WRITE(NFSTDO, NMSDCA)

      RCCAR  = RCCAR / 3.1536D+7  !! [/yr] --> [/s]

      END SUBROUTINE SDCACO3SET

      SUBROUTINE DISCAR( &
     &      SCAR, CCO32, SATCO32, PHI, &
     &      RCAR )

      double precision, intent(in) :: SCAR(NXY, SDZDIM)
      double precision, intent(in) :: CCO32(NXY, SDZDIM)
      double precision, intent(in) :: SATCO32(NXY)
      double precision, intent(in) :: PHI(SDZDIM)
      double precision, intent(out) :: RCAR(NXY, SDZDIM)

      integer :: IJ, K

      DO K = SKST, SKEN
      DO IJ = IJTSTR, IJTEND
         IF ( SCAR(IJ, K) /= 0.0D0 .and. &
     &        CCO32(IJ, K) < SATCO32(IJ) ) THEN
            RCAR(IJ, K) = RCCAR &
     &          * ( ( 1.0D0 - ( CCO32(IJ, K) / SATCO32(IJ) ) ) &
     &          ** NCAR )
         ELSE
            RCAR(IJ, K) = 0.0D0
         END IF
      END DO
      END DO

      END SUBROUTINE DISCAR

      SUBROUTINE SATCAR( &
     &      TEM, SAL, BDEP,&
     &      CSAT )

      double precision, intent(in) :: TEM(NXY), SAL(NXY)
      double precision, intent(in) :: BDEP(NXY)
      double precision, intent(out) :: CSAT(NXY)

      integer :: IJ, K

      DO IJ = IJTSTR, IJTEND
         PRES  =  BDEP(IJ) / 1.0D+3  !! [cm] --> [bar]
         KPRES = LOG(KPRIME) &
     &      - DELV / ( RR * ( TEM(IJ) + kelvin ) )        * (PRES)     & !! Effect of Pressure
     &      + 0.50D0 * DK / ( RR * ( TEM(IJ) + kelvin ) ) * (PRES)**2    !! Effect of Pressure
         KPRES = EXP(KPRES)
         CSAT(IJ)  = ( KPRES / 1.028D-2 ) * 1.0D+6                     !! Millero [1982]
      END DO

      END SUBROUTINE SATCAR

end module sdcaco3_module
