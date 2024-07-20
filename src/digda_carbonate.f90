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

  double precision :: rccar = 365.d0 !! 1.0 [1/day]
  double precision :: ncar = 4.5d0
  namelist /nmsdca/ rccar, ncar

 ! Millero [1995]
  double precision :: kprime = 4.750d-7 !! 4.750d-7 [mol2 / kg2]
  double precision :: delv = -4.400d+1  !! -4.400d+1 [cm3 / mol]
  double precision :: rr = 8.314d+1     !! 8.314d+1 [cm3 bar / K mol]
  double precision :: dk = -1.330d-2    !! -1.330d-2 [cm3 / bar mol]
  double precision :: pres, kpres

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

    write(nfstdo, *) '*** Sediment : set_digda_carbonate ***'
    rewind(nfparm)
    read(nfparm, nmsdca, end=99)
99 continue
    write(nfstdo, nmsdca)

    rccar  = rccar / 3.1536d+7  !! [/yr] --> [/s]

  end subroutine set_digda_carbonate

  subroutine calc_carbonate_dissolution( scar, cco32, satco32, phi, & !! IN
                                       & rcar                       ) !! OUT

    double precision, intent(in) :: scar(nxy, sdzdim)
    double precision, intent(in) :: cco32(nxy, sdzdim)
    double precision, intent(in) :: satco32(nxy)
    double precision, intent(in) :: phi(sdzdim)
    double precision, intent(out) :: rcar(nxy, sdzdim)

    integer :: ij, k

    do k = skst, sken
      do ij = ijtstr, ijtend
        if ( scar(ij, k) /= 0.0d0 .and. &
   &         cco32(ij, k) < satco32(ij) ) then
          rcar(ij, k) = rccar &
   &        * ( ( 1.0d0 - ( cco32(ij, k) / satco32(ij) ) ) &
   &        ** ncar )
        else
          rcar(ij, k) = 0.0d0
        end if
      end do
    end do

  end subroutine calc_carbonate_dissolution

  subroutine calc_carbonate_saturation( temp, salt, bdep, & !! IN
                                      & csat              ) !! OUT

    double precision, intent(in) :: temp(nxy), salt(nxy)
    double precision, intent(in) :: bdep(nxy)
    double precision, intent(out) :: csat(nxy)

    integer :: ij, k

    do ij = ijtstr, ijtend
      pres  =  bdep(ij) / 1.0d+3  !! [cm] --> [bar]
      kpres = log(kprime) &
   &        - delv / ( rr * ( temp(ij) + kelvin ) )        * (pres)    & !! Effect of Pressure
   &        + 0.50d0 * dk / ( rr * ( temp(ij) + kelvin ) ) * (pres)**2   !! Effect of Pressure
      kpres = exp(kpres)
      csat(ij)  = ( kpres / 1.028d-2 ) * 1.0d+6                         !! Millero [1982]
    end do

  end subroutine calc_carbonate_saturation

end module module_digda_carbonate
