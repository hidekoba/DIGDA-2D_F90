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

  ! Constants for carbonate processes [Millero, 1995]
  double precision :: kprime = 4.750d-7 !! 4.750d-7 [mol2 / kg2]
  double precision :: delv = -44.0      !! -44.0    [cm3 / mol]
  double precision :: rr = 83.14        !! 83.14    [cm3 bar / K mol]
  double precision :: dk = -0.0133      !! -0.0133  [cm3 / bar mol]

  ! Model parameters
  double precision :: rccar = 365.d0    !! Reaction rate constant [1/day]
  double precision :: ncar = 4.5d0      !! Reaction order
  namelist /nmsdca/ rccar, ncar

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

    rccar  = rccar / 3.1536d+7  !! Convert from [/yr] to [/s]

  end subroutine set_digda_carbonate

  subroutine calc_carbonate_dissolution( scar, cco32, satco32, phi, &
                                       & rcar                       )

    double precision, intent(in) :: scar(nxy, sdzdim)
    double precision, intent(in) :: cco32(nxy, sdzdim)
    double precision, intent(in) :: satco32(nxy)
    double precision, intent(in) :: phi(sdzdim)
    double precision, intent(out) :: rcar(nxy, sdzdim)

    integer :: ij, k
    logical :: valid_scar, undersaturated

    do k = skst, sken
      do ij = ijtstr, ijtend
        valid_scar = scar(ij, k) /= 0.0d0
        undersaturated = cco32(ij, k) < satco32(ij)
        if ( valid_scar .and. undersaturated ) then
          rcar(ij, k) = rccar &
   &        * ( ( 1.0d0 - ( cco32(ij, k) / satco32(ij) ) ) &
   &        ** ncar )
        else
          rcar(ij, k) = 0.0d0
        end if
      end do
    end do

  end subroutine calc_carbonate_dissolution

  subroutine calc_carbonate_saturation( temp, salt, bdep, &
                                      & csat              )

    double precision, intent(in) :: temp(nxy), salt(nxy)
    double precision, intent(in) :: bdep(nxy)
    double precision, intent(out) :: csat(nxy)

    integer :: ij, k
    double precision :: pres, t_kelvin, kpres

    do ij = ijtstr, ijtend
      t_kelvin = temp(ij) + kelvin
      pres  =  bdep(ij) / 1.0d+3  !! Convert depth to pressure [bar]

      kpres = log(kprime) &
   &        - delv / ( rr * ( t_kelvin ) ) * (pres) &
   &        + 0.50d0 * dk / ( rr * ( t_kelvin ) ) * (pres)**2
      kpres = exp(kpres)
      csat(ij)  = ( kpres / 1.028d-2 ) * 1.0d+6 !! Convert to [mMol/m3]
    end do

  end subroutine calc_carbonate_saturation

end module module_digda_carbonate
