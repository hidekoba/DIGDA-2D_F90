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
  use module_digda_timestep
  implicit none

  double precision :: kox = 2.00d-9  !! 2.00d-9 [/s]
  double precision :: kax = 7.90d-11 !! 7.90d-11 [/s]
  namelist /nmsdpm/ kox, kax

end module globals_digda_pom

! =====================================================================

module module_digda_pom
  use globals_digda_pom
  implicit none
contains

  subroutine set_digda_pom

    write(nfstdo, *) '*** Sediment : set_digda_pom ***'
    rewind(nfparm)
    read(nfparm, nmsdpm, end=99)
99 continue
    write(nfstdo, nmsdpm)

  end subroutine set_digda_pom

  subroutine calc_pom_dissolution( spom, cdo2, cp, op, phi, &
                                 & rpom                     )

    double precision, intent(in) :: spom(nxy, sdzdim)
    double precision, intent(in) :: cdo2(nxy, sdzdim)
    double precision, intent(in) :: cp, op
    double precision, intent(in) :: phi(sdzdim)
    double precision, intent(out) :: rpom(nxy, sdzdim)

    double precision :: rpommax

    integer :: ij, k

    do k = skst, sken
      do ij = ijtstr, ijtend
        if ( spom(ij, k) /= 0.0d0 ) then
          if ( cdo2(ij, k) > 0.0d0 ) then
            rpom(ij, k) = kox
          else if ( cdo2(ij, k) == 0.0d0 ) then
            rpom(ij, k) = kax
          end if
          rpommax = ( phi(k) / ( 1.0d0 - phi(k) ) ) &
         &        * ( cp / op ) * cdo2(ij, k) / spom(ij, k) / tssec
          rpom(ij, k) = max( min( rpom(ij, k), rpommax ), kax )
        else
          rpom(ij, k) = 0.0d0
        end if
      end do
    end do

  end subroutine calc_pom_dissolution
end module module_digda_pom
