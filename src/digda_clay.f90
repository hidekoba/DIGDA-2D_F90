!######################################################################

! --- information -----------------------------------------------------
!
!  Sediment model : Clay
!
!  HISTORY
!     '15.01.09  H.Kobayashi :
!
! ---------------------------------------------------------------------

module globals_digda_clay
  use module_zcfile
  use module_digda_grid
  implicit none

end module globals_digda_clay

! =====================================================================

module module_digda_clay
  use globals_digda_clay
  implicit none
contains

  subroutine set_digda_clay

    write(nfstdo, *) '*** Sediment : set_digda_clay ***'

  end subroutine set_digda_clay

  subroutine calc_clay_dissolution( scla, &
                                  & rcla  )

    double precision, intent(in) :: scla(nxy, sdzdim)
    double precision, intent(out) :: rcla(nxy, sdzdim)

    integer :: ij, k

    do k = skst, sken
      do ij = ijtstr, ijtend
        rcla(ij, k)  = 0.0d0
      end do
    end do

  end subroutine calc_clay_dissolution

end module module_digda_clay
