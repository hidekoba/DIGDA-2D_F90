!######################################################################

! --- information -----------------------------------------------------
!
!  Sediment model : Opal
!
!  HISTORY
!     '15.01.09  H.Kobayashi :
!
! ---------------------------------------------------------------------

module globals_digda_silicate
  use module_zcfile
  use module_digda_grid
  implicit none

  double precision :: rcopa = 5.00d-3 !! 5.0*10-3 [cm3/mol/s]
  logical :: roprid = .false. !! Ridgwell's DEI formulation ?
  logical :: roparc = .true.  !! Archer's DEI formulation ?
  namelist /nmsdop/ rcopa, roprid, roparc

end module globals_digda_silicate

! =====================================================================

module module_digda_silicate
  use globals_digda_silicate
  implicit none
contains

  subroutine set_digda_silicate

    write(nfstdo, *) '*** Sediment : set_digda_silicate ***'
    rewind(nfparm)
    read(nfparm, nmsdop, end=99)
99 continue
    write(nfstdo, nmsdop)

    rcopa  = rcopa * 1.0d-9  !! [cm3/mol/s] --> [m3/mMol/s]

  end subroutine set_digda_silicate

  subroutine calc_silicate_dissolution( sopa, csil, pcla, popa, temp, &
                                      & ropa, satsil                  )

    double precision, intent(in) :: sopa(nxy, sdzdim)
    double precision, intent(in) :: csil(nxy, sdzdim)
    double precision, intent(in) :: pcla(nxy, sdzdim)
    double precision, intent(in) :: popa(nxy, sdzdim)
    double precision, intent(in) :: temp(nxy)
    double precision, intent(out) :: ropa(nxy, sdzdim)
    double precision, intent(out) :: satsil(nxy, sdzdim)

    double precision :: gammal(nxy, sdzdim)
    double precision :: pco(nxy, sdzdim)
    double precision :: uo(nxy, sdzdim)
    double precision :: eta1(nxy, sdzdim), eta2(nxy, sdzdim)
    double precision :: ko(nxy, sdzdim)

    integer :: ij, k

    do k = skst, sken
      do ij = ijtstr, ijtend
        if ( sopa(ij, k) == 0.0d0 ) then
          gammal(ij, k) = 0.2d0
        else
          pco(ij, k) = pcla(ij, k) / popa(ij, k)
          if ( pco(ij, k) > 1.5d+1 ) then
            gammal(ij, k) = 0.2d0
          else
            gammal(ij, k) = 1.0d0 - ( ( 4.5d-2 * pco(ij, k) ) &
   &                      ** 5.8d-1 )
          end if
        end if
        satsil(ij, k) = ( gammal(ij, k) * 1.0d+1 &
   &      ** ( 6.44d0 - ( 9.68d+2 / ( 2.7315d+2 + temp(ij) ) ) ) )
      end do
    end do

    do k = skst, sken
      do ij = ijtstr, ijtend
        if ( sopa(ij, k) /= 0.0d0 .and. &
   &         csil(ij, k) < satsil(ij, k) ) then
          if (roprid) then
            uo(ij, k)   = 1.0d0 - ( csil(ij, k) / satsil(ij, k) )
            eta1(ij, k) = 2.25d-1 * ( 1.0d0 + ( temp(ij) / 1.5d+1 ) ) &
   &                    * uo(ij, k) +  7.75d-1 &
   &                    * ( ( ( ( 1.0d0 + ( temp(ij) / 4.0d+2 ) ) &
   &                    ** 4.0d0 ) * uo(ij, k) ) ** 9.25d0 )
            eta2(ij, k) = 2.6d-1 + 7.4d-1 &
   &                    * exp( - lz(k) / 7.0d-2 )
            ko(ij, k)   = 5.0d-2 + 5.5d-2 &
   &                    * ( ( 1.64d-2 + pco(ij, k) ) &
   &                    ** ( - 7.5d-1 ) )
            ropa(ij, k) = eta1(ij, k) * eta2(ij, k) * ko(ij, k) &
   &                    / 3.1536d+7
          else if (roparc) then
            ropa(ij, k) = rcopa &
   &          * ( satsil(ij, k) - csil(ij, k) )
          end if
        else
          ropa(ij, k) = 0.0d0
        end if
      end do
    end do

  end subroutine calc_silicate_dissolution

end module module_digda_silicate
