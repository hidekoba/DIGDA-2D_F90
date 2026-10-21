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
  use module_digda_timestep
  implicit none

  double precision :: dflx(nxy, sdzdim+1)

  double precision :: ea(sdzdim), eb(nxy, sdzdim), ec(sdzdim)
  double precision :: eal(sdzdim), ebe(sdzdim), ega(sdzdim), eep(nxy, sdzdim)
  double precision :: ed(nxy, sdzdim),  ex(nxy, sdzdim)

  logical :: cdcn = .false. !! Crank nicolson ?
  logical :: cdii = .true.  !! Implicit ?
  namelist /nmsdli/ cdcn, cdii

contains

  logical function is_missing_value(val)
    double precision, intent(in) :: val
    is_missing_value = (abs(val - rmiss) < 1.0d-6) .or. &
     &                 (abs(val - 998.9999d0) < 1.0d-6) .or. &
     &                 (.not. (val == val)) .or. &
     &                 (abs(val) > 1.0d30)
  end function is_missing_value

end module globals_digda_step_liquid

! =====================================================================

module module_digda_step_liquid
  use globals_digda_step_liquid
  implicit none
contains

  subroutine set_digda_step_liquid

    integer :: ij, k

    write(nfstdo, *) '*** Sediment : set_digda_step_liquid ***'
    rewind(nfparm)
    read(nfparm, nmsdli, end=99)
99 continue
    write(nfstdo, nmsdli)

  end subroutine set_digda_step_liquid

! =====================================================================

      subroutine calc_step_liquid( bw, rlt, phi, phii, dei, c, &
                                 & cn                          )

      double precision, intent(in)  :: bw(nxy)
      double precision, intent(in)  :: rlt(nxy,sdzdim)
      double precision, intent(in)  :: phi(sdzdim)
      double precision, intent(in)  :: phii(sdzdim+1)
      double precision, intent(in)  :: dei(sdzdim+1)
      double precision, intent(in)  :: c(nxy,sdzdim)
      double precision, intent(out) :: cn(nxy,sdzdim)

      integer :: ij, k

!     -- Diffusion : Crank nicolson --

      if (cdcn) then

      do k = skst, sken
         eal(k) = ( phi(k) / tssec )
         ebe(k) = ( ( phii(k)   * dei(k)   ) &
     &            / ( 2.0d0 * dlz(k) * dlzi(k)   ) )
         ega(k) = ( ( phii(k+1) * dei(k+1) ) &
     &            / ( 2.0d0 * dlz(k) * dlzi(k+1) ) )
      end do

      do k = skst+1, sken-1
         ea(k)         = - ebe(k)
         ec(k)         = - ega(k)
      do ij = ijtstr, ijtend
         if ( .not. is_missing_value(bw(ij)) ) then
            eb(ij, k) = eal(k) + ebe(k) + ega(k)
         end if
      end do
      end do

         ea(skst)      = - ebe(skst)
         ec(skst)      = - ega(skst)
         ea(sken)      = - ebe(sken)
         ec(sken)      = - ega(sken)
      do ij = ijtstr, ijtend
         if ( .not. is_missing_value(bw(ij)) ) then
            eb(ij, skst) = eal(skst) + ega(skst) &
     &         + ( phii(skst) * dei(skst) / dlzi(skst) / dlzi(skst) )
            eb(ij, sken) = eal(sken) + ebe(sken)
         end if
      end do

      do ij = ijtstr, ijtend
         if ( .not. is_missing_value(bw(ij)) ) then
         ed(ij, skst) = &
     &      + ( eal(skst) - ega(skst) ) * c(ij, skst) &
     &      + ega(skst) * c(ij, skst+1) &
     &      + rlt(ij, skst) &
     &   + ( phii(skst) * dei(skst) * bw(ij) / dlzi(skst) / dlzi(skst) )
         end if
      end do

      do k = skst+1, sken-1
      do ij = ijtstr, ijtend
         if ( .not. is_missing_value(bw(ij)) ) then
         ed(ij, k) = &
     &      + ebe(k) * c(ij, k-1) &
     &      + ( eal(k) - ebe(k) - ega(k) ) * c(ij, k) &
     &      + ega(k) * c(ij, k+1) &
     &      + rlt(ij, k)
         end if
      end do
      end do

      do ij = ijtstr, ijtend
         if ( .not. is_missing_value(bw(ij)) ) then
         ed(ij, sken) = &
     &      + ebe(sken) * c(ij, sken-1) &
     &      + ( eal(sken) - ega(sken) ) * c(ij, sken) &
     &      + rlt(ij, sken)
         end if
      end do

      call calc_tridiagonal_liquid( bw, ea, eb, ec, ed, ex )

      end if

!     -- Diffusion : Implicit --

      if (cdii) then

      do k = skst, sken
         eal(k) = ( phi(k) / tssec )
         ebe(k) = ( ( phii(k)   * dei(k)   ) &
     &            / ( dlz(k) * dlzi(k)   ) )
         ega(k) = ( ( phii(k+1) * dei(k+1) ) &
     &            / ( dlz(k) * dlzi(k+1) ) )
      end do

      do k = skst+1, sken-1
         ea(k)         = - ebe(k)
         ec(k)         = - ega(k)
      do ij = ijtstr, ijtend
         if ( .not. is_missing_value(bw(ij)) ) then
         eb(ij, k)     = eal(k) + ebe(k) + ega(k)
         end if
      end do
      end do

         ea(skst)      = - ebe(skst)
         ec(skst)      = - ega(skst)
         ea(sken)      = - ebe(sken)
         ec(sken)      = - ega(sken)
      do ij = ijtstr, ijtend
         if ( .not. is_missing_value(bw(ij)) ) then
         eb(ij, skst)  = eal(skst) + ega(skst) &
     &   + ( phii(skst) * dei(skst) / dlzi(skst) / dlzi(skst) )
         eb(ij, sken)  = eal(sken) + ebe(sken)
         end if
      end do

      do k = skst, sken
      do ij = ijtstr, ijtend
         if ( .not. is_missing_value(bw(ij)) ) then
         ed(ij, k) = &
     &      + eal(k) * c(ij, k) &
     &      + rlt(ij, k)
         end if
      end do
      end do

      do ij = ijtstr, ijtend
         if ( .not. is_missing_value(bw(ij)) ) then
         ed(ij, skst) = ed(ij, skst) &
     &   + ( phii(skst) * dei(skst) * bw(ij) / dlzi(skst) / dlzi(skst) )
         end if
      end do

      call calc_tridiagonal_liquid( bw, ea, eb, ec, ed, ex )

      end if

!     -- t = t + dt --

      do k = skst, sken
      do ij = ijtstr, ijtend
         if ( .not. is_missing_value(bw(ij)) ) then
            cn(ij, k) = max( ex(ij, k), 0.0d0 )
            ex(ij, k) = 0.0d0
         else
            cn(ij, k) = rmiss
            ex(ij, k) = 0.0d0
         end if
      end do
      end do

      end subroutine calc_step_liquid

! =====================================================================

      subroutine calc_step_liquid_o2( bw, rlt, phi, phii, dei, c, &
                                    & cn                          )

      double precision, intent(in)  :: bw(nxy)
      double precision, intent(inout) :: rlt(nxy,sdzdim)
      double precision, intent(in)  :: phi(sdzdim)
      double precision, intent(in)  :: phii(sdzdim+1)
      double precision, intent(in)  :: dei(sdzdim+1)
      double precision, intent(in)  :: c(nxy,sdzdim)
      double precision, intent(out) :: cn(nxy,sdzdim)

      integer :: ij, k

      do k = skst, sken
      do ij = ijtstr, ijtend
         if ( .not. (c(ij, k) == 0.0d0) .and. &
     &       (c(ij, k) == c(ij, k)) .and. &
     &       (abs(c(ij, k)) < 1.0d30) ) then
            rlt(ij, k) = rlt(ij, k) / c(ij, k)
         else
            rlt(ij, k) = 0.0d0
         end if
      end do
      end do

!     -- Diffusion : Crank nicolson --

      if (cdcn) then

      do k = skst, sken
         eal(k) = ( phi(k) / tssec )
         ebe(k) = ( ( phii(k)   * dei(k)   ) &
     &            / ( 2.0d0 * dlz(k) * dlzi(k)   ) )
         ega(k) = ( ( phii(k+1) * dei(k+1) ) &
     &            / ( 2.0d0 * dlz(k) * dlzi(k+1) ) )
      do ij = ijtstr, ijtend
         if ( .not. is_missing_value(bw(ij)) ) then
            eep(ij, k) = - 0.50d0 * rlt(ij, k)
         end if
      end do
      end do

      do k = skst+1, sken-1
         ea(k)         = - ebe(k)
         ec(k)         = - ega(k)
      do ij = ijtstr, ijtend
         if ( .not. is_missing_value(bw(ij)) ) then
            eb(ij, k) = eal(k) + ebe(k) + ega(k) + eep(ij, k)
         end if
      end do
      end do

         ea(skst)      = - ebe(skst)
         ec(skst)      = - ega(skst)
         ea(sken)      = - ebe(sken)
         ec(sken)      = - ega(sken)
      do ij = ijtstr, ijtend
         if ( .not. is_missing_value(bw(ij)) ) then
            eb(ij, skst) = eal(skst) + ega(skst) + eep(ij, skst) &
     &         + ( phii(skst) * dei(skst) / dlzi(skst) / dlzi(skst) )
            eb(ij, sken) = eal(sken) + ebe(sken) + eep(ij, sken)
         end if
      end do

      do ij = ijtstr, ijtend
         if ( .not. is_missing_value(bw(ij)) ) then
            ed(ij, skst) = &
     &         + ( eal(skst) - ega(skst) ) * c(ij, skst) &
     &         + ega(skst) * c(ij, skst+1) &
     &         + 0.50d0 * rlt(ij, skst) * c(ij, skst) &
     &         + ( phii(skst) * dei(skst) * bw(ij) / dlzi(skst) / dlzi(skst) )
         end if
      end do

      do k = skst+1, sken-1
      do ij = ijtstr, ijtend
         if ( .not. is_missing_value(bw(ij)) ) then
            ed(ij, k) = &
     &         + ebe(k) * c(ij, k-1) &
     &         + ( eal(k) - ebe(k) - ega(k) ) * c(ij, k) &
     &         + ega(k) * c(ij, k+1) &
     &         + 0.50d0 * rlt(ij, k) * c(ij, k)
         end if
      end do
      end do

      do ij = ijtstr, ijtend
         if ( .not. is_missing_value(bw(ij)) ) then
            ed(ij, sken) = &
     &         + ebe(sken) * c(ij, sken-1) &
     &         + ( eal(sken) - ega(sken) ) * c(ij, sken) &
     &         + 0.50d0 * rlt(ij, sken) * c(ij, sken)
         end if
      end do

      call calc_tridiagonal_liquid( bw, ea, eb, ec, ed, ex )

      end if

!     -- Diffusion : Implicit --

      if (cdii) then

      do k = skst, sken
         eal(k) = ( phi(k) / tssec )
         ebe(k) = ( ( phii(k)   * dei(k)   ) &
     &            / ( dlz(k) * dlzi(k)   ) )
         ega(k) = ( ( phii(k+1) * dei(k+1) ) &
     &            / ( dlz(k) * dlzi(k+1) ) )
      do ij = ijtstr, ijtend
         if ( .not. is_missing_value(bw(ij)) ) then
            eep(ij, k) = - rlt(ij, k)
         end if
      end do
      end do

      do k = skst+1, sken-1
         ea(k)         = - ebe(k)
         ec(k)         = - ega(k)
      do ij = ijtstr, ijtend
         if ( .not. is_missing_value(bw(ij)) ) then
            eb(ij, k) = eal(k) + ebe(k) + ega(k) + eep(ij, k)
         end if
      end do
      end do

         ea(skst)      = - ebe(skst)
         ec(skst)      = - ega(skst)
         ea(sken)      = - ebe(sken)
         ec(sken)      = - ega(sken)
      do ij = ijtstr, ijtend
         if ( .not. is_missing_value(bw(ij)) ) then
            eb(ij, skst) = eal(skst) + ega(skst) + eep(ij, skst) &
     &         + ( phii(skst) * dei(skst) / dlzi(skst) / dlzi(skst) )
            eb(ij, sken) = eal(sken) + ebe(sken) + eep(ij, sken)
         end if
      end do

      do k = skst, sken
      do ij = ijtstr, ijtend
         if ( .not. is_missing_value(bw(ij)) ) then
            ed(ij, k) = &
     &         + eal(k) * c(ij, k)
         end if
      end do
      end do

      do ij = ijtstr, ijtend
         if ( .not. is_missing_value(bw(ij)) ) then
            ed(ij, skst) = ed(ij, skst) &
     &         + ( phii(skst) * dei(skst) * bw(ij) / dlzi(skst) / dlzi(skst) )
         end if
      end do

      call calc_tridiagonal_liquid( bw, ea, eb, ec, ed, ex )

      end if

!     -- t = t + dt --

      do k = skst, sken
      do ij = ijtstr, ijtend
         if ( .not. is_missing_value(bw(ij)) ) then
            cn(ij, k) = max( ex(ij, k), 0.0d0 )
            ex(ij, k) = 0.0d0
         else
            cn(ij, k) = rmiss
            ex(ij, k) = 0.0d0
         end if
      end do
      end do

      end subroutine calc_step_liquid_o2

! =====================================================================

      subroutine calc_tridiagonal_liquid( bw, ea, eb, ec, ed, ex )

      double precision, intent(in)  :: ea(sdzdim), eb(nxy,sdzdim), ec(sdzdim)
      double precision, intent(in)  :: ed(nxy,sdzdim)
      double precision, intent(in)  :: bw(nxy)
      double precision, intent(out) :: ex(nxy,sdzdim)

      double precision :: ew(nxy,sdzdim)
      double precision :: bet(nxy)

      integer :: ij, k

      do ij = ijtstr, ijtend
         if ( .not. is_missing_value(bw(ij)) ) then
            bet(ij) = eb(ij, skst)
            if ( (bet(ij) == bet(ij)) .and. (abs(bet(ij)) > 1.0d-30) ) then
               ex(ij, skst) = ed(ij, skst) / bet(ij)
            else
               ex(ij, skst) = rmiss
            end if
         end if
      end do

      do k = skst+1, sken
      do ij = ijtstr, ijtend
         if ( .not. is_missing_value(bw(ij)) ) then
            if ( (bet(ij) == bet(ij)) .and. (abs(bet(ij)) > 1.0d-30) ) then
               ew(ij, k) = ec(k-1) / bet(ij)
               bet(ij)   = eb(ij, k) - ea(k) * ew(ij, k)
               if ( (bet(ij) == bet(ij)) .and. (abs(bet(ij)) > 1.0d-30) ) then
                  ex(ij, k) = (ed(ij, k) - ea(k) * ex(ij, k-1) ) / bet(ij)
               else
                  ex(ij, k) = rmiss
               end if
            else
               ex(ij, k) = rmiss
            end if
         else
            ex(ij, k) = rmiss
         end if
      end do
      end do

      do k = sken-1, skst, -1
      do ij = ijtstr, ijtend
         if ( .not. is_missing_value(bw(ij)) ) then
            ex(ij, k) = ex(ij, k) - ew(ij, k+1) * ex(ij, k+1)
         else
            ex(ij, k) = rmiss
         end if
      end do
      end do

      end subroutine calc_tridiagonal_liquid

end module module_digda_step_liquid
