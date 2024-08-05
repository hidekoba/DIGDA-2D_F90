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

      subroutine calc_step_liquid( bw, rlt, phi, phii, dei, c, & !! IN
                                 & cn                          ) !! OUT

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
         eb(ij, k)     = eal(k) + ebe(k) + ega(k)
      end do
      end do

         ea(skst)      = - ebe(skst)
         ec(skst)      = - ega(skst)
         ea(sken)      = - ebe(sken)
         ec(sken)      = - ega(sken)
      do ij = ijtstr, ijtend
         eb(ij, skst)  = eal(skst) + ega(skst) &
     &   + ( phii(skst) * dei(skst) / dlzi(skst) / dlzi(skst) )
         eb(ij, sken)  = eal(sken) + ebe(sken)
      end do

      do ij = ijtstr, ijtend
         ed(ij, skst) = &
     &      + ( eal(skst) - ega(skst) ) * c(ij, skst) &
     &      + ega(skst) * c(ij, skst+1) &
     &      + rlt(ij, skst) &
     &   + ( phii(skst) * dei(skst) * bw(ij) / dlzi(skst) / dlzi(skst) )
      end do

      do k = skst+1, sken-1
      do ij = ijtstr, ijtend
         ed(ij, k) = &
     &      + ebe(k) * c(ij, k-1) &
     &      + ( eal(k) - ebe(k) - ega(k) ) * c(ij, k) &
     &      + ega(k) * c(ij, k+1) &
     &      + rlt(ij, k)
      end do
      end do

      do ij = ijtstr, ijtend
         ed(ij, sken) = &
     &      + ebe(sken) * c(ij, sken-1) &
     &      + ( eal(sken) - ega(sken) ) * c(ij, sken) &
     &      + rlt(ij, sken)
      end do

      call calc_tridiagonal_liquid( ea, eb, ec, ed, ex )

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
         eb(ij, k)     = eal(k) + ebe(k) + ega(k)
      end do
      end do

         ea(skst)      = - ebe(skst)
         ec(skst)      = - ega(skst)
         ea(sken)      = - ebe(sken)
         ec(sken)      = - ega(sken)
      do ij = ijtstr, ijtend
         eb(ij, skst)  = eal(skst) + ega(skst) &
     &   + ( phii(skst) * dei(skst) / dlzi(skst) / dlzi(skst) )
         eb(ij, sken)  = eal(sken) + ebe(sken)
      end do

      do k = skst, sken
      do ij = ijtstr, ijtend
         ed(ij, k) = &
     &      + eal(k) * c(ij, k) &
     &      + rlt(ij, k)
      end do
      end do

      do ij = ijtstr, ijtend
         ed(ij, skst) = ed(ij, skst) &
     &   + ( phii(skst) * dei(skst) * bw(ij) / dlzi(skst) / dlzi(skst) )
      end do

      call calc_tridiagonal_liquid( ea, eb, ec, ed, ex )

      end if

!     -- t = t + dt --

      do k = skst, sken
      do ij = ijtstr, ijtend
         cn(ij, k) = max( ex(ij, k), 0.0d0 )
         ex(ij, k) = 0.0d0
      end do
      end do

      end subroutine calc_step_liquid

! =====================================================================

      subroutine calc_step_liquid_o2( bw, rlt, phi, phii, dei, c, & !! IN
                                    & cn                          ) !! OUT

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
         if ( c(ij, k) /= 0.0d0 ) then
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
         eep(ij, k) = - 0.50d0 * rlt(ij, k)
      end do
      end do

      do k = skst+1, sken-1
         ea(k)         = - ebe(k)
         ec(k)         = - ega(k)
      do ij = ijtstr, ijtend
         eb(ij, k)     = eal(k) + ebe(k) + ega(k) + eep(ij, k)
      end do
      end do

         ea(skst)      = - ebe(skst)
         ec(skst)      = - ega(skst)
         ea(sken)      = - ebe(sken)
         ec(sken)      = - ega(sken)
      do ij = ijtstr, ijtend
         eb(ij, skst)  = eal(skst) + ega(skst) + eep(ij, skst) &
     &   + ( phii(skst) * dei(skst) / dlzi(skst) / dlzi(skst) )
         eb(ij, sken)  = eal(sken) + ebe(sken) + eep(ij, sken)
      end do

      do ij = ijtstr, ijtend
         ed(ij, skst) = &
     &      + ( eal(skst) - ega(skst) ) * c(ij, skst) &
     &      + ega(skst) * c(ij, skst+1) &
     &      + 0.50d0 * rlt(ij, skst) * c(ij, skst) &
     &   + ( phii(skst) * dei(skst) * bw(ij) / dlzi(skst) / dlzi(skst) )
      end do

      do k = skst+1, sken-1
      do ij = ijtstr, ijtend
         ed(ij, k) = &
     &      + ebe(k) * c(ij, k-1) &
     &      + ( eal(k) - ebe(k) - ega(k) ) * c(ij, k) &
     &      + ega(k) * c(ij, k+1) &
     &      + 0.50d0 * rlt(ij, k) * c(ij, k)
      end do
      end do

      do ij = ijtstr, ijtend
         ed(ij, sken) = &
     &      + ebe(sken) * c(ij, sken-1) &
     &      + ( eal(sken) - ega(sken) ) * c(ij, sken) &
     &      + 0.50d0 * rlt(ij, sken) * c(ij, sken)
      end do

      call calc_tridiagonal_liquid( ea, eb, ec, ed, ex )

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
         eep(ij, k) = - rlt(ij, k)
      end do
      end do

      do k = skst+1, sken-1
         ea(k)         = - ebe(k)
         ec(k)         = - ega(k)
      do ij = ijtstr, ijtend
         eb(ij, k)     = eal(k) + ebe(k) + ega(k) + eep(ij, k)
      end do
      end do

         ea(skst)      = - ebe(skst)
         ec(skst)      = - ega(skst)
         ea(sken)      = - ebe(sken)
         ec(sken)      = - ega(sken)
      do ij = ijtstr, ijtend
         eb(ij, skst)  = eal(skst) + ega(skst) + eep(ij, skst) &
     &   + ( phii(skst) * dei(skst) / dlzi(skst) / dlzi(skst) )
         eb(ij, sken)  = eal(sken) + ebe(sken) + eep(ij, sken)
      end do

      do k = skst, sken
      do ij = ijtstr, ijtend
         ed(ij, k) = &
     &      + eal(k) * c(ij, k)
      end do
      end do

      do ij = ijtstr, ijtend
         ed(ij, skst) = ed(ij, skst) &
     &   + ( phii(skst) * dei(skst) * bw(ij) / dlzi(skst) / dlzi(skst) )
      end do

      call calc_tridiagonal_liquid( ea, eb, ec, ed, ex )

      end if

!     -- t = t + dt --

      do k = skst, sken
      do ij = ijtstr, ijtend
         cn(ij, k) = max( ex(ij, k), 0.0d0 )
         ex(ij, k) = 0.0d0
      end do
      end do

      end subroutine calc_step_liquid_o2

! =====================================================================

      subroutine calc_tridiagonal_liquid( ea, eb, ec, ed, ex )

      double precision, intent(in)  :: ea(sdzdim), eb(nxy,sdzdim), ec(sdzdim)
      double precision, intent(in)  :: ed(nxy,sdzdim)
      double precision, intent(out) :: ex(nxy,sdzdim)

      double precision :: ew(nxy,sdzdim)
      double precision :: bet(nxy)

      integer :: ij, k

      do ij = ijtstr, ijtend
         bet(ij)      = eb(ij, skst)
         ex(ij, skst) = ed(ij, skst) / bet(ij)
      end do

      do k = skst+1, sken
      do ij = ijtstr, ijtend
         ew(ij, k) = ec(k-1) / bet(ij)
         bet(ij)   = eb(ij, k) - ea(k) * ew(ij, k)
         ex(ij, k) = (ed(ij, k) - ea(k) * ex(ij, k-1) ) / bet(ij)
      end do
      end do

      do k = sken-1, skst, -1
      do ij = ijtstr, ijtend
         ex(ij, k) = ex(ij, k) - ew(ij, k+1) * ex(ij, k+1)
      end do
      end do

      end subroutine calc_tridiagonal_liquid

end module module_digda_step_liquid
