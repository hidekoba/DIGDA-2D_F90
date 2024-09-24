!######################################################################

! --- information -----------------------------------------------------
!
!  Sediment model : Particles
!
!  HISTORY
!     '15.01.09  H.Kobayashi :
!
! ---------------------------------------------------------------------

module globals_digda_step_solid
  use module_zcfile
  use module_digda_grid
  use module_digda_timestep
  implicit none

  double precision :: si(nxy, sdzdim+1)

  double precision :: c0(nxy, sdzdim+1)
  double precision :: c1(nxy, sdzdim+1)
  double precision :: c2(nxy, sdzdim+1)
  double precision :: flx(nxy, sdzdim+1), dflx(nxy, sdzdim+1)

  double precision :: ea(sdzdim), eb(sdzdim), ec(sdzdim)
  double precision :: eal(sdzdim), ebe(sdzdim), ega(sdzdim), eep(nxy, sdzdim)
  double precision :: ed(nxy, sdzdim),  ex(nxy, sdzdim)
  double precision :: flgam(sdzdim+1),  fldel(sdzdim+1)
  double precision :: fltmp

  integer :: iunit = 1
  double precision ::  aunit

  ! Advection choice
  logical :: sadu = .false. !! Upwind ?
  logical :: sadq = .true.  !! Quickest ?
  logical :: sadqf = .true. !!   + Flux limiter ?
  ! Diffusion choice
  logical :: sdcn = .false. !! Crank nicolson ?
  logical :: sdii = .true.  !! Implicit ?
  namelist /nmsdso/ sadu, sadq, sadqf, sdcn, sdii

end module globals_digda_step_solid

! =====================================================================

module module_digda_step_solid
  use globals_digda_step_solid
  implicit none
contains

  subroutine set_digda_step_solid

    integer :: ij, k

    write(nfstdo, *) '*** Sediment : set_digda_step_solid ***'
    rewind(nfparm)
    read(nfparm, nmsdso, end=99)
99 continue
    write(nfstdo, nmsdso)

    if (iunit.eq.1) then
      aunit = 1.0d0 / 8.64d+4 / 365.0d0
    else if (iunit.eq.2) then
      aunit = 1.0d0 / 8.64d+4 / 30.0d0
    else if (iunit.eq.3) then
      aunit = 1.0d0 / 8.64d+4
    else if (iunit.eq.4) then
      aunit = 1.0d0 / 3600.0d0
    else if (iunit.eq.5) then
      aunit = 1.0d0 / 60.0d0
    else if (iunit.eq.6) then
      aunit = 1.0d0
    else
      aunit = 0.0d0
    end if

  end subroutine set_digda_step_solid

! =====================================================================

      subroutine calc_step_solid( ff, rst, w, phi, phii, dbi, s, & !! IN
                                & sn, bfz                        ) !! OUT

      double precision, intent(in)  :: ff(nxy)
      double precision, intent(in)  :: rst(nxy,sdzdim)
      double precision, intent(in)  :: w(nxy,sdzdim+1)
      double precision, intent(in)  :: phi(sdzdim)
      double precision, intent(in)  :: phii(sdzdim+1)
      double precision, intent(in)  :: dbi(sdzdim+1)
      double precision, intent(in)  :: s(nxy,sdzdim)
      double precision, intent(out) :: sn(nxy,sdzdim)
      double precision, intent(out) :: bfz(nxy)

      integer :: ij, k

      do ij = ijtstr, ijtend
         flx(ij, skst) = ff(ij)
      end do

!     -- Advection : Upwind --

      if (sadu) then

      do k = skst+1, sken
      do ij = ijtstr, ijtend
         if ( w(ij, k) >= 0.0d0 ) then
            si(ij, k) = s(ij, k-1)
         else if ( w(ij, k) < 0.0d0 ) then
            si(ij, k) = s(ij, k)
         end if
         flx(ij, k)   = ( 1.0d0 - phii(k) ) * w(ij, k) * si(ij, k)
      end do
      end do

      do ij = ijtstr, ijtend
      if ( w(ij, sken+1) >= 0.0d0 ) then
         si(ij, sken+1) = s(ij, sken)
      else if ( w(ij, sken+1) < 0.0d0 ) then
         si(ij, sken+1) = 0.0d0
      end if
         flx(ij, sken+1) = ( 1.0d0 - phii(sken+1) ) &
     &                   * w(ij, sken+1) * si(ij, sken+1)
      end do

      end if

!     -- Advection : Quickest --

      if (sadq) then

      do k = skst+2, sken-1
      do ij = ijtstr, ijtend
         if ( w(ij, k) >= 0.0d0 ) then
            c2(ij, k) = ( 1.0d0 / ( dlzi(k-1) + dlzi(k) ) ) &
     &         * ( ( ( s(ij, k-2) - s(ij, k-1) ) / dlzi(k-1) ) &
     &           - ( ( s(ij, k-1) - s(ij, k)   ) / dlzi(k)   ) )
         else if ( w(ij, k) < 0.0d0 ) then
            c2(ij, k) = ( 1.0d0 / ( dlzi(k) + dlzi(k+1) ) ) &
     &         * ( ( ( s(ij, k-1) - s(ij, k) ) / dlzi(k)   ) &
     &           - ( ( s(ij, k) - s(ij, k+1) ) / dlzi(k+1) ) )
         end if
            c1(ij, k) = ( ( s(ij, k) - s(ij, k-1) ) / dlzi(k) ) &
     &         + ( ( ( dlz(k-1) - dlz(k) ) / 2.0d0 ) * c2(ij, k) )
            c0(ij, k) = ( ( dlz(k) * s(ij, k-1) + dlz(k-1) * s(ij, k) ) &
     &                / ( 2.0d0 * dlzi(k) ) ) &
     &         - ( ( ( dlz(k-1) * dlz(k) ) / 4.0d0 ) * c2(ij, k) )
      end do
      end do

      do ij = ijtstr, ijtend
      if ( w(ij, skst+1) >= 0.0d0 ) then
         c2(ij, skst+1) = ( 1.0d0 / ( dlzi(skst) + dlzi(skst+1) ) ) &
     &      * ( ( ( 0.0d0       - s(ij, skst)   ) / dlzi(skst) ) &
     &        - ( ( s(ij, skst) - s(ij, skst+1) ) / dlzi(skst+1) ) )
      else if ( w(ij, skst+1) < 0.0d0 ) then
         c2(ij, skst+1) = ( 1.0d0 / ( dlzi(skst+1) + dlzi(skst+2) ) ) &
     &      * ( ( ( s(ij, skst)   - s(ij, skst+1) ) / dlzi(skst+1) ) &
     &        - ( ( s(ij, skst+1) - s(ij, skst+2) ) / dlzi(skst+2) ) )
      end if
         c1(ij, skst+1) = ( ( s(ij, skst+1) - s(ij, skst) ) &
     &                      / dlzi(skst+1) ) &
     &      + ( ( ( dlz(skst) - dlz(skst+1) ) / 2.0d0 ) &
     &          * c2(ij, skst+1) )
         c0(ij, skst+1) = ( ( dlz(skst+1) * s(ij, skst) &
     &                      + dlz(skst)   * s(ij, skst+1) ) &
     &                  / ( 2.0d0 * dlzi(skst+1) ) ) &
     &      - ( ( ( dlz(skst) * dlz(skst+1) ) / 4.0d0 ) &
     &          * c2(ij, skst+1) )
      end do

      do ij = ijtstr, ijtend
      if ( w(ij, sken) >= 0.0d0 ) then
         c2(ij, sken) = ( 1.0d0 / ( dlzi(sken-1) + dlzi(sken) ) ) &
     &      * ( ( ( s(ij, sken-2) - s(ij, sken-1) ) / dlzi(sken-1) ) &
     &        - ( ( s(ij, sken-1) - s(ij, sken)   ) / dlzi(sken)   ) )
      else if ( w(ij, sken) < 0.0d0 ) then
         c2(ij, sken) = ( 1.0d0 / ( dlzi(sken) + dlzi(sken+1) ) ) &
     &      * ( ( ( s(ij, sken-1) - s(ij, sken) ) / dlzi(sken)   ) &
     &        - ( ( s(ij, sken)   - 0.0d0       ) / dlzi(sken+1) ) )
      end if
         c1(ij, sken) = ( ( s(ij, sken) - s(ij, sken-1) ) / dlzi(sken) ) &
     &      + ( ( ( dlz(sken-1) - dlz(sken) ) / 2.0d0 ) &
     &          * c2(ij, sken) )
         c0(ij, sken) = ( ( dlz(sken)   * s(ij, sken-1) &
     &                    + dlz(sken-1) * s(ij, sken)   ) &
     &                / ( 2.0d0 * dlzi(sken) ) ) &
     &      - ( ( ( dlz(sken-1) * dlz(sken) ) / 4.0d0 ) &
     &          * c2(ij, sken) )
      end do

      do ij = ijtstr, ijtend
      if ( w(ij, sken+1) >= 0.0d0 ) then
         c2(ij, sken+1) = ( 1.0d0 / ( dlzi(sken) + dlzi(sken+1) ) ) &
     &      * ( ( ( s(ij, sken-1) - s(ij, sken) ) / dlzi(sken)   ) &
     &        - ( ( s(ij, sken)   - s(ij, sken) ) / dlzi(sken+1) ) )
      else if ( w(ij, sken+1) < 0.0d0 ) then
         c2(ij, sken) = ( 1.0d0 / ( dlzi(sken+1) + dlzi(sken+1) ) ) &
     &      * ( ( ( s(ij, sken) - 0.0d0 ) / dlzi(sken+1) ) &
     &        - ( ( 0.0d0       - 0.0d0 ) / dlzi(sken+1) ) )
      end if
      if ( w(ij, sken+1) >= 0.0d0 ) then
         c1(ij, sken+1) = ( ( s(ij, sken) - s(ij, sken) ) &
     &                    / dlzi(sken+1) )
         c0(ij, sken+1) = ( ( dlz(sken) * s(ij, sken) &
     &                      + dlz(sken) * s(ij, sken) ) &
     &                  / ( 2.0d0 * dlzi(sken+1) ) ) &
     &      - ( ( ( dlz(sken) * dlz(sken) ) / 4.0d0 ) &
     &          * c2(ij, sken+1) )
      else if ( w(ij, sken+1) < 0.0d0 ) then
         c1(ij, sken+1) = ( ( 0.0d0 - s(ij, sken) ) &
     &                    / dlzi(sken+1) )
         c0(ij, sken+1) = ( ( dlz(sken) * s(ij, sken) &
     &                      + dlz(sken) * 0.0d0 ) &
     &                  / ( 2.0d0 * dlzi(sken+1) ) ) &
     &      - ( ( ( dlz(sken) * dlz(sken) ) / 4.0d0 ) &
     &          * c2(ij, sken+1) )
      end if
      end do

      do k = skst+1, sken+1
      do ij = ijtstr, ijtend
         si(ij, k) = c0(ij, k) &
     &   - ( c1(ij, k) * ( ( w(ij, k) * tssec ) / 2.0d0 ) ) &
     &   + ( c2(ij, k) * ( ( ( ( w(ij, k) * tssec ) ** 2.0d0 ) / 3.0d0 ) &
     &             - ( ( ( dlzi(k) ** 2.0d0 ) / 12.0d0 ) ) ) )
      end do
      end do

!     -- Advection : Quickest Flux limiter --

      if (sadqf) then

      do ij = ijtstr, ijtend
      if ( w(ij, skst+1) >= 0.0d0 ) then
      fldel(skst+1) = s(ij, skst+1) - 0.0d0
      flgam(skst+1) = ( ( s(ij, skst+1) - s(ij, skst) ) &
     &                  / dlzi(skst+1) ) &
     &              - ( ( s(ij, skst)   - 0.0d0       ) &
     &                  / dlzi(skst)   )
      if ( abs( flgam(skst+1)*dlz(skst) ) >= abs( fldel(skst+1) ) ) then
         si(ij, skst+1) = s(ij, skst)
      else if ( fldel(skst+1) >= 0 ) then
         fltmp = 0.0d0 &
     &         + ( s(ij, skst) - 0.0d0 ) / w(ij, skst+1) * dlzi(skst+1)
         si(ij, skst+1) = &
     &      min( si(ij, skst+1), min( fltmp, s(ij, skst+1) ) )
         si(ij, skst+1) = max( si(ij, skst+1), s(ij, skst) )
      else if ( fldel(skst+1) <  0 ) then
         fltmp = 0.0d0 &
     &         + ( s(ij, skst) - 0.0d0 ) / w(ij, skst+1) * dlzi(skst+1)
         si(ij, skst+1) = &
     &      max( si(ij, skst+1), max( fltmp, s(ij, skst+1) ) )
         si(ij, skst+1) = min( si(ij, skst+1), s(ij, skst) )
      end if
      else if ( w(ij, skst+1) < 0.0d0 ) then
      fldel(skst+1) = s(ij, skst) - s(ij, skst+2)
      flgam(skst+1) = ( ( s(ij, skst)   - s(ij, skst+1) ) &
     &                  / dlzi(skst+1) ) &
     &              - ( ( s(ij, skst+1) - s(ij, skst+2) ) &
     &                  / dlzi(skst+2) )
      if ( abs( flgam(skst+1)*dlz(skst+1) ) &
     &  >= abs( fldel(skst+1) ) ) then
         si(ij, skst+1) = s(ij, skst+1)
      else if ( fldel(skst+1) >= 0 ) then
         fltmp = s(ij, skst+2) &
     &- ( s(ij, skst+1) - s(ij, skst+2) ) / w(ij, skst+1) * dlzi(skst+1)
         si(ij, skst+1) = &
     &      min( si(ij, skst+1), min( fltmp, s(ij, skst) ) )
         si(ij, skst+1) = max( si(ij, skst+1), s(ij, skst+1) )
      else if ( fldel(skst+1) <  0 ) then
         fltmp = s(ij, skst+2) &
     &- ( s(ij, skst+1) - s(ij, skst+2) ) / w(ij, skst+1) * dlzi(skst+1)
         si(ij, skst+1) = &
     &      max( si(ij, skst+1), max( fltmp, s(ij, skst) ) )
         si(ij, skst+1) = min( si(ij, skst+1), s(ij, skst+1) )
      end if
      end if
      end do

      do k = skst+2, sken-1
      do ij = ijtstr, ijtend
      if ( w(ij, k) >= 0.0d0 ) then
         fldel(k) = s(ij, k) - s(ij, k-2)
         flgam(k) = ( ( s(ij, k)   - s(ij, k-1) ) / dlzi(k)   ) &
     &            - ( ( s(ij, k-1) - s(ij, k-2) ) / dlzi(k-1) )
         if ( abs( flgam(k)*dlz(k-1) ) >= abs( fldel(k) ) ) then
            si(ij, k) = s(ij, k-1)
         else if ( fldel(k) >= 0 ) then
            fltmp = s(ij, k-2) &
     &            + ( s(ij, k-1) - s(ij, k-2) ) / w(ij, k) * dlzi(k)
            si(ij, k) = &
     &         min( si(ij, k), min( fltmp, s(ij, k) ) )
            si(ij, k) = max( si(ij, k), s(ij, k-1) )
         else if ( fldel(k) <  0 ) then
            fltmp = s(ij, k-2) &
     &            + ( s(ij, k-1) - s(ij, k-2) ) / w(ij, k) * dlzi(k)
            si(ij, k) = &
     &         max( si(ij, k), max( fltmp, s(ij, k) ) )
            si(ij, k) = min( si(ij, k), s(ij, k-1) )
         end if
      else if ( w(ij, k) < 0.0d0 ) then
         fldel(k) = s(ij, k-1) - s(ij, k+1)
         flgam(k) = ( ( s(ij, k-1) - s(ij, k)   ) / dlzi(k)   ) &
     &            - ( ( s(ij, k)   - s(ij, k+1) ) / dlzi(k+1) )
         if ( abs( flgam(k)*dlz(k) ) >= abs( fldel(k) ) ) then
            si(ij, k) = s(ij, k)
         else if ( fldel(k) >= 0 ) then
            fltmp = s(ij, k+1) &
     &            - ( s(ij, k) - s(ij, k+1) ) / w(ij, k) * dlzi(k)
            si(ij, k) = &
     &         min( si(ij, k), min( fltmp, s(ij, k-1) ) )
            si(ij, k) = max( si(ij, k), s(ij, k) )
         else if ( fldel(k) <  0 ) then
            fltmp = s(ij, k+1) &
     &            - ( s(ij, k) - s(ij, k+1) ) / w(ij, k) * dlzi(k)
            si(ij, k) = &
     &         max( si(ij, k), max( fltmp, s(ij, k-1) ) )
            si(ij, k) = min( si(ij, k), s(ij, k) )
         end if
      end if
      end do
      end do

      do ij = ijtstr, ijtend
      if ( w(ij, sken) >= 0.0d0 ) then
      fldel(sken) = s(ij, sken) - s(ij, sken-2)
      flgam(sken) = ( ( s(ij, sken)   - s(ij, sken-1) ) / dlzi(sken)   ) &
     &            - ( ( s(ij, sken-1) - s(ij, sken-2) ) / dlzi(sken-1) )
      if ( abs( flgam(sken)*dlz(sken-1) ) >= abs( fldel(sken) ) ) then
         si(ij, sken) = s(ij, sken-1)
      else if ( fldel(sken) >= 0 ) then
         fltmp = s(ij, sken-2) &
     &+ ( s(ij, sken-1) - s(ij, sken-2) ) / w(ij, sken) * dlzi(sken)
         si(ij, sken) = &
     &      min( si(ij, sken), min( fltmp, s(ij, sken) ) )
         si(ij, sken) = max( si(ij, sken), s(ij, sken-1) )
      else if ( fldel(sken) <  0 ) then
         fltmp = s(ij, sken-2) &
     &+ ( s(ij, sken-1) - s(ij, sken-2) ) / w(ij, sken) * dlzi(sken)
         si(ij, sken) = &
     &      max( si(ij, sken), max( fltmp, s(ij, sken) ) )
         si(ij, sken) = min( si(ij, sken), s(ij, sken-1) )
      end if
      else if ( w(ij, sken) < 0.0d0 ) then
      fldel(sken) = s(ij, sken-1) - 0.0d0
      flgam(sken) = ( ( s(ij, sken-1) - s(ij, sken) ) / dlzi(sken)   ) &
     &            - ( ( s(ij, sken)   - 0.0d0       ) / dlzi(sken+1) )
      if ( abs( flgam(sken)*dlz(sken) ) >= abs( fldel(sken) ) ) then
         si(ij, sken) = s(ij, sken)
      else if ( fldel(sken) >= 0 ) then
         fltmp = 0.0d0 &
     &- ( s(ij, sken) - 0.0d0 ) / w(ij, sken) * dlzi(sken)
         si(ij, sken) = &
     &      min( si(ij, sken), min( fltmp, s(ij, sken-1) ) )
         si(ij, sken) = max( si(ij, sken), s(ij, sken) )
      else if ( fldel(sken) <  0 ) then
         fltmp = 0.0d0 &
     &- ( s(ij, sken) - 0.0d0 ) / w(ij, sken) * dlzi(sken)
         si(ij, sken) = &
     &      max( si(ij, sken), max( fltmp, s(ij, sken-1) ) )
         si(ij, sken) = min( si(ij, sken), s(ij, sken) )
      end if
      end if
      end do

      do ij = ijtstr, ijtend
      if ( w(ij, sken+1) >= 0.0d0 ) then
      fldel(sken+1) = s(ij, sken) - s(ij, sken-1)
      flgam(sken+1) = ( ( s(ij, sken) - s(ij, sken)   ) / dlzi(sken+1) ) &
     &              - ( ( s(ij, sken) - s(ij, sken-1) ) / dlzi(sken)   )
      if ( abs( flgam(sken+1)*dlz(sken) ) >= abs( fldel(sken+1) ) ) then
         si(ij, sken+1) = s(ij, sken)
      else if ( fldel(sken+1) >= 0 ) then
         fltmp = s(ij, sken-1) &
     &+ ( s(ij, sken) - s(ij, sken-1) ) / w(ij, sken+1) * dlzi(sken+1)
         si(ij, sken+1) = &
     &      min( si(ij, sken+1), min( fltmp, s(ij, sken) ) )
         si(ij, sken+1) = max( si(ij, sken+1), s(ij, sken) )
      else if ( fldel(sken+1) <  0 ) then
         fltmp = s(ij, sken-1) &
     &+ ( s(ij, sken) - s(ij, sken-1) ) / w(ij, sken+1) * dlzi(sken+1)
         si(ij, sken+1) = &
     &      max( si(ij, sken+1), max( fltmp, s(ij, sken) ) )
         si(ij, sken+1) = min( si(ij, sken+1), s(ij, sken) )
      end if
      else if ( w(ij, sken+1) < 0.0d0 ) then
      fldel(sken+1) = s(ij, sken) - 0.0d0
      flgam(sken+1) = ( ( s(ij, sken) - 0.0d0 ) / dlzi(sken+1) ) &
     &              - ( ( 0.0d0       - 0.0d0 ) / dlzi(sken+1) )
      if ( abs( flgam(sken+1)*dlz(sken) ) >= abs( fldel(sken+1) ) ) then
         si(ij, sken+1) = 0.0d0
      else if ( fldel(sken+1) >= 0 ) then
         fltmp = 0.0d0 &
     &- ( 0.0d0 - 0.0d0 ) / w(ij, sken+1) * dlzi(sken+1)
         si(ij, sken+1) = &
     &      min( si(ij, sken+1), min( fltmp, s(ij, sken) ) )
         si(ij, sken+1) = max( si(ij, sken+1), 0.0d0 )
      else if ( fldel(sken+1) <  0 ) then
         fltmp = 0.0d0 &
     &- ( 0.0d0 - 0.0d0 ) / w(ij, sken+1) * dlzi(sken+1)
         si(ij, sken+1) = &
     &      max( si(ij, sken+1), max( fltmp, s(ij, sken) ) )
         si(ij, sken+1) = min( si(ij, sken+1), 0.0d0 )
      end if
      end if
      end do

      end if

      do k = skst+1, sken+1
      do ij = ijtstr, ijtend
         flx(ij, k) = ( 1.0d0 - phii(k) ) * w(ij, k) * si(ij, k)
      end do
      end do

      end if

!     -- Diffusion : Crank nicolson

      if (sdcn) then

      do k = skst, sken
         eal(k) = ( ( 1.0d0 - phi(k) ) / tssec )
         ebe(k) = ( ( ( 1.0d0 - phii(k) )   * dbi(k) ) &
     &            / ( 2.0d0 * dlz(k) * dlzi(k)   ) )
         ega(k) = ( ( ( 1.0d0 - phii(k+1) ) * dbi(k+1) ) &
     &            / ( 2.0d0 * dlz(k) * dlzi(k+1) ) )
      end do

      do k = skst+1, sken-1
         ea(k)     = - ebe(k)
         eb(k)     = eal(k) + ebe(k) + ega(k)
         ec(k)     = - ega(k)
      end do

         ea(skst)  = - ebe(skst)
         eb(skst)  = eal(skst) + ega(skst)
         ec(skst)  = - ega(skst)
         ea(sken)  = - ebe(sken)
         eb(sken)  = eal(sken) + ebe(sken)
         ec(sken)  = - ega(sken)

      do ij = ijtstr, ijtend
         ed(ij, skst) = &
     &      + ( eal(skst) - ega(skst) ) * s(ij, skst) &
     &      + ega(skst) * s(ij, skst+1) &
     &      + ( ( flx(ij, skst) - flx(ij, skst+1) ) / dlz(skst) ) &
     &      - rst(ij, skst)
      end do

      do k = skst+1, sken-1
      do ij = ijtstr, ijtend
         ed(ij, k) = &
     &      + ebe(k) * s(ij, k-1) &
     &      + ( eal(k) - ebe(k) - ega(k) ) * s(ij, k) &
     &      + ega(k) * s(ij, k+1) &
     &      + ( ( flx(ij, k) - flx(ij, k+1) ) / dlz(k) ) &
     &      - rst(ij, k)
      end do
      end do

      do ij = ijtstr, ijtend
         ed(ij, sken) = &
     &      + ebe(sken) * s(ij, sken-1) &
     &      + ( eal(sken) - ebe(sken) ) * s(ij, sken) &
     &      + ( ( flx(ij, sken) - flx(ij, sken+1) ) / dlz(sken) ) &
     &      - rst(ij, sken)
      end do

      call calc_tridiagonal_solid( ea, eb, ec, ed, ex )

      end if

!     -- Diffusion : Implicit

      if (sdii) then

      do k = skst, sken
         eal(k) = ( ( 1.0d0 - phi(k) ) / tssec )
         ebe(k) = ( ( ( 1.0d0 - phii(k) )   * dbi(k) ) &
     &            / ( dlz(k) * dlzi(k)   ) )
         ega(k) = ( ( ( 1.0d0 - phii(k+1) ) * dbi(k+1) ) &
     &            / ( dlz(k) * dlzi(k+1) ) )
      end do

      do k = skst+1, sken-1
         ea(k)     = - ebe(k)
         eb(k)     = eal(k) + ebe(k) + ega(k)
         ec(k)     = - ega(k)
      end do

         ea(skst)  = - ebe(skst)
         eb(skst)  = eal(skst) + ega(skst)
         ec(skst)  = - ega(skst)
         ea(sken)  = - ebe(sken)
         eb(sken)  = eal(sken) + ebe(sken)
         ec(sken)  = - ega(sken)

      do k = skst, sken
      do ij = ijtstr, ijtend
         ed(ij, k) = &
     &      + eal(k) * s(ij, k) &
     &      + ( ( flx(ij, k) - flx(ij, k+1) ) / dlz(k) ) &
     &      - rst(ij, k)
      end do
      end do

      call calc_tridiagonal_solid( ea, eb, ec, ed, ex )

      end if

!     -- t = t + dt --

      do k = skst, sken
      do ij = ijtstr, ijtend
         sn(ij, k) = max( ex(ij, k), 0.0d0 )
         ex(ij, k) = 0.0d0
      end do
      end do

      do ij = ijtstr, ijtend
         bfz(ij) = flx(ij, sken+1)
      end do

      end subroutine calc_step_solid

! =====================================================================

      subroutine calc_step_solid_clay( ff, rst, w, phi, phii, dbi, S, rhocla, & !! IN
                                     & sn, bfz                                ) !! OUT

      double precision, intent(in)  :: ff(nxy)
      double precision, intent(in)  :: rst(nxy,sdzdim)
      double precision, intent(in)  :: w(nxy,sdzdim+1)
      double precision, intent(in)  :: phi(sdzdim)
      double precision, intent(in)  :: phii(sdzdim+1)
      double precision, intent(in)  :: dbi(sdzdim+1)
      double precision, intent(in)  :: s(nxy,sdzdim)
      double precision, intent(in)  :: rhocla
      double precision, intent(out) :: sn(nxy,sdzdim)
      double precision, intent(out) :: bfz(nxy)

      integer :: ij, k

      do ij = ijtstr, ijtend
         flx(ij, skst) = ff(ij)
      end do

!     -- Advection : Upwind --

      if (sadu) then

      do k = skst+1, sken
      do ij = ijtstr, ijtend
         if ( w(ij, k) >= 0.0d0 ) then
            si(ij, k) = s(ij, k-1)
         else if ( w(ij, k) < 0.0d0 ) then
            si(ij, k) = s(ij, k)
         end if
         flx(ij, k)   = ( 1.0d0 - phii(k) ) * w(ij, k) * si(ij, k)
      end do
      end do

      do ij = ijtstr, ijtend
      if ( w(ij, sken+1) >= 0.0d0 ) then
         si(ij, sken+1) = s(ij, sken)
      else if ( w(ij, sken+1) < 0.0d0 ) then
         si(ij, sken+1) = 1.0d+6 * rhocla
      end if
         flx(ij, sken+1) = ( 1.0d0 - phii(sken+1) ) &
     &                   * w(ij, sken+1) * si(ij, sken+1)
      end do

      end if

!     -- Advection : Quickest --

      if (sadq) then

      do k = skst+2, sken-1
      do ij = ijtstr, ijtend
         if ( w(ij, k) >= 0.0d0 ) then
            c2(ij, k) = ( 1.0d0 / ( dlzi(k-1) + dlzi(k) ) ) &
     &         * ( ( ( s(ij, k-2) - s(ij, k-1) ) / dlzi(k-1) ) &
     &           - ( ( s(ij, k-1) - s(ij, k)   ) / dlzi(k)   ) )
         else if ( w(ij, k) < 0.0d0 ) then
            c2(ij, k) = ( 1.0d0 / ( dlzi(k) + dlzi(k+1) ) ) &
     &         * ( ( ( s(ij, k-1) - s(ij, k) ) / dlzi(k)   ) &
     &           - ( ( s(ij, k) - s(ij, k+1) ) / dlzi(k+1) ) )
         end if
            c1(ij, k) = ( ( s(ij, k) - s(ij, k-1) ) / dlzi(k) ) &
     &         + ( ( ( dlz(k-1) - dlz(k) ) / 2.0d0 ) * c2(ij, k) )
            c0(ij, k) = ( ( dlz(k) * s(ij, k-1) + dlz(k-1) * s(ij, k) ) &
     &                / ( 2.0d0 * dlzi(k) ) ) &
     &         - ( ( ( dlz(k-1) * dlz(k) ) / 4.0d0 ) * c2(ij, k) )
      end do
      end do

      do ij = ijtstr, ijtend
      if ( w(ij, skst+1) >= 0.0d0 ) then
         c2(ij, skst+1) = ( 1.0d0 / ( dlzi(skst) + dlzi(skst+1) ) ) &
     &      * ( ( ( 0.0d0       - s(ij, skst)   ) / dlzi(skst) ) &
     &        - ( ( s(ij, skst) - s(ij, skst+1) ) / dlzi(skst+1) ) )
      else if ( w(ij, skst+1) < 0.0d0 ) then
         c2(ij, skst+1) = ( 1.0d0 / ( dlzi(skst+1) + dlzi(skst+2) ) ) &
     &      * ( ( ( s(ij, skst)   - s(ij, skst+1) ) / dlzi(skst+1) ) &
     &        - ( ( s(ij, skst+1) - s(ij, skst+2) ) / dlzi(skst+2) ) )
      end if
         c1(ij, skst+1) = ( ( s(ij, skst+1) - s(ij, skst) ) &
     &                      / dlzi(skst+1) ) &
     &      + ( ( ( dlz(skst) - dlz(skst+1) ) / 2.0d0 ) &
     &          * c2(ij, skst+1) )
         c0(ij, skst+1) = ( ( dlz(skst+1) * s(ij, skst) &
     &                      + dlz(skst)   * s(ij, skst+1) ) &
     &                  / ( 2.0d0 * dlzi(skst+1) ) ) &
     &      - ( ( ( dlz(skst) * dlz(skst+1) ) / 4.0d0 ) &
     &          * c2(ij, skst+1) )
      end do

      do ij = ijtstr, ijtend
      if ( w(ij, sken) >= 0.0d0 ) then
         c2(ij, sken) = ( 1.0d0 / ( dlzi(sken-1) + dlzi(sken) ) ) &
     &      * ( ( ( s(ij, sken-2) - s(ij, sken-1) ) / dlzi(sken-1) ) &
     &        - ( ( s(ij, sken-1) - s(ij, sken)   ) / dlzi(sken)   ) )
      else if ( w(ij, sken) < 0.0d0 ) then
         c2(ij, sken) = ( 1.0d0 / ( dlzi(sken) + dlzi(sken+1) ) ) &
     &      * ( ( ( s(ij, sken-1) - s(ij, sken) ) / dlzi(sken)   ) &
     &        - ( ( s(ij, sken)   - 1.0d+6 * rhocla ) / dlzi(sken+1) ) )
      end if
         c1(ij, sken) = ( ( s(ij, sken) - s(ij, sken-1) ) / dlzi(sken) ) &
     &      + ( ( ( dlz(sken-1) - dlz(sken) ) / 2.0d0 ) &
     &          * c2(ij, sken) )
         c0(ij, sken) = ( ( dlz(sken)   * s(ij, sken-1) &
     &                    + dlz(sken-1) * s(ij, sken)   ) &
     &                / ( 2.0d0 * dlzi(sken) ) ) &
     &      - ( ( ( dlz(sken-1) * dlz(sken) ) / 4.0d0 ) &
     &          * c2(ij, sken) )
      end do

      do ij = ijtstr, ijtend
      if ( w(ij, sken+1) >= 0.0d0 ) then
         c2(ij, sken+1) = ( 1.0d0 / ( dlzi(sken) + dlzi(sken+1) ) ) &
     &      * ( ( ( s(ij, sken-1) - s(ij, sken) ) / dlzi(sken)   ) &
     &        - ( ( s(ij, sken)   - s(ij, sken) ) / dlzi(sken+1) ) )
      else if ( w(ij, sken+1) < 0.0d0 ) then
         c2(ij, sken) = ( 1.0d0 / ( dlzi(sken+1) + dlzi(sken+1) ) ) &
     &      * ( ( ( s(ij, sken) - 1.0d+6 * rhocla ) / dlzi(sken+1) ) &
     &    - ( ( 1.0d+6 * rhocla - 1.0d+6 * rhocla ) / dlzi(sken+1) ) )
      end if
      if ( w(ij, sken+1) >= 0.0d0 ) then
         c1(ij, sken+1) = ( ( s(ij, sken) - s(ij, sken) ) &
     &                    / dlzi(sken+1) )
         c0(ij, sken+1) = ( ( dlz(sken) * s(ij, sken) &
     &                      + dlz(sken) * s(ij, sken) ) &
     &                  / ( 2.0d0 * dlzi(sken+1) ) ) &
     &      - ( ( ( dlz(sken) * dlz(sken) ) / 4.0d0 ) &
     &          * c2(ij, sken+1) )
      else if ( w(ij, sken+1) < 0.0d0 ) then
         c1(ij, sken+1) = ( ( 1.0d+6 * rhocla - s(ij, sken) ) &
     &                    / dlzi(sken+1) )
         c0(ij, sken+1) = ( ( dlz(sken) * s(ij, sken) &
     &                      + dlz(sken) * 1.0d+6 * rhocla ) &
     &                  / ( 2.0d0 * dlzi(sken+1) ) ) &
     &      - ( ( ( dlz(sken) * dlz(sken) ) / 4.0d0 ) &
     &          * c2(ij, sken+1) )
      end if
      end do

      do k = skst+1, sken+1
      do ij = ijtstr, ijtend
         si(ij, k) = c0(ij, k) &
     &   - ( c1(ij, k) * ( ( w(ij, k) * tssec ) / 2.0d0 ) ) &
     &   + ( c2(ij, k) * ( ( ( ( w(ij, k) * tssec ) ** 2.0d0 ) / 3.0d0 ) &
     &             - ( ( ( dlzi(k) ** 2.0d0 ) / 12.0d0 ) ) ) )
      end do
      end do

!     -- Advection : Quickest Flux limiter --

      if (sadqf) then

      do ij = ijtstr, ijtend
      if ( w(ij, skst+1) >= 0.0d0 ) then
      fldel(skst+1) = s(ij, skst+1) - 0.0d0
      flgam(skst+1) = ( ( s(ij, skst+1) - s(ij, skst) ) &
     &                  / dlzi(skst+1) ) &
     &              - ( ( s(ij, skst)   - 0.0d0       ) &
     &                  / dlzi(skst)   )
      if ( abs( flgam(skst+1)*dlz(skst) ) >= abs( fldel(skst+1) ) ) then
         si(ij, skst+1) = s(ij, skst)
      else if ( fldel(skst+1) >= 0 ) then
         fltmp = 0.0d0 &
     &         + ( s(ij, skst) - 0.0d0 ) / w(ij, skst+1) * dlzi(skst+1)
         si(ij, skst+1) = &
     &      min( si(ij, skst+1), min( fltmp, s(ij, skst+1) ) )
         si(ij, skst+1) = max( si(ij, skst+1), s(ij, skst) )
      else if ( fldel(skst+1) <  0 ) then
         fltmp = 0.0d0 &
     &         + ( s(ij, skst) - 0.0d0 ) / w(ij, skst+1) * dlzi(skst+1)
         si(ij, skst+1) = &
     &      max( si(ij, skst+1), max( fltmp, s(ij, skst+1) ) )
         si(ij, skst+1) = min( si(ij, skst+1), s(ij, skst) )
      end if
      else if ( w(ij, skst+1) < 0.0d0 ) then
      fldel(skst+1) = s(ij, skst) - s(ij, skst+2)
      flgam(skst+1) = ( ( s(ij, skst)   - s(ij, skst+1) ) &
     &                  / dlzi(skst+1) ) &
     &              - ( ( s(ij, skst+1) - s(ij, skst+2) ) &
     &                  / dlzi(skst+2) )
      if ( abs( flgam(skst+1)*dlz(skst+1) ) &
     &  >= abs( fldel(skst+1) ) ) then
         si(ij, skst+1) = s(ij, skst+1)
      else if ( fldel(skst+1) >= 0 ) then
         fltmp = s(ij, skst+2) &
     &- ( s(ij, skst+1) - s(ij, skst+2) ) / w(ij, skst+1) * dlzi(skst+1)
         si(ij, skst+1) = &
     &      min( si(ij, skst+1), min( fltmp, s(ij, skst) ) )
         si(ij, skst+1) = max( si(ij, skst+1), s(ij, skst+1) )
      else if ( fldel(skst+1) <  0 ) then
         fltmp = s(ij, skst+2) &
     &- ( s(ij, skst+1) - s(ij, skst+2) ) / w(ij, skst+1) * dlzi(skst+1)
         si(ij, skst+1) = &
     &      max( si(ij, skst+1), max( fltmp, s(ij, skst) ) )
         si(ij, skst+1) = min( si(ij, skst+1), s(ij, skst+1) )
      end if
      end if
      end do

      do k = skst+2, sken-1
      do ij = ijtstr, ijtend
      if ( w(ij, k) >= 0.0d0 ) then
         fldel(k) = s(ij, k) - s(ij, k-2)
         flgam(k) = ( ( s(ij, k)   - s(ij, k-1) ) / dlzi(k)   ) &
     &            - ( ( s(ij, k-1) - s(ij, k-2) ) / dlzi(k-1) )
         if ( abs( flgam(k)*dlz(k-1) ) >= abs( fldel(k) ) ) then
            si(ij, k) = s(ij, k-1)
         else if ( fldel(k) >= 0 ) then
            fltmp = s(ij, k-2) &
     &            + ( s(ij, k-1) - s(ij, k-2) ) / w(ij, k) * dlzi(k)
            si(ij, k) = &
     &         min( si(ij, k), min( fltmp, s(ij, k) ) )
            si(ij, k) = max( si(ij, k), s(ij, k-1) )
         else if ( fldel(k) <  0 ) then
            fltmp = s(ij, k-2) &
     &            + ( s(ij, k-1) - s(ij, k-2) ) / w(ij, k) * dlzi(k)
            si(ij, k) = &
     &         max( si(ij, k), max( fltmp, s(ij, k) ) )
            si(ij, k) = min( si(ij, k), s(ij, k-1) )
         end if
      else if ( w(ij, k) < 0.0d0 ) then
         fldel(k) = s(ij, k-1) - s(ij, k+1)
         flgam(k) = ( ( s(ij, k-1) - s(ij, k)   ) / dlzi(k)   ) &
     &            - ( ( s(ij, k)   - s(ij, k+1) ) / dlzi(k+1) )
         if ( abs( flgam(k)*dlz(k) ) >= abs( fldel(k) ) ) then
            si(ij, k) = s(ij, k)
         else if ( fldel(k) >= 0 ) then
            fltmp = s(ij, k+1) &
     &            - ( s(ij, k) - s(ij, k+1) ) / w(ij, k) * dlzi(k)
            si(ij, k) = &
     &         min( si(ij, k), min( fltmp, s(ij, k-1) ) )
            si(ij, k) = max( si(ij, k), s(ij, k) )
         else if ( fldel(k) <  0 ) then
            fltmp = s(ij, k+1) &
     &            - ( s(ij, k) - s(ij, k+1) ) / w(ij, k) * dlzi(k)
            si(ij, k) = &
     &         max( si(ij, k), max( fltmp, s(ij, k-1) ) )
            si(ij, k) = min( si(ij, k), s(ij, k) )
         end if
      end if
      end do
      end do

      do ij = ijtstr, ijtend
      if ( w(ij, sken) >= 0.0d0 ) then
      fldel(sken) = s(ij, sken) - s(ij, sken-2)
      flgam(sken) = ( ( s(ij, sken)   - s(ij, sken-1) ) / dlzi(sken)   ) &
     &            - ( ( s(ij, sken-1) - s(ij, sken-2) ) / dlzi(sken-1) )
      if ( abs( flgam(sken)*dlz(sken-1) ) >= abs( fldel(sken) ) ) then
         si(ij, sken) = s(ij, sken-1)
      else if ( fldel(sken) >= 0 ) then
         fltmp = s(ij, sken-2) &
     &+ ( s(ij, sken-1) - s(ij, sken-2) ) / w(ij, sken) * dlzi(sken)
         si(ij, sken) = &
     &      min( si(ij, sken), min( fltmp, s(ij, sken) ) )
         si(ij, sken) = max( si(ij, sken), s(ij, sken-1) )
      else if ( fldel(sken) <  0 ) then
         fltmp = s(ij, sken-2) &
     &+ ( s(ij, sken-1) - s(ij, sken-2) ) / w(ij, sken) * dlzi(sken)
         si(ij, sken) = &
     &      max( si(ij, sken), max( fltmp, s(ij, sken) ) )
         si(ij, sken) = min( si(ij, sken), s(ij, sken-1) )
      end if
      else if ( w(ij, sken) < 0.0d0 ) then
      fldel(sken) = s(ij, sken-1) - 1.0d+6 * rhocla
      flgam(sken) = ( ( s(ij, sken-1) - s(ij, sken) ) / dlzi(sken)   ) &
     &            - ( ( s(ij, sken) - 1.0d+6 * rhocla ) / dlzi(sken+1) )
      if ( abs( flgam(sken)*dlz(sken) ) >= abs( fldel(sken) ) ) then
         si(ij, sken) = s(ij, sken)
      else if ( fldel(sken) >= 0 ) then
         fltmp = 1.0d+6 * rhocla &
     &- ( s(ij, sken) - 1.0d+6 * rhocla ) / w(ij, sken) * dlzi(sken)
         si(ij, sken) = &
     &      min( si(ij, sken), min( fltmp, s(ij, sken-1) ) )
         si(ij, sken) = max( si(ij, sken), s(ij, sken) )
      else if ( fldel(sken) <  0 ) then
         fltmp = 1.0d+6 * rhocla &
     &- ( s(ij, sken) - 1.0d+6 * rhocla ) / w(ij, sken) * dlzi(sken)
         si(ij, sken) = &
     &      max( si(ij, sken), max( fltmp, s(ij, sken-1) ) )
         si(ij, sken) = min( si(ij, sken), s(ij, sken) )
      end if
      end if
      end do

      do ij = ijtstr, ijtend
      if ( w(ij, sken+1) >= 0.0d0 ) then
      fldel(sken+1) = s(ij, sken) - s(ij, sken-1)
      flgam(sken+1) = ( ( s(ij, sken) - s(ij, sken)   ) / dlzi(sken+1) ) &
     &              - ( ( s(ij, sken) - s(ij, sken-1) ) / dlzi(sken)   )
      if ( abs( flgam(sken+1)*dlz(sken) ) >= abs( fldel(sken+1) ) ) then
         si(ij, sken+1) = s(ij, sken)
      else if ( fldel(sken+1) >= 0 ) then
         fltmp = s(ij, sken-1) &
     &+ ( s(ij, sken) - s(ij, sken-1) ) / w(ij, sken+1) * dlzi(sken+1)
         si(ij, sken+1) = &
     &      min( si(ij, sken+1), min( fltmp, s(ij, sken) ) )
         si(ij, sken+1) = max( si(ij, sken+1), s(ij, sken) )
      else if ( fldel(sken+1) <  0 ) then
         fltmp = s(ij, sken-1) &
     &+ ( s(ij, sken) - s(ij, sken-1) ) / w(ij, sken+1) * dlzi(sken+1)
         si(ij, sken+1) = &
     &      max( si(ij, sken+1), max( fltmp, s(ij, sken) ) )
         si(ij, sken+1) = min( si(ij, sken+1), s(ij, sken) )
      end if
      else if ( w(ij, sken+1) < 0.0d0 ) then
      fldel(sken+1) = s(ij, sken) - 1.0d+6 * rhocla
      flgam(sken+1) &
     &= ( ( s(ij, sken)     - 1.0d+6 * rhocla ) / dlzi(sken+1) ) &
     &- ( ( 1.0d+6 * rhocla - 1.0d+6 * rhocla ) / dlzi(sken+1) )
      if ( abs( flgam(sken+1)*dlz(sken) ) >= abs( fldel(sken+1) ) ) then
         si(ij, sken+1) = 1.0d+6 * rhocla
      else if ( fldel(sken+1) >= 0 ) then
         fltmp = 1.0d+6 * rhocla &
     &         - ( 1.0d+6 * rhocla - 1.0d+6 * rhocla ) &
     &         / w(ij, sken+1) * dlzi(sken+1)
         si(ij, sken+1) = &
     &      min( si(ij, sken+1), min( fltmp, s(ij, sken) ) )
         si(ij, sken+1) = max( si(ij, sken+1), 1.0d+6 * rhocla )
      else if ( fldel(sken+1) <  0 ) then
         fltmp = 1.0d+6 * rhocla &
     &         - ( 1.0d+6 * rhocla - 1.0d+6 * rhocla ) &
     &         / w(ij, sken+1) * dlzi(sken+1)
         si(ij, sken+1) = &
     &      max( si(ij, sken+1), max( fltmp, s(ij, sken) ) )
         si(ij, sken+1) = min( si(ij, sken+1), 1.0d+6 * rhocla )
      end if
      end if
      end do

      end if

      do k = skst+1, sken+1
      do ij = ijtstr, ijtend
         flx(ij, k) = ( 1.0d0 - phii(k) ) * w(ij, k) * si(ij, k)
      end do
      end do

      end if

!     -- Diffusion : Crank nicolson

      if (sdcn) then

      do k = skst, sken
         eal(k) = ( ( 1.0d0 - phi(k) ) / tssec )
         ebe(k) = ( ( ( 1.0d0 - phii(k) )   * dbi(k) ) &
     &            / ( 2.0d0 * dlz(k) * dlzi(k)   ) )
         ega(k) = ( ( ( 1.0d0 - phii(k+1) ) * dbi(k+1) ) &
     &            / ( 2.0d0 * dlz(k) * dlzi(k+1) ) )
      end do

      do k = skst+1, sken-1
         ea(k)     = - ebe(k)
         eb(k)     = eal(k) + ebe(k) + ega(k)
         ec(k)     = - ega(k)
      end do

         ea(skst)  = - ebe(skst)
         eb(skst)  = eal(skst) + ega(skst)
         ec(skst)  = - ega(skst)
         ea(sken)  = - ebe(sken)
         eb(sken)  = eal(sken) + ebe(sken)
         ec(sken)  = - ega(sken)

      do ij = ijtstr, ijtend
         ed(ij, skst) = &
     &      + ( eal(skst) - ega(skst) ) * s(ij, skst) &
     &      + ega(skst) * s(ij, skst+1) &
     &      + ( ( flx(ij, skst) - flx(ij, skst+1) ) / dlz(skst) ) &
     &      - rst(ij, skst)
      end do

      do k = skst+1, sken-1
      do ij = ijtstr, ijtend
         ed(ij, k) = &
     &      + ebe(k) * s(ij, k-1) &
     &      + ( eal(k) - ebe(k) - ega(k) ) * s(ij, k) &
     &      + ega(k) * s(ij, k+1) &
     &      + ( ( flx(ij, k) - flx(ij, k+1) ) / dlz(k) ) &
     &      - rst(ij, k)
      end do
      end do

      do ij = ijtstr, ijtend
         ed(ij, sken) = &
     &      + ebe(sken) * s(ij, sken-1) &
     &      + ( eal(sken) - ebe(sken) ) * s(ij, sken) &
     &      + ( ( flx(ij, sken) - flx(ij, sken+1) ) / dlz(sken) ) &
     &      - rst(ij, sken)
      end do

      call calc_tridiagonal_solid( ea, eb, ec, ed, ex )

      end if

!     -- Diffusion : Implicit

      if (sdii) then

      do k = skst, sken
         eal(k) = ( ( 1.0d0 - phi(k) ) / tssec )
         ebe(k) = ( ( ( 1.0d0 - phii(k) )   * dbi(k) ) &
     &            / ( dlz(k) * dlzi(k)   ) )
         ega(k) = ( ( ( 1.0d0 - phii(k+1) ) * dbi(k+1) ) &
     &            / ( dlz(k) * dlzi(k+1) ) )
      end do

      do k = skst+1, sken-1
         ea(k)     = - ebe(k)
         eb(k)     = eal(k) + ebe(k) + ega(k)
         ec(k)     = - ega(k)
      end do

         ea(skst)  = - ebe(skst)
         eb(skst)  = eal(skst) + ega(skst)
         ec(skst)  = - ega(skst)
         ea(sken)  = - ebe(sken)
         eb(sken)  = eal(sken) + ebe(sken)
         ec(sken)  = - ega(sken)

      do k = skst, sken
      do ij = ijtstr, ijtend
         ed(ij, k) = &
     &      + eal(k) * s(ij, k) &
     &      + ( ( flx(ij, k) - flx(ij, k+1) ) / dlz(k) ) &
     &      - rst(ij, k)
      end do
      end do

      call calc_tridiagonal_solid( ea, eb, ec, ed, ex )

      end if

!     -- t = t + dt --

      do k = skst, sken
      do ij = ijtstr, ijtend
         sn(ij, k) = max( ex(ij, k), 0.0d0 )
         ex(ij, k) = 0.0d0
      end do
      end do

      do ij = ijtstr, ijtend
         bfz(ij) = flx(ij, sken+1)
      end do

      end subroutine calc_step_solid_clay

! =====================================================================
      subroutine calc_step_solid_age( &
     &      ff, w, phi, phii, dbi, S, &
     &      SN )

      double precision, intent(in)  :: ff(nxy)
      double precision, intent(in)  :: w(nxy,sdzdim+1)
      double precision, intent(in)  :: phi(sdzdim)
      double precision, intent(in)  :: phii(sdzdim+1)
      double precision, intent(in)  :: dbi(sdzdim+1)
      double precision, intent(in)  :: s(nxy,sdzdim)
      double precision, intent(out) :: sn(nxy,sdzdim)

      integer :: ij, k

      do ij = ijtstr, ijtend
         flx(ij, skst) = ff(ij)
      end do

!     -- Advection : Upwind --

      if (sadu) then

      do k = skst+1, sken
      do ij = ijtstr, ijtend
         if ( w(ij, k) >= 0.0d0 ) then
            si(ij, k) = s(ij, k-1)
         else if ( w(ij, k) < 0.0d0 ) then
            si(ij, k) = s(ij, k)
         end if
         flx(ij, k)   = ( 1.0d0 - phii(k) ) * w(ij, k) * si(ij, k)
      end do
      end do

      do ij = ijtstr, ijtend
      if ( w(ij, sken+1) >= 0.0d0 ) then
         si(ij, sken+1) = s(ij, sken)
      else if ( w(ij, sken+1) < 0.0d0 ) then
         si(ij, sken+1) = 0.0d0
      end if
         flx(ij, sken+1) = ( 1.0d0 - phii(sken+1) ) &
     &                   * w(ij, sken+1) * si(ij, sken+1)
      end do

      end if

!     -- Advection : Quickest --

      if (sadq) then

      do k = skst+2, sken-1
      do ij = ijtstr, ijtend
         if ( w(ij, k) >= 0.0d0 ) then
            c2(ij, k) = ( 1.0d0 / ( dlzi(k-1) + dlzi(k) ) ) &
     &         * ( ( ( s(ij, k-2) - s(ij, k-1) ) / dlzi(k-1) ) &
     &           - ( ( s(ij, k-1) - s(ij, k)   ) / dlzi(k)   ) )
         else if ( w(ij, k) < 0.0d0 ) then
            c2(ij, k) = ( 1.0d0 / ( dlzi(k) + dlzi(k+1) ) ) &
     &         * ( ( ( s(ij, k-1) - s(ij, k) ) / dlzi(k)   ) &
     &           - ( ( s(ij, k) - s(ij, k+1) ) / dlzi(k+1) ) )
         end if
            c1(ij, k) = ( ( s(ij, k) - s(ij, k-1) ) / dlzi(k) ) &
     &         + ( ( ( dlz(k-1) - dlz(k) ) / 2.0d0 ) * c2(ij, k) )
            c0(ij, k) = ( ( dlz(k) * s(ij, k-1) + dlz(k-1) * s(ij, k) ) &
     &                / ( 2.0d0 * dlzi(k) ) ) &
     &         - ( ( ( dlz(k-1) * dlz(k) ) / 4.0d0 ) * c2(ij, k) )
      end do
      end do

      do ij = ijtstr, ijtend
      if ( w(ij, skst+1) >= 0.0d0 ) then
         c2(ij, skst+1) = ( 1.0d0 / ( dlzi(skst) + dlzi(skst+1) ) ) &
     &      * ( ( ( 0.0d0       - s(ij, skst)   ) / dlzi(skst) ) &
     &        - ( ( s(ij, skst) - s(ij, skst+1) ) / dlzi(skst+1) ) )
      else if ( w(ij, skst+1) < 0.0d0 ) then
         c2(ij, skst+1) = ( 1.0d0 / ( dlzi(skst+1) + dlzi(skst+2) ) ) &
     &      * ( ( ( s(ij, skst)   - s(ij, skst+1) ) / dlzi(skst+1) ) &
     &        - ( ( s(ij, skst+1) - s(ij, skst+2) ) / dlzi(skst+2) ) )
      end if
         c1(ij, skst+1) = ( ( s(ij, skst+1) - s(ij, skst) ) &
     &                      / dlzi(skst+1) ) &
     &      + ( ( ( dlz(skst) - dlz(skst+1) ) / 2.0d0 ) &
     &          * c2(ij, skst+1) )
         c0(ij, skst+1) = ( ( dlz(skst+1) * s(ij, skst) &
     &                      + dlz(skst)   * s(ij, skst+1) ) &
     &                  / ( 2.0d0 * dlzi(skst+1) ) ) &
     &      - ( ( ( dlz(skst) * dlz(skst+1) ) / 4.0d0 ) &
     &          * c2(ij, skst+1) )
      end do

      do ij = ijtstr, ijtend
      if ( w(ij, sken) >= 0.0d0 ) then
         c2(ij, sken) = ( 1.0d0 / ( dlzi(sken-1) + dlzi(sken) ) ) &
     &      * ( ( ( s(ij, sken-2) - s(ij, sken-1) ) / dlzi(sken-1) ) &
     &        - ( ( s(ij, sken-1) - s(ij, sken)   ) / dlzi(sken)   ) )
      else if ( w(ij, sken) < 0.0d0 ) then
         c2(ij, sken) = ( 1.0d0 / ( dlzi(sken) + dlzi(sken+1) ) ) &
     &      * ( ( ( s(ij, sken-1) - s(ij, sken) ) / dlzi(sken)   ) &
     &        - ( ( s(ij, sken)   - 0.0d0       ) / dlzi(sken+1) ) )
      end if
         c1(ij, sken) = ( ( s(ij, sken) - s(ij, sken-1) ) / dlzi(sken) ) &
     &      + ( ( ( dlz(sken-1) - dlz(sken) ) / 2.0d0 ) &
     &          * c2(ij, sken) )
         c0(ij, sken) = ( ( dlz(sken)   * s(ij, sken-1) &
     &                    + dlz(sken-1) * s(ij, sken)   ) &
     &                / ( 2.0d0 * dlzi(sken) ) ) &
     &      - ( ( ( dlz(sken-1) * dlz(sken) ) / 4.0d0 ) &
     &          * c2(ij, sken) )
      end do

      do ij = ijtstr, ijtend
      if ( w(ij, sken+1) >= 0.0d0 ) then
         c2(ij, sken+1) = ( 1.0d0 / ( dlzi(sken) + dlzi(sken+1) ) ) &
     &      * ( ( ( s(ij, sken-1) - s(ij, sken) ) / dlzi(sken)   ) &
     &        - ( ( s(ij, sken)   - s(ij, sken) ) / dlzi(sken+1) ) )
      else if ( w(ij, sken+1) < 0.0d0 ) then
         c2(ij, sken) = ( 1.0d0 / ( dlzi(sken+1) + dlzi(sken+1) ) ) &
     &      * ( ( ( s(ij, sken) - 0.0d0 ) / dlzi(sken+1) ) &
     &        - ( ( 0.0d0       - 0.0d0 ) / dlzi(sken+1) ) )
      end if
      if ( w(ij, sken+1) >= 0.0d0 ) then
         c1(ij, sken+1) = ( ( s(ij, sken) - s(ij, sken) ) &
     &                    / dlzi(sken+1) )
         c0(ij, sken+1) = ( ( dlz(sken) * s(ij, sken) &
     &                      + dlz(sken) * s(ij, sken) ) &
     &                  / ( 2.0d0 * dlzi(sken+1) ) ) &
     &      - ( ( ( dlz(sken) * dlz(sken) ) / 4.0d0 ) &
     &          * c2(ij, sken+1) )
      else if ( w(ij, sken+1) < 0.0d0 ) then
         c1(ij, sken+1) = ( ( 0.0d0 - s(ij, sken) ) &
     &                    / dlzi(sken+1) )
         c0(ij, sken+1) = ( ( dlz(sken) * s(ij, sken) &
     &                      + dlz(sken) * 0.0d0 ) &
     &                  / ( 2.0d0 * dlzi(sken+1) ) ) &
     &      - ( ( ( dlz(sken) * dlz(sken) ) / 4.0d0 ) &
     &          * c2(ij, sken+1) )
      end if
      end do

      do k = skst+1, sken+1
      do ij = ijtstr, ijtend
         si(ij, k) = c0(ij, k) &
     &   - ( c1(ij, k) * ( ( w(ij, k) * tssec ) / 2.0d0 ) ) &
     &   + ( c2(ij, k) * ( ( ( ( w(ij, k) * tssec ) ** 2.0d0 ) / 3.0d0 ) &
     &             - ( ( ( dlzi(k) ** 2.0d0 ) / 12.0d0 ) ) ) )
      end do
      end do

!     -- Advection : Quickest Flux limiter --

      if (sadqf) then

      do ij = ijtstr, ijtend
      if ( w(ij, skst+1) >= 0.0d0 ) then
      fldel(skst+1) = s(ij, skst+1) - 0.0d0
      flgam(skst+1) = ( ( s(ij, skst+1) - s(ij, skst) ) &
     &                  / dlzi(skst+1) ) &
     &              - ( ( s(ij, skst)   - 0.0d0       ) &
     &                  / dlzi(skst)   )
      if ( abs( flgam(skst+1)*dlz(skst) ) >= abs( fldel(skst+1) ) ) then
         si(ij, skst+1) = s(ij, skst)
      else if ( fldel(skst+1) >= 0 ) then
         fltmp = 0.0d0 &
     &         + ( s(ij, skst) - 0.0d0 ) / w(ij, skst+1) * dlzi(skst+1)
         si(ij, skst+1) = &
     &      min( si(ij, skst+1), min( fltmp, s(ij, skst+1) ) )
         si(ij, skst+1) = max( si(ij, skst+1), s(ij, skst) )
      else if ( fldel(skst+1) <  0 ) then
         fltmp = 0.0d0 &
     &         + ( s(ij, skst) - 0.0d0 ) / w(ij, skst+1) * dlzi(skst+1)
         si(ij, skst+1) = &
     &      max( si(ij, skst+1), max( fltmp, s(ij, skst+1) ) )
         si(ij, skst+1) = min( si(ij, skst+1), s(ij, skst) )
      end if
      else if ( w(ij, skst+1) < 0.0d0 ) then
      fldel(skst+1) = s(ij, skst) - s(ij, skst+2)
      flgam(skst+1) = ( ( s(ij, skst)   - s(ij, skst+1) ) &
     &                  / dlzi(skst+1) ) &
     &              - ( ( s(ij, skst+1) - s(ij, skst+2) ) &
     &                  / dlzi(skst+2) )
      if ( abs( flgam(skst+1)*dlz(skst+1) ) &
     &  >= abs( fldel(skst+1) ) ) then
         si(ij, skst+1) = s(ij, skst+1)
      else if ( fldel(skst+1) >= 0 ) then
         fltmp = s(ij, skst+2) &
     &- ( s(ij, skst+1) - s(ij, skst+2) ) / w(ij, skst+1) * dlzi(skst+1)
         si(ij, skst+1) = &
     &      min( si(ij, skst+1), min( fltmp, s(ij, skst) ) )
         si(ij, skst+1) = max( si(ij, skst+1), s(ij, skst+1) )
      else if ( fldel(skst+1) <  0 ) then
         fltmp = s(ij, skst+2) &
     &- ( s(ij, skst+1) - s(ij, skst+2) ) / w(ij, skst+1) * dlzi(skst+1)
         si(ij, skst+1) = &
     &      max( si(ij, skst+1), max( fltmp, s(ij, skst) ) )
         si(ij, skst+1) = min( si(ij, skst+1), s(ij, skst+1) )
      end if
      end if
      end do

      do k = skst+2, sken-1
      do ij = ijtstr, ijtend
      if ( w(ij, k) >= 0.0d0 ) then
         fldel(k) = s(ij, k) - s(ij, k-2)
         flgam(k) = ( ( s(ij, k)   - s(ij, k-1) ) / dlzi(k)   ) &
     &            - ( ( s(ij, k-1) - s(ij, k-2) ) / dlzi(k-1) )
         if ( abs( flgam(k)*dlz(k-1) ) >= abs( fldel(k) ) ) then
            si(ij, k) = s(ij, k-1)
         else if ( fldel(k) >= 0 ) then
            fltmp = s(ij, k-2) &
     &            + ( s(ij, k-1) - s(ij, k-2) ) / w(ij, k) * dlzi(k)
            si(ij, k) = &
     &         min( si(ij, k), min( fltmp, s(ij, k) ) )
            si(ij, k) = max( si(ij, k), s(ij, k-1) )
         else if ( fldel(k) <  0 ) then
            fltmp = s(ij, k-2) &
     &            + ( s(ij, k-1) - s(ij, k-2) ) / w(ij, k) * dlzi(k)
            si(ij, k) = &
     &         max( si(ij, k), max( fltmp, s(ij, k) ) )
            si(ij, k) = min( si(ij, k), s(ij, k-1) )
         end if
      else if ( w(ij, k) < 0.0d0 ) then
         fldel(k) = s(ij, k-1) - s(ij, k+1)
         flgam(k) = ( ( s(ij, k-1) - s(ij, k)   ) / dlzi(k)   ) &
     &            - ( ( s(ij, k)   - s(ij, k+1) ) / dlzi(k+1) )
         if ( abs( flgam(k)*dlz(k) ) >= abs( fldel(k) ) ) then
            si(ij, k) = s(ij, k)
         else if ( fldel(k) >= 0 ) then
            fltmp = s(ij, k+1) &
     &            - ( s(ij, k) - s(ij, k+1) ) / w(ij, k) * dlzi(k)
            si(ij, k) = &
     &         min( si(ij, k), min( fltmp, s(ij, k-1) ) )
            si(ij, k) = max( si(ij, k), s(ij, k) )
         else if ( fldel(k) <  0 ) then
            fltmp = s(ij, k+1) &
     &            - ( s(ij, k) - s(ij, k+1) ) / w(ij, k) * dlzi(k)
            si(ij, k) = &
     &         max( si(ij, k), max( fltmp, s(ij, k-1) ) )
            si(ij, k) = min( si(ij, k), s(ij, k) )
         end if
      end if
      end do
      end do

      do ij = ijtstr, ijtend
      if ( w(ij, sken) >= 0.0d0 ) then
      fldel(sken) = s(ij, sken) - s(ij, sken-2)
      flgam(sken) = ( ( s(ij, sken)   - s(ij, sken-1) ) / dlzi(sken)   ) &
     &            - ( ( s(ij, sken-1) - s(ij, sken-2) ) / dlzi(sken-1) )
      if ( abs( flgam(sken)*dlz(sken-1) ) >= abs( fldel(sken) ) ) then
         si(ij, sken) = s(ij, sken-1)
      else if ( fldel(sken) >= 0 ) then
         fltmp = s(ij, sken-2) &
     &+ ( s(ij, sken-1) - s(ij, sken-2) ) / w(ij, sken) * dlzi(sken)
         si(ij, sken) = &
     &      min( si(ij, sken), min( fltmp, s(ij, sken) ) )
         si(ij, sken) = max( si(ij, sken), s(ij, sken-1) )
      else if ( fldel(sken) <  0 ) then
         fltmp = s(ij, sken-2) &
     &+ ( s(ij, sken-1) - s(ij, sken-2) ) / w(ij, sken) * dlzi(sken)
         si(ij, sken) = &
     &      max( si(ij, sken), max( fltmp, s(ij, sken) ) )
         si(ij, sken) = min( si(ij, sken), s(ij, sken-1) )
      end if
      else if ( w(ij, sken) < 0.0d0 ) then
      fldel(sken) = s(ij, sken-1) - 0.0d0
      flgam(sken) = ( ( s(ij, sken-1) - s(ij, sken) ) / dlzi(sken)   ) &
     &            - ( ( s(ij, sken)   - 0.0d0       ) / dlzi(sken+1) )
      if ( abs( flgam(sken)*dlz(sken) ) >= abs( fldel(sken) ) ) then
         si(ij, sken) = s(ij, sken)
      else if ( fldel(sken) >= 0 ) then
         fltmp = 0.0d0 &
     &- ( s(ij, sken) - 0.0d0 ) / w(ij, sken) * dlzi(sken)
         si(ij, sken) = &
     &      min( si(ij, sken), min( fltmp, s(ij, sken-1) ) )
         si(ij, sken) = max( si(ij, sken), s(ij, sken) )
      else if ( fldel(sken) <  0 ) then
         fltmp = 0.0d0 &
     &- ( s(ij, sken) - 0.0d0 ) / w(ij, sken) * dlzi(sken)
         si(ij, sken) = &
     &      max( si(ij, sken), max( fltmp, s(ij, sken-1) ) )
         si(ij, sken) = min( si(ij, sken), s(ij, sken) )
      end if
      end if
      end do

      do ij = ijtstr, ijtend
      if ( w(ij, sken+1) >= 0.0d0 ) then
      fldel(sken+1) = s(ij, sken) - s(ij, sken-1)
      flgam(sken+1) = ( ( s(ij, sken) - s(ij, sken)   ) / dlzi(sken+1) ) &
     &              - ( ( s(ij, sken) - s(ij, sken-1) ) / dlzi(sken)   )
      if ( abs( flgam(sken+1)*dlz(sken) ) >= abs( fldel(sken+1) ) ) then
         si(ij, sken+1) = s(ij, sken)
      else if ( fldel(sken+1) >= 0 ) then
         fltmp = s(ij, sken-1) &
     &+ ( s(ij, sken) - s(ij, sken-1) ) / w(ij, sken+1) * dlzi(sken+1)
         si(ij, sken+1) = &
     &      min( si(ij, sken+1), min( fltmp, s(ij, sken) ) )
         si(ij, sken+1) = max( si(ij, sken+1), s(ij, sken) )
      else if ( fldel(sken+1) <  0 ) then
         fltmp = s(ij, sken-1) &
     &+ ( s(ij, sken) - s(ij, sken-1) ) / w(ij, sken+1) * dlzi(sken+1)
         si(ij, sken+1) = &
     &      max( si(ij, sken+1), max( fltmp, s(ij, sken) ) )
         si(ij, sken+1) = min( si(ij, sken+1), s(ij, sken) )
      end if
      else if ( w(ij, sken+1) < 0.0d0 ) then
      fldel(sken+1) = s(ij, sken) - 0.0d0
      flgam(sken+1) = ( ( s(ij, sken) - 0.0d0 ) / dlzi(sken+1) ) &
     &              - ( ( 0.0d0       - 0.0d0 ) / dlzi(sken+1) )
      if ( abs( flgam(sken+1)*dlz(sken) ) >= abs( fldel(sken+1) ) ) then
         si(ij, sken+1) = 0.0d0
      else if ( fldel(sken+1) >= 0 ) then
         fltmp = 0.0d0 &
     &- ( 0.0d0 - 0.0d0 ) / w(ij, sken+1) * dlzi(sken+1)
         si(ij, sken+1) = &
     &      min( si(ij, sken+1), min( fltmp, s(ij, sken) ) )
         si(ij, sken+1) = max( si(ij, sken+1), 0.0d0 )
      else if ( fldel(sken+1) <  0 ) then
         fltmp = 0.0d0 &
     &- ( 0.0d0 - 0.0d0 ) / w(ij, sken+1) * dlzi(sken+1)
         si(ij, sken+1) = &
     &      max( si(ij, sken+1), max( fltmp, s(ij, sken) ) )
         si(ij, sken+1) = min( si(ij, sken+1), 0.0d0 )
      end if
      end if
      end do

      end if

      do k = skst+1, sken+1
      do ij = ijtstr, ijtend
         flx(ij, k) = ( 1.0d0 - phii(k) ) * w(ij, k) * si(ij, k)
      end do
      end do

      end if

!     -- Diffusion : Crank nicolson

      if (sdcn) then

      do k = skst, sken
         eal(k) = ( ( 1.0d0 - phi(k) ) / tssec )
         ebe(k) = ( ( ( 1.0d0 - phii(k) )   * dbi(k) ) &
     &            / ( 2.0d0 * dlz(k) * dlzi(k)   ) )
         ega(k) = ( ( ( 1.0d0 - phii(k+1) ) * dbi(k+1) ) &
     &            / ( 2.0d0 * dlz(k) * dlzi(k+1) ) )
      end do

      do k = skst+1, sken-1
         ea(k)     = - ebe(k)
         eb(k)     = eal(k) + ebe(k) + ega(k)
         ec(k)     = - ega(k)
      end do

         ea(skst)  = - ebe(skst)
         eb(skst)  = eal(skst) + ega(skst)
         ec(skst)  = - ega(skst)
         ea(sken)  = - ebe(sken)
         eb(sken)  = eal(sken) + ebe(sken)
         ec(sken)  = - ega(sken)

      do ij = ijtstr, ijtend
         ed(ij, skst) = &
     &      + ( eal(skst) - ega(skst) ) * s(ij, skst) &
     &      + ega(skst) * s(ij, skst+1) &
     &      + ( ( flx(ij, skst) - flx(ij, skst+1) ) / dlz(skst) ) &
     &      + ( ( 1.0d0 - phi(skst) ) / 3.1536d+7 )
      end do

      do k = skst+1, sken-1
      do ij = ijtstr, ijtend
         ed(ij, k) = &
     &      + ebe(k) * s(ij, k-1) &
     &      + ( eal(k) - ebe(k) - ega(k) ) * s(ij, k) &
     &      + ega(k) * s(ij, k+1) &
     &      + ( ( flx(ij, k) - flx(ij, k+1) ) / dlz(k) ) &
     &      + ( ( 1.0d0 - phi(k) ) / 3.1536d+7 )
      end do
      end do

      do ij = ijtstr, ijtend
         ed(ij, sken) = &
     &      + ebe(sken) * s(ij, sken-1) &
     &      + ( eal(sken) - ebe(sken) ) * s(ij, sken) &
     &      + ( ( flx(ij, sken) - flx(ij, sken+1) ) / dlz(sken) ) &
     &      + ( ( 1.0d0 - phi(sken) ) / 3.1536d+7 )
      end do

      call calc_tridiagonal_solid( ea, eb, ec, ed, ex )

      end if

!     -- Diffusion : Implicit

      if (sdii) then

      do k = skst, sken
         eal(k) = ( ( 1.0d0 - phi(k) ) / tssec )
         ebe(k) = ( ( ( 1.0d0 - phii(k) )   * dbi(k) ) &
     &            / ( dlz(k) * dlzi(k)   ) )
         ega(k) = ( ( ( 1.0d0 - phii(k+1) ) * dbi(k+1) ) &
     &            / ( dlz(k) * dlzi(k+1) ) )
      end do

      do k = skst+1, sken-1
         ea(k)     = - ebe(k)
         eb(k)     = eal(k) + ebe(k) + ega(k)
         ec(k)     = - ega(k)
      end do

         ea(skst)  = - ebe(skst)
         eb(skst)  = eal(skst) + ega(skst)
         ec(skst)  = - ega(skst)
         ea(sken)  = - ebe(sken)
         eb(sken)  = eal(sken) + ebe(sken)
         ec(sken)  = - ega(sken)

      do k = skst, sken
      do ij = ijtstr, ijtend
         ed(ij, k) = &
     &      + eal(k) * s(ij, k) &
     &      + ( ( flx(ij, k) - flx(ij, k+1) ) / dlz(k) ) &
     &      + ( ( 1.0d0 - phi(k) ) / 3.1536d+7 )
      end do
      end do

      call calc_tridiagonal_solid( ea, eb, ec, ed, ex )

      end if

!     -- t = t + dt --

      do k = skst, sken
      do ij = ijtstr, ijtend
         sn(ij, k) = max( ex(ij, k), 0.0d0 )
         ex(ij, k) = 0.0d0
      end do
      end do

      end subroutine calc_step_solid_age

! =====================================================================

      subroutine calc_tridiagonal_solid( ea, eb, ec, ed, ex )

      double precision, intent(in)  :: ea(sdzdim), eb(sdzdim), ec(sdzdim)
      double precision, intent(in)  :: ed(nxy,sdzdim)
      double precision, intent(out) :: ex(nxy,sdzdim)

      double precision :: ew(nxy,sdzdim)
      double precision :: bet(nxy)

      integer :: ij, k

      do ij = ijtstr, ijtend
         bet(ij)      = eb(skst)
         ex(ij, skst) = ed(ij, skst) / bet(ij)
      end do

      do k = skst+1, sken
      do ij = ijtstr, ijtend
         ew(ij, k) = ec(k-1) / bet(ij)
         bet(ij)   = eb(k) - ea(k) * ew(ij, k)
         ex(ij, k) = (ed(ij, k) - ea(k) * ex(ij, k-1) ) / bet(ij)
      end do
      end do

      do k = sken-1, skst, -1
      do ij = ijtstr, ijtend
         ex(ij, k) = ex(ij, k) - ew(ij, k+1) * ex(ij, k+1)
      end do
      end do

      end subroutine calc_tridiagonal_solid

end module module_digda_step_solid
