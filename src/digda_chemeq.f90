!######################################################################

! --- information -----------------------------------------------------
!
!  Calculate [co2(g)] from chemical equilibrium in ocean
!
!  HISTORY
!     '04.04.23  A.Oka :
!     '06.06.03  A.Oka : ocean density [mMol/m^3] --> [Mol/kg]
!
! ---------------------------------------------------------------------

module globals_digda_chemeq
  use module_zcfile
  use module_zcphys
  use module_digda_grid
  use module_zchem
  implicit none

  double precision :: phini = 8.0
  double precision :: dphini = 8.0
  double precision :: dph = 3.0
  double precision :: alkerr = 0.1
  double precision :: rocn = 1.0245 !! ocean density [kg/L]
  integer :: itrmax = 10000

  namelist /nmsdchem/ phini, dphini, dph, alkerr, rocn, itrmax

end module globals_digda_chemeq

! =====================================================================

module module_digda_chemeq
  use globals_digda_chemeq
  implicit none

  double precision :: ph(nxy, sdzdim)
  double precision :: cu, rcu, pm
  logical :: ofirst

contains

  subroutine set_digda_chemeq

    integer :: ij, k

    write(nfstdo, *) '*** Sediment : set_digda_chemeq ***'
    rewind(nfparm)
    read(nfparm, nmsdchem, end=99)
99 continue
    write(nfstdo, nmsdchem)

    do ij = 1, nxy
      do k = 1, sdzdim
        ph(ij, k) = phini
      end do
    end do

    pm  = 1.0d+6          !! [atm] --> [ppmv]

    cu  = 1.0d-6 / rocn   !! [mMol/m^3] --> [Mol/kg]
    rcu = 1.0d+6 * rocn   !! [Mol/kg] --> [mMol/m^3]
    alkerr = alkerr * cu

    ofirst = .true.

  end subroutine set_digda_chemeq

  subroutine calc_digda_chemeq( temp, salt, ctin, atin, & !! IN
                              & sitin, ptin,            & !! IN
                              & co32                    ) !! OUT

    double precision, intent(in)  :: temp(nxy)         !! Temperature
    double precision, intent(in)  :: salt(nxy)         !! Salinity
    double precision, intent(in)  :: ctin(nxy, sdzdim) !! Total carbon [mMol/m^3]
    double precision, intent(in)  :: atin(nxy, sdzdim) !! Alkalinity [mMol/m^3]
    double precision, intent(in)  :: sitin(nxy)        !! Total silicon [mMol/m^3]
    double precision, intent(in)  :: ptin(nxy)         !! Total phosphorous [mMol/m^3]
    double precision, intent(out) :: co32(nxy, sdzdim) !! [co32-] [mMol/m^3]

    double precision :: pco2o(nxy)        !! [co2(g)] [ppmv]
    double precision :: alp(nxy)          !! Solubility of co2 [mMol/m^3/ppmv]
    double precision :: ct(nxy, sdzdim)   !! Total carbon [Mol/kg]
    double precision :: at(nxy, sdzdim)   !! Alkalinity [Mol/kg]
    double precision :: sit(nxy)          !! Total silicon [Mol/kg]
    double precision :: pt(nxy)           !! Total phosphorous [Mol/kg]
    double precision :: bt(nxy)           !! Total boric acid [Mol/kg]
    double precision :: so4t(nxy)         !! Total SO4 [Mol/kg]
    double precision :: flt(nxy)          !! Total fluorosis [Mol/kg]
    double precision :: co2(nxy, sdzdim)  !! [co2]   [mMol/m^3]
    double precision :: hco3(nxy, sdzdim) !! [hco3-] [mMol/m^3]

    double precision :: k0(nxy),  k1(nxy),   k2(nxy)
    double precision :: kb(nxy),  kw(nxy)
    double precision :: k1p(nxy), k2p(nxy),  k3p(nxy)
    double precision :: ksi(nxy), kso4(nxy), kf(nxy)

    double precision :: h, x, xlw, xhi, f, delt

    integer :: ij, k, itr

    do ij = ijtstr, ijtend
    ! -- [mMol/m^3]-->[Mol/kg] --
      do k = skst, sken
        ct(ij, k)  = ctin(ij, k)  * cu
        at(ij, k)  = atin(ij, k)  * cu
      end do
      sit(ij) = sitin(ij) * cu
      pt(ij)  = ptin(ij)  * cu
    end do

    do ij = ijtstr, ijtend
   !if ( temp(ij) /= 0.0d0 .and. salt(ij) /= 0.0d0 ) then
    if ( temp(ij) /= rmiss .and. salt(ij) /= rmiss ) then
    ! -- Chemical equilibrium constant --
      k0(ij)   = fk0( temp(ij)+kelvin, salt(ij) )
      k1(ij)   = fk1( temp(ij)+kelvin, salt(ij) )
      k2(ij)   = fk2( temp(ij)+kelvin, salt(ij) )
      kb(ij)   = fkb( temp(ij)+kelvin, salt(ij) )
      kw(ij)   = fkw( temp(ij)+kelvin, salt(ij) )
      k1p(ij)  = fk1p( temp(ij)+kelvin, salt(ij) )
      k2p(ij)  = fk2p( temp(ij)+kelvin, salt(ij) )
      k3p(ij)  = fk3p( temp(ij)+kelvin, salt(ij) )
      ksi(ij)  = fksi( temp(ij)+kelvin, salt(ij) )
      kso4(ij) = fkso4( temp(ij)+kelvin, salt(ij) )
      kf(ij)   = fkf( temp(ij)+kelvin, salt(ij), kso4(ij) )
!     -- Total [Mol/kg] --
      bt(ij)   = fbt( salt(ij) )
      so4t(ij) = fso4t( salt(ij) )
      flt(ij)  = fflt( salt(ij) )
    end if
    end do

    if (ofirst) then
      delt = dphini
      ofirst = .false.
    else
      delt = dph
    end if

    do k = skst, sken
    do ij = ijtstr, ijtend
   !if ( temp(ij) /= 0.0d0 .and. salt(ij) /= 0.0d0 ) then
    if ( temp(ij) /= rmiss .and. salt(ij) /= rmiss ) then

      xlw = ph(ij, k) - delt
      xhi = ph(ij, k) + delt

      !-- Iteration : bisection method --
      do itr = 1, itrmax
        x = 0.5d0 * ( xlw + xhi )
        h = 10.0d0 ** ( - x )
       !if ( amskt(ij, kstr). eq. 0.0d0 ) goto 999
       !if ( temp(ij) .eq. 0.0d0 .and. salt(ij) .eq. 0.0d0 ) goto 999
        if ( temp(ij) .eq. rmiss .and. salt(ij) .eq. rmiss ) goto 999
        f = falk( &
     &      h, &
     &      k1(ij),  k2(ij), &
     &      kb(ij),  kw(ij), &
     &      k1p(ij), k2p(ij),  k3p(ij), &
     &      ksi(ij), kso4(ij), kf(ij), &
     &      ct(ij, k),  bt(ij),   pt(ij), &
     &      sit(ij), so4t(ij), flt(ij) )
        if ( abs(f-at(ij, k)).lt.alkerr ) goto 999
        if ( ( f - at(ij, k) ).gt.0.d0 ) then
          xhi = x
        else
          xlw = x
        end if
      end do

      write(nfstdo,*) ' ITERATION COUNT OVER (SEDIMENT) :', IJ
      write(nfstdo,*) '    PH  = ', x
      write(nfstdo,*) '    F   = ', f * rcu
      write(nfstdo,*) '    ALK = ', atin(ij, k)
      write(nfstdo,*) '    F - ALK = ', f * rcu - atin(ij, k)
      write(nfstdo,*) '            < ALKERR = ', alkerr * rcu
      write(nfstdo,*) '   [Info] '
     !write(nfstdo,*) '    at  I = ', mod(IJ,nxdim) - istr + 1
     !write(nfstdo,*) '        J = ', IJ/nxdim - jstr + 2
     !write(nfstdo,*) '        K = ', K - skst + 1
      write(nfstdo,*) '    TEM = ', temp(ij)
      write(nfstdo,*) '    SAL = ', salt(ij)
      write(nfstdo,*) '    DIC = ', ctin(ij, k)

      xlw = ph(ij, k) - delt
      xhi = ph(ij, k) + delt
      do itr = 1, itrmax
        x = 0.5d0 * ( xlw + xhi )
        h = 10.0d0 ** ( - x )
        f = falk( &
     &      h, &
     &      k1(ij),  k2(ij), &
     &      kb(ij),  kw(ij), &
     &      k1p(ij), k2p(ij),  k3p(ij), &
     &      ksi(ij), kso4(ij), kf(ij), &
     &      ct(ij, k),  bt(ij),   pt(ij), &
     &      sit(ij), so4t(ij), flt(ij) )
        if ( abs(f-at(ij, k)).lt.alkerr ) goto 999
        if ( ( f - at(ij, k) ).gt.0.d0 ) then
          xhi = x
        else
          xlw = x
        end if
      end do

!     STOP

999 continue
      ph(ij, k) = x

    end if
    end do
    end do

    do k = skst, sken
    do ij = ijtstr, ijtend
   !if ( temp(ij) /= 0.0d0 .and. salt(ij) /= 0.0d0 ) then
    if ( temp(ij) /= rmiss .and. salt(ij) /= rmiss ) then
      x = 10.d0 ** ( - ph(ij, k) )
!     pco2o(ij) = ct(ij) *x*x &
!  &            / ( x*x + x*k1(ij) + k1(ij)*k2(ij) ) &
!  &            / k0(ij) &
!  &            * pm      !! [atm]-->[ppmv]
!     alp(ij) = k0(ij)      !! [Mol/kg/atm] = [mMol/m^3/ppmv]
      ! -- co2 --
      co2(ij, k)  = ct(ij, k) *x*x &
   &              / ( x*x + x*k1(ij) + k1(ij)*k2(ij) ) / cu
      ! -- BI-CARBONATE ION --
      hco3(ij, k) = k1(ij) * co2(ij, k) / x
      ! -- CARBONATE ION --
      co32(ij, k) = k2(ij) * hco3(ij, k) / x
    end if
    end do
    end do

  end subroutine calc_digda_chemeq

end module module_digda_chemeq
