!######################################################################

! --- information -----------------------------------------------------
!
!  Three-dimensional sediment model
!
!  HISTORY
!     '14.08.01  H.Kobayashi :
!
! ---------------------------------------------------------------------

module module_digda_cprdc
  use netcdf_utilities
  use gtool_utilities
  use module_zcfile
  use module_ucloc
  use module_digda_grid
  use module_digda_mpi
  use module_digda_parameter
  use module_digda_timestep
  use module_digda_clay
  use module_digda_pom
  use module_digda_carbonate
  use module_digda_silicate
  use module_digda_step_solid
  use module_digda_step_liquid
  use module_digda_chemeq
  use module_digda_time_mean
  implicit none

  logical :: sdrest = .true.
  logical :: sdofirst = .true.
  namelist /nmsdrst/ sdrest

  logical :: deiarc = .true.  !! Archer's DEI formulation ?
  namelist /nmsddear/ deiarc

  logical :: calc13 = .true. !! calculate carbon-13 ?
  logical :: calc14 = .true. !! calculate carbon-14 ?
  namelist /nmsdiso/ calc13, calc14

  ! etc
  character(16) :: sedtnm(sedtdm)
  character(128) :: cfile, citem, ctitle, cunit
  character(len=16) :: chead(64)

contains

  subroutine digda_cprdc( &
 &  bdepg, &
 &  bfrcg, &
 &  omaskg, &
 &  sdtx, &
 &  lon_values, &
 &  lat_values, &
 &  depth_values )

    ! --------------------------------------------------
    !   Forcing of Particles and Dissolved material
    !      Clay, POM, CaCO3, Opal, agetracer
    !      Phosphate, Dissolved Inorganic Carbon,
    !      Alkalinity, Dissolved Oxygen, Silicate
    ! --------------------------------------------------

    double precision, intent(in) :: bdepg(nxydim)
    double precision, intent(in) :: bfrcg(nxydim, benfdm)
    double precision, intent(in) :: omaskg(nxdim, nydim)

    double precision, intent(inout) :: sdtx(nxydim, sdzdim, sedtdm)
    double precision :: sdtxbuf(nxdim, nydim, sdzdim, sedtdm)

    ! --------------------------------------------------
    !   Particles
    !     Clay, POM, CaCO3, Opal, agetracer
    ! --------------------------------------------------

    ! Concentration
    !   clay [g/m^3], pom, caco3, opal [mMol/m^3], agetracer [yr]
    double precision :: scla(nxy, sdz),   sncla(nxy, sdz)
    double precision :: spom(nxy, sdz),   snpom(nxy, sdz)
    double precision :: scar(nxy, sdz),   sncar(nxy, sdz)
    double precision :: sopa(nxy, sdz),   snopa(nxy, sdz)
    double precision :: sage(nxy, sdz),   snage(nxy, sdz)
    double precision :: s13pom(nxy, sdz), sn13pom(nxy, sdz)
    double precision :: s14pom(nxy, sdz), sn14pom(nxy, sdz)
    double precision :: s13car(nxy, sdz), sn13car(nxy, sdz)
    double precision :: s14car(nxy, sdz), sn14car(nxy, sdz)
    double precision :: sclag(nxydim, sdzdim)
    double precision :: spomg(nxydim, sdzdim)
    double precision :: scarg(nxydim, sdzdim)
    double precision :: sopag(nxydim, sdzdim)
    double precision :: sageg(nxydim, sdzdim)
    double precision :: s13pomg(nxydim, sdzdim)
    double precision :: s14pomg(nxydim, sdzdim)
    double precision :: s13carg(nxydim, sdzdim)
    double precision :: s14carg(nxydim, sdzdim)

    ! Weight percent [%]
    double precision :: pcla(nxy, sdz)
    double precision :: ppom(nxy, sdz)
    double precision :: pcar(nxy, sdz)
    double precision :: popa(nxy, sdz)
    double precision :: wt(nxy, sdz)
!!!!double precision :: pclag(nxydim, sdzdim)
!!!!double precision :: ppomg(nxydim, sdzdim)
!!!!double precision :: pcarg(nxydim, sdzdim)
!!!!double precision :: popag(nxydim, sdzdim)
!!!!double precision :: wtg(nxydim, sdzdim)

    ! Dissolution rate [1/s], Dissolution term [mMol/m^3/s]
    double precision :: rcla(nxy, sdz),   trcla(nxy, sdz)
    double precision :: rpom(nxy, sdz),   trpom(nxy, sdz)
    double precision :: rcar(nxy, sdz),   trcar(nxy, sdz)
    double precision :: ropa(nxy, sdz),   tropa(nxy, sdz)
    double precision :: r13pom(nxy, sdz), tr13pom(nxy, sdz)
    double precision :: r14pom(nxy, sdz), tr14pom(nxy, sdz)
    double precision :: r13car(nxy, sdz), tr13car(nxy, sdz)
    double precision :: r14car(nxy, sdz), tr14car(nxy, sdz)

    ! --------------------------------------------------
    !   Dissolved material
    !     Phosphate, Dissolved Inorganic Carbon,
    !     Alkalinity, Dissolved Oxygen, Silicate
    ! --------------------------------------------------

    ! Concentration [mMol/m^3]
    double precision :: cnit(nxy, sdz),   cnnit(nxy, sdz)
    double precision :: cdic(nxy, sdz),   cndic(nxy, sdz)
    double precision :: calk(nxy, sdz),   cnalk(nxy, sdz)
    double precision :: cdo2(nxy, sdz),   cndo2(nxy, sdz)
    double precision :: csil(nxy, sdz),   cnsil(nxy, sdz)
    double precision :: c13dic(nxy, sdz), cn13dic(nxy, sdz)
    double precision :: c14dic(nxy, sdz), cn14dic(nxy, sdz)
    double precision :: cdc13(nxy, sdz),  cdc14(nxy, sdz)
    double precision :: cco32(nxy, sdz)
    double precision :: satco32(nxy)
    double precision :: satsil(nxy, sdz)
    double precision :: cnitg(nxydim, sdzdim)
    double precision :: cdicg(nxydim, sdzdim)
    double precision :: calkg(nxydim, sdzdim)
    double precision :: cdo2g(nxydim, sdzdim)
    double precision :: csilg(nxydim, sdzdim)
    double precision :: c13dicg(nxydim, sdzdim)
    double precision :: c14dicg(nxydim, sdzdim)
    double precision :: cdc13g(nxydim, sdzdim)
    double precision :: cdc14g(nxydim, sdzdim)
    double precision :: cco32g(nxydim, sdzdim)
    double precision :: satco32g(nxydim)
    double precision :: satsilg(nxydim, sdzdim)

    ! Concentration (bottom water)
    !   temperature [C], salinity [psu], the others [mMol/m^3]
    double precision :: bwtem(nxy)
    double precision :: bwsal(nxy)
    double precision :: bwnit(nxy)
    double precision :: bwdic(nxy)
    double precision :: bwalk(nxy)
    double precision :: bwdo2(nxy)
    double precision :: bwsil(nxy)
    double precision :: bw13dic(nxy)
    double precision :: bw14dic(nxy)
    double precision :: bwtemg(nxydim)
    double precision :: bwsalg(nxydim)
    double precision :: bwnitg(nxydim)
    double precision :: bwdicg(nxydim)
    double precision :: bwalkg(nxydim)
    double precision :: bwdo2g(nxydim)
    double precision :: bwsilg(nxydim)
    double precision :: bw13dicg(nxydim)
    double precision :: bw14dicg(nxydim)

    ! Concentration (total)
    !   total phosphorous [mMol/m^3], total silicon [mMol/m^3]
    double precision :: pt(nxy)
    double precision :: sit(nxy)

    ! Dissolution term [mMol/m^3/s]
    double precision :: trnit(nxy, sdz)
    double precision :: trdic(nxy, sdz)
    double precision :: tralk(nxy, sdz)
    double precision :: trdo2(nxy, sdz)
    double precision :: trsil(nxy, sdz)
    double precision :: tr13dic(nxy, sdz)
    double precision :: tr14dic(nxy, sdz)

    ! --------------------------------------------------
    !   Fluxes
    ! --------------------------------------------------

    ! Sedimentation rate
    double precision :: w(nxy, sdz+1)
    double precision :: wg(nxydim, sdzdim+1)

    ! Sea-bottom Flux
    !   clay [g/m^3*cm/s], the others [mMol/m^3*cm/s]
    double precision :: ficla(nxy),   focla(nxy)
    double precision :: fipom(nxy),   fopom(nxy)
    double precision :: ficar(nxy),   focar(nxy)
    double precision :: fiopa(nxy),   foopa(nxy)
    double precision :: fiage(nxy)
    double precision :: fi13pom(nxy), fo13pom(nxy)
    double precision :: fi13car(nxy), fo13car(nxy)
    double precision :: fi14pom(nxy), fo14pom(nxy)
    double precision :: fi14car(nxy), fo14car(nxy)

    double precision :: fiocla(nxy)
    double precision :: fiopom(nxy)
    double precision :: fiocar(nxy)
    double precision :: fioopa(nxy)

    double precision :: ficlag(nxydim),   foclag(nxydim)
    double precision :: fipomg(nxydim),   fopomg(nxydim)
    double precision :: ficarg(nxydim),   focarg(nxydim)
    double precision :: fiopag(nxydim),   foopag(nxydim)
    double precision :: fiageg(nxydim),   foageg(nxydim)
    double precision :: fi13pomg(nxydim), fo13pomg(nxydim)
    double precision :: fi13carg(nxydim), fo13carg(nxydim)
    double precision :: fi14pomg(nxydim), fo14pomg(nxydim)
    double precision :: fi14carg(nxydim), fo14carg(nxydim)

    double precision :: fioclag(nxydim)
    double precision :: fiopomg(nxydim)
    double precision :: fiocarg(nxydim)
    double precision :: fioopag(nxydim)

    ! Sea-bottom Flux [mMol/m^3*cm/s]
    double precision :: finit(nxy),   fonit(nxy)
    double precision :: fidic(nxy),   fodic(nxy)
    double precision :: fialk(nxy),   foalk(nxy)
    double precision :: fido2(nxy),   fodo2(nxy)
    double precision :: fisil(nxy),   fosil(nxy)
    double precision :: fi13dic(nxy), fo13dic(nxy)
    double precision :: fi14dic(nxy), fo14dic(nxy)

    ! Sediment-bottom Flux
    !   clay [g/m^3*cm/s], the others [mMol/m^3*cm/s]
    double precision :: bfcla(nxy)
    double precision :: bfpom(nxy)
    double precision :: bfcar(nxy)
    double precision :: bfopa(nxy)
    double precision :: bf13pom(nxy)
    double precision :: bf14pom(nxy)
    double precision :: bf13car(nxy)
    double precision :: bf14car(nxy)
    double precision :: bfclag(nxydim)
    double precision :: bfpomg(nxydim)
    double precision :: bfcarg(nxydim)
    double precision :: bfopag(nxydim)
    double precision :: bf13pomg(nxydim)
    double precision :: bf14pomg(nxydim)
    double precision :: bf13carg(nxydim)
    double precision :: bf14carg(nxydim)

    ! --------------------------------------------------
    !   Output
    ! --------------------------------------------------

    type(time_mean_type), allocatable :: mean_data(:)
    integer :: num_vars
    double precision :: scla_avg(nxy, sdz)
    double precision :: spom_avg(nxy, sdz)
    double precision :: scar_avg(nxy, sdz)
    double precision :: sopa_avg(nxy, sdz)
    double precision :: sage_avg(nxy, sdz)
    double precision :: cnit_avg(nxy, sdz)
    double precision :: cdic_avg(nxy, sdz)
    double precision :: calk_avg(nxy, sdz)
    double precision :: cdo2_avg(nxy, sdz)
    double precision :: csil_avg(nxy, sdz)
    double precision :: ficla_avg(nxy)
    double precision :: fipom_avg(nxy)
    double precision :: ficar_avg(nxy)
    double precision :: fiopa_avg(nxy)
    double precision :: focla_avg(nxy)
    double precision :: fopom_avg(nxy)
    double precision :: focar_avg(nxy)
    double precision :: foopa_avg(nxy)
    double precision :: fiocla_avg(nxy)
    double precision :: fiopom_avg(nxy)
    double precision :: fiocar_avg(nxy)
    double precision :: fioopa_avg(nxy)
    double precision :: bfcla_avg(nxy)
    double precision :: bfpom_avg(nxy)
    double precision :: bfcar_avg(nxy)
    double precision :: bfopa_avg(nxy)
   !double precision :: w_avg(nxy, sdz+1)
   !double precision :: cdc13_avg(nxy, sdz)
   !double precision :: cdc14_avg(nxy, sdz)

    ! --------------------------------------------------
    !   Others
    ! --------------------------------------------------

    ! MASK
    ! ocean depth [m]
    double precision :: bdep(nxy)

    ! Carbon isotope metrics
    double precision :: feoph, feosh, r13, r14, adcay, thalf, dltc13, dltc14

    ! Parallel computing
    integer, allocatable :: nrcvct2d(:),  ndispt2d(:)
    integer, allocatable :: nrcvct3d(:),  ndispt3d(:)
    integer, allocatable :: nrcvcw3d(:),  ndispw3d(:)
    double precision, allocatable :: vardt3d(:), vardt2d(:), vardw3d(:)
    double precision, allocatable :: varavgdt(:,:)

    ! etc
    character(len=64) :: filename
    character(len=16) :: varname, varunit
    double precision :: lon_values(nxdim)
    double precision :: lat_values(nydim)
    double precision :: depth_values(sdzdim)
    double precision, allocatable :: var2d_values(:,:,:)
    double precision, allocatable :: var3d_values(:,:,:,:)

    double precision :: integral_value, mean_value

    integer :: i, j, ij, k, l, ic
    integer :: tt, yy
!   integer, dimension(:),  :: yy_array
    integer(kind=4), dimension(:), allocatable :: yy_array
    integer :: ios

    !-- SEDIMENT TRACERS : setup --

    call clcstr('PRDC_INIT')

    if (use_mpi) then
      do irank = 0, nprocs-1
        call set_parallel_range( 1, nxydim, nprocs, irank,  &
                               & jjstr(irank), jjend(irank) )
      end do
    end if

    write(nfstdo, *) "'--------------------------------------------------'"
    write(nfstdo, *) "'-- DIGDA Prediction:                            --'"
    write(nfstdo, *) "'-- Started setting up model parameters          --'"
    write(nfstdo, *) "'--------------------------------------------------'"
    rewind(nfparm)
    read(nfparm, NML=nmsdhi, IOSTAT=ios)
    write(nfstdo, nmsdhi)
    rewind(nfparm)
    read(nfparm, NML=nmsdrst, IOSTAT=ios)
    write(nfstdo, nmsdrst)
    rewind(nfparm)
    read(nfparm, NML=nmsddear, IOSTAT=ios)
    write(nfstdo, nmsddear)
    rewind(nfparm)
    read(nfparm, NML=nmsdiso, IOSTAT=ios)
    write(nfstdo, nmsdiso)

    call initialize_timestep()

    sedtnm(1)  = 'scla'
    sedtnm(2)  = 'spom'
    sedtnm(3)  = 'scar'
    sedtnm(4)  = 'sopa'
    sedtnm(5)  = 'sage'
    sedtnm(6)  = 'cnit'
    sedtnm(7)  = 'cdic'
    sedtnm(8)  = 'calk'
    sedtnm(9)  = 'csil'
    sedtnm(10) = 'cdo2'
    if ( sedtdm == 16 ) then
    sedtnm(11) = 's13pom'
    sedtnm(12) = 's13car'
    sedtnm(13) = 'c13dic'
    sedtnm(14) = 's14pom'
    sedtnm(15) = 's14car'
    sedtnm(16) = 'c14dic'
    endif

    call set_digda_clay()
    call set_digda_pom()
    call set_digda_carbonate()
    call set_digda_silicate()
    call set_digda_chemeq()
    call set_digda_step_solid()
    call set_digda_step_liquid()

    write(nfstdo, *) "'--------------------------------------------------'"
    write(nfstdo, *) "'-- DIGDA Prediction:                            --'"
    write(nfstdo, *) "'-- Finished setting up model parameters         --'"
    write(nfstdo, *) "'--------------------------------------------------'"
    call flush(nfstdo)

    if (hstavg) then
      !-- Initialize time_mean --
      num_vars = 40  ! Number of variables to track
      allocate(mean_data(num_vars))
      ! Variable 1--10: 3D grid
      ! scla, spom, scar, sopa, sage
      ! cnit, cdic, calk, csil, cdo2
      do i = 1, 10
         call initialize_time_mean(mean_data(i), nxy, sdz, 3)
      end do
      ! Variable 11--26: 2D grid (nz = 1)
      ! ficla,  fipom,  ficar,  fiopa
      ! focla,  fopom,  focar,  foopa
      ! fiocla, fiopom, fiocar, fioopa
      ! bfcla,  bfpom,  bfcar,  bfopa
      do i = 11, 26
         call initialize_time_mean(mean_data(i), nxy, 1, 2)
      end do
      ! Variable 27--28: 3D grid
      ! cco32, satsil
      do i = 27, 28
         call initialize_time_mean(mean_data(i), nxy, sdz, 3)
      end do
      ! Variable 29: 2D grid
      ! satco32
      do i = 29, 29
         call initialize_time_mean(mean_data(i), nxy, 1, 2)
      end do
      ! Variable 30: 3D grid
      ! w
      do i = 30, 30
         call initialize_time_mean(mean_data(i), nxy, sdz+1, 3)
      end do
      if ( calc13 .and. sedtdm == 16 ) then
      ! Variable 31--34: 3D grid
      ! s13pom, s13car, c13dic, cdc13
      do i = 31, 34
         call initialize_time_mean(mean_data(i), nxy, sdz, 3)
      end do
      ! Variable 35--40: 2D grid
      ! fi13pom, fi13car, fo13pom, fo13car, bf13pom, bf13car
      do i = 35, 40
         call initialize_time_mean(mean_data(i), nxy, 1, 2)
      end do
      end if
      if ( calc14 .and. sedtdm == 16 ) then
      ! Variable 41--44: 3D grid
      ! s14pom, s14car, c14dic, cdc14
      do i = 41, 44
         call initialize_time_mean(mean_data(i), nxy, sdz, 3)
      end do
      ! Variable 45--50: 2D grid
      ! fi14pom, fi14car, fo14pom, fo14car, bf14pom, bf14car
      do i = 45, 50
         call initialize_time_mean(mean_data(i), nxy, 1, 2)
      end do
      end if
    end if

    !-- SEDIMENT TRACERS : setup (end) --

    ! -- Ocean depth --
    if (use_mpi) then
      do ij = 1, nxy
         bdep(ij) = max( bdepg(myrank*nxy+ij), 0.0d0 )
      end do
    else
      do ij = 1, nxy
         bdep(ij) = max( bdepg(ij), 0.0d0 )
      end do
    end if

    ! -- Convert unit --

    z0  = z0  * 1.0d+2                       !! [m] --> [cm]

    db0 = db0 * 1.0d+4 / 3.1536d+7           !! [m^2/yr] --> [cm^2/s]
    zb  = zb  * 1.0d+2                       !! [m] --> [cm]

    do k = skst, sken
       lz(k)   = lz(k)   * 1.0d+2            !! [m] --> [cm]
       dlz(k)  = dlz(k)  * 1.0d+2            !! [m] --> [cm]
    end do
    do k = skst, sken+1
       lzi(k)  = lzi(k)  * 1.0d+2            !! [m] --> [cm]
       dlzi(k) = dlzi(k) * 1.0d+2            !! [m] --> [cm]
    end do

    dswnit   = dswnit   * 1.0d+4             !! [m^2/s] --> [cm^2/s]
    dswdic   = dswdic   * 1.0d+4             !! [m^2/s] --> [cm^2/s]
    dswalk   = dswalk   * 1.0d+4             !! [m^2/s] --> [cm^2/s]
    dswdo2   = dswdo2   * 1.0d+4             !! [m^2/s] --> [cm^2/s]
    dswsil   = dswsil   * 1.0d+4             !! [m^2/s] --> [cm^2/s]
    if ( calc13 .and. sedtdm == 16 ) then
    dsw13dic = dsw13dic * 1.0d+4             !! [m^2/s] --> [cm^2/s]
    end if
    if ( calc14 .and. sedtdm == 16 ) then
    dsw14dic = dsw14dic * 1.0d+4             !! [m^2/s] --> [cm^2/s]
    end if
    dsw    = dsw        * 1.0d+4             !! [m^2/s] --> [cm^2/s]

    mcla   = 1.0d-6 /  rhocla                !! [m^3/g]
    mpom   = 1.0d-9 / ( rhopom / molpom )    !! [m^3/mMol]
    mcar   = 1.0d-9 / ( rhocar / molcar )    !! [m^3/mMol]
    mopa   = 1.0d-9 / ( rhoopa / molopa )    !! [m^3/mMol]
    if ( calc13 .and. sedtdm == 16 ) then
    m13pom = 1.0d-9 / ( rhopom / mol13pom )  !! [m^3/mMol]
    m13car = 1.0d-9 / ( rhocar / mol13car )  !! [m^3/mMol]
    end if
    if ( calc14 .and. sedtdm == 16 ) then
    m14pom = 1.0d-9 / ( rhopom / mol14pom )  !! [m^3/mMol]
    m14car = 1.0d-9 / ( rhocar / mol14car )  !! [m^3/mMol]
    thalf  = thalfyr * 3.1536d+7             !! [yr] --> [s]
    adcay  = log( 2.0d0 ) / thalf / 365.0d0 / 8.64d+4
    end if
    if ( calc13 .or. calc14 ) then
    feoph  = 1.0d0 + eoph * 1.0d-3
    feosh  = 1.0d0 + eosh * 1.0d-3
    end if

    do ij = 1, nxy
      bdep(ij) = bdep(ij) * 1.0d+2           !! [m] --> [cm]
    end do

!   -- Porosity, Diffusivity --

    do k = skst, sken
      phi(k) = phi0 * exp( - lz(k)  / z0 )
    end do

    do k = skst, sken+1
      phii(k) = phi0 * exp( - lzi(k) / z0 )
      dbi(k)  = ( db0 / ( 1.0d0 + exp( 2.0d0 * ( lzi(k) - zb ) ) ) )
    end do

    do k = skst, sken+1
       if (deiarc) then
          deinit(k)   = dswnit   * ( phii(k) ** 3.0d0 )
          deidic(k)   = dswdic   * ( phii(k) ** 3.0d0 )
          deialk(k)   = dswalk   * ( phii(k) ** 3.0d0 )
          deido2(k)   = dswdo2   * ( phii(k) ** 3.0d0 )
          deisil(k)   = dswsil   * ( phii(k) ** 3.0d0 )
          if ( calc13 .and. sedtdm == 16 ) then
          dei13dic(k) = dsw13dic * ( phii(k) ** 3.0d0 )
          end if
          if ( calc14 .and. sedtdm == 16 ) then
          dei14dic(k) = dsw14dic * ( phii(k) ** 3.0d0 )
          end if
       else
          deinit(k)   = ( phii(k) ** ( nopa - 1.0d0 ) ) * dsw
          deidic(k)   = ( phii(k) ** ( nopa - 1.0d0 ) ) * dsw
          deialk(k)   = ( phii(k) ** ( nopa - 1.0d0 ) ) * dsw
          deido2(k)   = ( phii(k) ** ( nopa - 1.0d0 ) ) * dsw
          deisil(k)   = ( phii(k) ** ( nopa - 1.0d0 ) ) * dsw
          if ( calc13 .and. sedtdm == 16 ) then
          dei13dic(k) = ( phii(k) ** ( nopa - 1.0d0 ) ) * dsw
          end if
          if ( calc14 .and. sedtdm == 16 ) then
          dei14dic(k) = ( phii(k) ** ( nopa - 1.0d0 ) ) * dsw
          end if
       end if
    end do

!   -- Set boundary conditions --

    do i = 1, nxydim
      ficlag(i) = bfrcg(i, nbfcla)
      fipomg(i) = bfrcg(i, nbfpom)
      ficarg(i) = bfrcg(i, nbfcar)
      fiopag(i) = bfrcg(i, nbfopa)
      fiageg(i) = bfrcg(i, nbfage)
      bwtemg(i) = bfrcg(i, nbftem)
      bwsalg(i) = bfrcg(i, nbfsal)
      bwnitg(i) = bfrcg(i, nbfnit)
      bwdicg(i) = bfrcg(i, nbfdic)
      bwalkg(i) = bfrcg(i, nbfalk)
      bwsilg(i) = bfrcg(i, nbfsil)
      bwdo2g(i) = bfrcg(i, nbfdo2)
      if ( calc13 .and. sedtdm == 16 ) then
      bw13dicg(i) = bfrcg(i, nbf13dic)
      end if
      if ( calc14 .and. sedtdm == 16 ) then
      bw14dicg(i) = bfrcg(i, nbf14dic)
      end if
    end do

if (use_mpi) then
    do ij = 1, nxy
       ficla(ij) = max( ficlag(myrank*nxy+ij), 0.0d0 )
       fipom(ij) = max( fipomg(myrank*nxy+ij), 0.0d0 )
       ficar(ij) = max( ficarg(myrank*nxy+ij), 0.0d0 )
       fiopa(ij) = max( fiopag(myrank*nxy+ij), 0.0d0 )
       fiage(ij) = 0.0d0
    end do
else
    do ij = 1, nxy
       ficla(ij) = max( ficlag(ij), 0.0d0 )
       fipom(ij) = max( fipomg(ij), 0.0d0 )
       ficar(ij) = max( ficarg(ij), 0.0d0 )
       fiopa(ij) = max( fiopag(ij), 0.0d0 )
       fiage(ij) = 0.0d0
    end do
end if
if (use_mpi) then
    do ij = 1, nxy
       bwtem(ij)   = bwtemg(myrank*nxy+ij)
       bwsal(ij)   = bwsalg(myrank*nxy+ij)
       bwnit(ij)   = max( bwnitg(myrank*nxy+ij), 0.0d0 )
       bwdic(ij)   = max( bwdicg(myrank*nxy+ij), 0.0d0 )
       bwalk(ij)   = max( bwalkg(myrank*nxy+ij), 0.0d0 )
       bwdo2(ij)   = max( bwdo2g(myrank*nxy+ij), 0.0d0 )
       bwsil(ij)   = max( bwsilg(myrank*nxy+ij), 0.0d0 )
    end do
else
    do ij = 1, nxy
       bwtem(ij)   = bwtemg(ij)
       bwsal(ij)   = bwsalg(ij)
       bwnit(ij)   = max( bwnitg(ij), 0.0d0 )
       bwdic(ij)   = max( bwdicg(ij), 0.0d0 )
       bwalk(ij)   = max( bwalkg(ij), 0.0d0 )
       bwdo2(ij)   = max( bwdo2g(ij), 0.0d0 )
       bwsil(ij)   = max( bwsilg(ij), 0.0d0 )
    end do
end if

    if ( calc13 .and. sedtdm == 16 ) then
    do ij = ijtstr, ijtend
       fi13pom(ij) = feoph * pdb  * fipom(ij)
       fi13car(ij) = feosh * pdb  * ficar(ij)
    end do
if (use_mpi) then
    do ij = 1, nxy
       bw13dic(ij) = max( bw13dicg(myrank*nxy+ij), 0.0d0 )
    end do
else
    do ij = 1, nxy
       bw13dic(ij) = max( bw13dicg(ij), 0.0d0 )
    end do
end if
    end if

    if ( calc14 .and. sedtdm == 16 ) then
    do ij = ijtstr, ijtend
       fi14pom(ij) = feoph * dstd * fipom(ij)
       fi14car(ij) = feosh * dstd * ficar(ij)
    end do
if (use_mpi) then
    do ij = 1, nxy
       bw14dic(ij) = max( bw14dicg(myrank*nxy+ij), 0.0d0 )
    end do
else
    do ij = 1, nxy
       bw14dic(ij) = max( bw14dicg(ij), 0.0d0 )
    end do
end if
    end if

    do ij = 1, nxy
       sit(ij) = 0.51d0
       pt(ij)  = 7.70d0
    end do

!   do l = 1, sedtdm
!      do k = 1, sdzdim
!         do j = 1, nydim
!            do i = 1, nxdim
!               ij = i + (j - 1) * nxdim
!               sdtxbuf(i, j, k, l) = sdtx(ij, k, l)
!            end do
!         end do
!      end do
!   end do
!   do l = 1, sedtdm
!      do k = skst, sken
!         do j = jstr, jend
!            do i = 1, istr-1
!               sdtxbuf(i, j, k, l) = sdtxbuf(iend-istr+1+i, j, k, l)
!            end do
!            do i = iend+1, nxdim
!               sdtxbuf(i, j, k, l) = sdtxbuf(istr-1+i-iend, j, k, l)
!            end do
!         end do
!      end do
!   end do
!   do l = 1, sedtdm
!      do k = 1, sdzdim
!         do j = 1, nydim
!            do i = 1, nxdim
!               ij = i + (j - 1) * nxdim
!               sdtx(ij, k, l) = sdtxbuf(i, j, k, l)
!            end do
!         end do
!      end do
!   end do
!   do l = 1, sedtdm
!      do k = 1, skst-1
!         do ij = 1, nxy
!            sdtx(ij, k, l) = sdtx(ij, skst, l)
!         end do
!      end do
!   end do

      call clcend('PRDC_INIT')

!     -- Time Integration --

      write(nfstdo, *) "'--------------------------------------------------'"
      write(nfstdo, *) "'-- DIGDA Prediction:                            --'"
      write(nfstdo, *) "'-- Started time integration                     --'"
      write(nfstdo, *) "'--------------------------------------------------'"

      call clcstr('PRDC_STEP')

      do tt = tstr, tend

         yy = int( sdtsyr * tt )

      if (.not. sdrest) then
         do k = skst, sken
            do ij = ijtstr, ijtend
               scla(ij, k) = rhocla * 1.0d+6
               spom(ij, k) = 0.0d0
               scar(ij, k) = 0.0d0
               sopa(ij, k) = 0.0d0
               sage(ij, k) = 0.0d0
               cnit(ij, k) = bwnit(ij)
               cdic(ij, k) = bwdic(ij)
               calk(ij, k) = bwalk(ij)
               cdo2(ij, k) = bwdo2(ij)
               csil(ij, k) = bwsil(ij)
               if ( calc13 .and. sedtdm == 16 ) then
               s13pom(ij, k) = 0.0d0
               s13car(ij, k) = 0.0d0
               c13dic(ij, k) = bw13dic(ij)
               end if
               if ( calc14 .and. sedtdm == 16 ) then
               s14pom(ij, k) = 0.0d0
               s14car(ij, k) = 0.0d0
               c14dic(ij, k) = bw14dic(ij)
               end if
               pcla(ij, k) = 1.00+2
               ppom(ij, k) = 0.0d0
               pcar(ij, k) = 0.0d0
               popa(ij, k) = 0.0d0
            end do
         end do
         sdrest = .false.
      else
if (use_mpi) then
         do k = skst, sken
            do ij = 1, nxy
               scla(ij, k) = sdtx(myrank*nxy+ij, k,  1)
               spom(ij, k) = sdtx(myrank*nxy+ij, k,  2)
               scar(ij, k) = sdtx(myrank*nxy+ij, k,  3)
               sopa(ij, k) = sdtx(myrank*nxy+ij, k,  4)
               sage(ij, k) = sdtx(myrank*nxy+ij, k,  5)
               cnit(ij, k) = sdtx(myrank*nxy+ij, k,  6)
               cdic(ij, k) = sdtx(myrank*nxy+ij, k,  7)
               calk(ij, k) = sdtx(myrank*nxy+ij, k,  8)
               csil(ij, k) = sdtx(myrank*nxy+ij, k,  9)
               cdo2(ij, k) = sdtx(myrank*nxy+ij, k, 10)
            end do
         end do
         if ( calc13 .and. sedtdm == 16 ) then
         do k = skst, sken
            do ij = 1, nxy
               s13pom(ij, k) = sdtx(myrank*nxy+ij, k, 11)
               s13car(ij, k) = sdtx(myrank*nxy+ij, k, 12)
               c13dic(ij, k) = sdtx(myrank*nxy+ij, k, 13)
            end do
         end do
         end if
         if ( calc14 .and. sedtdm == 16 ) then
         do k = skst, sken
            do ij = 1, nxy
               s14pom(ij, k) = sdtx(myrank*nxy+ij, k, 14)
               s14car(ij, k) = sdtx(myrank*nxy+ij, k, 15)
               c14dic(ij, k) = sdtx(myrank*nxy+ij, k, 16)
            end do
         end do
         end if
else
         do k = skst, sken
            do ij = 1, nxy
               scla(ij, k) = sdtx(ij, k,  1)
               spom(ij, k) = sdtx(ij, k,  2)
               scar(ij, k) = sdtx(ij, k,  3)
               sopa(ij, k) = sdtx(ij, k,  4)
               sage(ij, k) = sdtx(ij, k,  5)
               cnit(ij, k) = sdtx(ij, k,  6)
               cdic(ij, k) = sdtx(ij, k,  7)
               calk(ij, k) = sdtx(ij, k,  8)
               csil(ij, k) = sdtx(ij, k,  9)
               cdo2(ij, k) = sdtx(ij, k, 10)
            end do
         end do
         if ( calc13 .and. sedtdm == 16 ) then
         do k = skst, sken
            do ij = 1, nxy
               s13pom(ij, k) = sdtx(ij, k, 11)
               s13car(ij, k) = sdtx(ij, k, 12)
               c13dic(ij, k) = sdtx(ij, k, 13)
            end do
         end do
         end if
         if ( calc14 .and. sedtdm == 16 ) then
         do k = skst, sken
            do ij = 1, nxy
               s14pom(ij, k) = sdtx(ij, k, 14)
               s14car(ij, k) = sdtx(ij, k, 15)
               c14dic(ij, k) = sdtx(ij, k, 16)
            end do
         end do
         end if
end if
      end if

!     -- Saturated carbonate ion concentration --

      call calc_digda_chemeq( &
     &    bwtem, bwsal, cdic, calk, &
     &    sit, pt, &
     &    cco32 )

      call calc_carbonate_saturation( &
     &    bwtem, bwsal, bdep, &
     &    satco32 )

!     -- Dissolution term --

      call calc_clay_dissolution( &
     &    scla, &
     &    rcla )

      call calc_pom_dissolution( &
     &    spom, cdo2, cp, op, phi, &
     &    rpom )

      call calc_carbonate_dissolution( &
     &    scar, cco32, satco32, phi, &
     &    rcar )

      call calc_silicate_dissolution( &
     &    sopa, csil, pcla, popa, bwtem, &
     &    ropa, satsil )

!     -- Dissolution rate --
!        Clay, POM, CaCO3, SiO2

      do k = skst, sken
         do ij = ijtstr, ijtend
            trcla(ij, k)   = ( 1.0d0 - phi(k) ) * rcla(ij, k) &
     &                     * scla(ij, k)
            trpom(ij, k)   = ( 1.0d0 - phi(k) ) * rpom(ij, k) &
     &                     * spom(ij, k)
            trcar(ij, k)   = ( 1.0d0 - phi(k) ) * rcar(ij, k) &
     &                     * scar(ij, k)
            tropa(ij, k)   = ( 1.0d0 - phi(k) ) * ropa(ij, k) &
     &                     * sopa(ij, k)
         end do
      end do

      if ( calc13 .and. sedtdm == 16 ) then
      do k = skst, sken
         do ij = ijtstr, ijtend
            tr13pom(ij, k) = ( 1.0d0 - phi(k) ) * rpom(ij, k) &
     &                     * s13pom(ij, k)
            tr13car(ij, k) = ( 1.0d0 - phi(k) ) * rcar(ij, k) &
     &                     * s13car(ij, k)
         end do
      end do
      end if

      if ( calc14 .and. sedtdm == 16 ) then
      do k = skst, sken
         do ij = ijtstr, ijtend
            tr14pom(ij, k) = ( 1.0d0 - phi(k) ) * rpom(ij, k) &
     &                     * s14pom(ij, k)
            tr14car(ij, k) = ( 1.0d0 - phi(k) ) * rcar(ij, k) &
     &                     * s14car(ij, k)
         end do
      end do
      end if

!     -- Sedimentation rate --

      do ij = ijtstr, ijtend
         w(ij, skst) = ( ( ( ficla(ij) * mcla ) &
     &                   + ( fipom(ij) * mpom ) &
     &                   + ( ficar(ij) * mcar ) &
     &                   + ( fiopa(ij) * mopa ) ) &
     &               / ( 1.0d0 - phii(skst) ) )
         do k = skst, sken
            w(ij, k+1) = ( ( 1.0d0 - phii(k) ) / ( 1.0d0 - phii(k+1) ) &
     &                   * w(ij, k) ) &
     &                 - ( dlz(k) * ( ( trpom(ij, k) * mpom ) &
     &                              + ( trcar(ij, k) * mcar ) &
     &                              + ( tropa(ij, k) * mopa ) ) &
     &                   / ( 1.0d0 - phii(k+1) ) )
         end do
      end do

!     -- Particles (t --> t+dt) --
!        Clay, POM, CaCO3, SiO2

      call calc_step_solid_clay( &
     &    ficla, trcla, w, phi, phii, dbi, scla, rhocla, &
     &    sncla, bfcla )

      call calc_step_solid( &
     &    fipom, trpom, w, phi, phii, dbi, spom, &
     &    snpom, bfpom )

      call calc_step_solid( &
     &    ficar, trcar, w, phi, phii, dbi, scar, &
     &    sncar, bfcar )

      call calc_step_solid( &
     &    fiopa, tropa, w, phi, phii, dbi, sopa, &
     &    snopa, bfopa )

      call calc_step_solid_age( &
     &    fiage, w, phi, phii, dbi, sage, &
     &    snage )

      if ( calc13 .and. sedtdm == 16 ) then
      call calc_step_solid( &
     &    fi13pom, tr13pom, w, phi, phii, dbi, s13pom, &
     &    sn13pom, bf13pom )

      call calc_step_solid( &
     &    fi13car, tr13car, w, phi, phii, dbi, s13car, &
     &    sn13car, bf13car )
      end if

      if ( calc14 .and. sedtdm == 16 ) then
      call calc_step_solid( &
     &    fi14pom, tr14pom, w, phi, phii, dbi, s14pom, &
     &    sn14pom, bf14pom )

      call calc_step_solid( &
     &    fi14car, tr14car, w, phi, phii, dbi, s14car, &
     &    sn14car, bf14car )
      do k = skst, sken
         do ij = ijtstr, ijtend
            sn14pom(ij, k) = sn14pom(ij, k) &
     &                     - adcay * tssec * s14pom(ij, k)
            sn14car(ij, k) = sn14car(ij, k) &
     &                     - adcay * tssec * s14car(ij, k)
         end do
      end do
      end if

!     -- Dissolution rate --
!        Phosphate, DIC, Alkalinity, Oxygen, Silicate

      do k = skst, sken
         do ij = ijtstr, ijtend
            trnit(ij, k)   = ( ( pp / cp ) * trpom(ij, k) )
            trdic(ij, k)   = ( trpom(ij, k) + trcar(ij, k) )
            tralk(ij, k)   = ( - ( np / cp ) * trpom(ij, k) &
     &                         + ( 2.0d0 * trcar(ij, k) ) )
            trdo2(ij, k)   = - ( ( op / cp ) * trpom(ij, k) )
            if ( cdo2(ij, k) <= 0.0d0 ) then
               trdo2(ij, k) = 0.0d0
            end if
            trsil(ij, k)   = tropa(ij, k)
         end do
      end do

      if ( calc13 .and. sedtdm == 16 ) then
      do k = skst, sken
         do ij = ijtstr, ijtend
            tr13dic(ij, k) = ( tr13pom(ij, k) + tr13car(ij, k) )
         end do
      end do
      end if

      if ( calc14 .and. sedtdm == 16 ) then
      do k = skst, sken
         do ij = ijtstr, ijtend
            tr14dic(ij, k) = ( tr14pom(ij, k) + tr14car(ij, k) )
         end do
      end do
      end if

!     -- Dissolved material (t --> t+dt) --
!        Phosphate, DIC, Alkalinity, Oxygen, Silicate

      call calc_step_liquid( &
     &    bwnit, trnit, phi, phii, deinit, cnit, &
     &    cnnit )

      call calc_step_liquid( &
     &    bwdic, trdic, phi, phii, deidic, cdic, &
     &    cndic )

      call calc_step_liquid( &
     &    bwalk, tralk, phi, phii, deialk, calk, &
     &    cnalk )

      call calc_step_liquid_o2( &
     &    bwdo2, trdo2, phi, phii, deido2, cdo2, &
     &    cndo2 )

      call calc_step_liquid( &
     &    bwsil, trsil, phi, phii, deisil, csil, &
     &    cnsil )

      if ( calc13 .and. sedtdm == 16 ) then
      call calc_step_liquid( &
     &    bw13dic, tr13dic, phi, phii, dei13dic, c13dic, &
     &    cn13dic )
      end if

      if ( calc14 .and. sedtdm == 16 ) then
      call calc_step_liquid( &
     &    bw14dic, tr14dic, phi, phii, dei14dic, c14dic, &
     &    cn14dic )
      end if
      do k = skst, sken
         do ij = ijtstr, ijtend
            cn14dic(ij, k) = cn14dic(ij, k) &
     &                     - adcay * tssec * c14dic(ij, k)
         end do
      end do

!     -- Particles weight percent --

      do k = skst, sken
         do ij = ijtstr, ijtend
            wt(ij, k) = sncla(ij, k) &
     &                + ( snpom(ij, k) * molpom * 1.0d-3 ) &
     &                + ( sncar(ij, k) * molcar * 1.0d-3 ) &
     &                + ( snopa(ij, k) * molopa * 1.0d-3 )
            if ( wt(ij, k) == 0.0d0 ) then
               pcla(ij, k) = 0.0d0
               ppom(ij, k) = 0.0d0
               pcar(ij, k) = 0.0d0
               popa(ij, k) = 0.0d0
            else
               pcla(ij, k) = &
     &            ( sncla(ij, k) / wt(ij, k) ) &
     &            * 1.0d+2
               ppom(ij, k) = &
     &            ( snpom(ij, k) * molpom * 1.0d-3 / wt(ij, k) ) &
     &            * 1.0d+2
               pcar(ij, k) = &
     &            ( sncar(ij, k) * molcar * 1.0d-3 / wt(ij, k) ) &
     &            * 1.0d+2
               popa(ij, k) = &
     &            ( snopa(ij, k) * molopa * 1.0d-3 / wt(ij, k) ) &
     &            * 1.0d+2
            end if
         end do
      end do

!     -- Isotopic fractionation --

      if ( calc13 .and. sedtdm == 16 ) then
      do k = skst, sken
         do ij = ijtstr, ijtend
            if ( cndic(ij, k) /= 0.0d0 ) then
               r13 = cn13dic(ij, k) / max( cndic(ij, k), eps )
               cdc13(ij, k) = &
     &            ( 1.0d0 / 0.0112372d0 &
     &              * r13 / ( 1.0d0 - r13 ) - 1.0d0 ) * 1.0d+3
            else
               cdc13(ij, k) = 0.0d0
            end if
         end do
      end do
      end if

      if ( calc14 .and. sedtdm == 16 ) then
      do k = skst, sken
         do ij = ijtstr, ijtend
            if ( cndic(ij, k) /= 0.0d0 ) then
               r14 = cn14dic(ij, k) / max( cndic(ij, k), eps )
               cdc14(ij, k) = &
     &            ( 1.0d0 / 1.1892d-12 &
     &              * r14 / ( 1.0d0 - r14 ) - 1.0d0 ) * 1.0d+3
            else
               cdc14(ij, k) = 0.0d0
            end if
         end do
      end do
      do k = skst, sken
         do ij = ijtstr, ijtend
            dltc13 = cdc13(ij, k) * 1.0d-3  !! [permil]
            dltc14 = cdc14(ij, k) * 1.0d-3  !! [permil]
            cdc14(ij, k) = &
     &         ( dltc14 - 2.0d0 * ( dltc13 + 0.025d0 ) &
     &                          * ( 1.0d0 + dltc14 ) &
     &         ) * 1.0d+3   !! [permil]
         end do
      end do
      end if

!     -- Sea-bottom flux --

      do ij = ijtstr, ijtend
         fonit(ij) = - ( phii(skst) * deinit(skst) &
     &      * ( bwnit(ij) - cnnit(ij, skst) ) / dlzi(skst) )
         fodic(ij) = - ( phii(skst) * deidic(skst) &
     &      * ( bwdic(ij) - cndic(ij, skst) ) / dlzi(skst) )
         foalk(ij) = - ( phii(skst) * deialk(skst) &
     &      * ( bwalk(ij) - cnalk(ij, skst) ) / dlzi(skst) )
         fodo2(ij) = - ( phii(skst) * deido2(skst) &
     &      * ( bwdo2(ij) - cndo2(ij, skst) ) / dlzi(skst) )
         fosil(ij) = - ( phii(skst) * deisil(skst) &
     &      * ( bwsil(ij) - cnsil(ij, skst) ) / dlzi(skst) )
         finit(ij)  = fipom(ij) * ( pp / cp )
         fidic(ij)  = fipom(ij) + ficar(ij)
         fialk(ij)  = - fipom(ij) * ( np / cp ) + 2.0d0 * ficar(ij)
         fido2(ij)  = fipom(ij) * ( op / cp )
         fisil(ij)  = fiopa(ij)
         if ( calc13 .and. sedtdm == 16 ) then
         fo13dic(ij) = - ( phii(skst) * dei13dic(skst) &
     &      * ( bw13dic(ij) - cn13dic(ij, skst) ) / dlzi(skst) )
         fi13dic(ij)  = fi13pom(ij) + fi13car(ij)
         end if
         if ( calc14 .and. sedtdm == 16 ) then
         fo14dic(ij) = - ( phii(skst) * dei14dic(skst) &
     &      * ( bw14dic(ij) - cn14dic(ij, skst) ) / dlzi(skst) )
         fi14dic(ij)  = fi14pom(ij) + fi14car(ij)
         end if
      end do

      do ij = ijtstr, ijtend
         fopom(ij) = cp * fonit(ij)
         focar(ij) = 0.5d0 * ( + np * fonit(ij) + foalk(ij) )
!         fopom(ij) = ( cp / ( np + 2.0d0 * cp ) ) &
!     &      * ( 2.0d0 * fodic(ij) - foalk(ij) )
!         focar(ij) = ( cp / ( np + 2.0d0 * cp ) ) &
!     &      * ( ( np / cp ) * fodic(ij) + foalk(ij) )
         foopa(ij) = fosil(ij)
      end do

      do ij = ijtstr, ijtend
         fiopom(ij) = fipom(ij) - fopom(ij)
         fiocar(ij) = ficar(ij) - focar(ij)
         fioopa(ij) = fiopa(ij) - foopa(ij)
      end do

!     -- Copy the array for the next time step --

if (use_mpi) then
      do k = skst, sken
         do ij = 1, nxy
            sdtx(myrank*nxy+ij, k,  1) = sncla(ij, k)
            sdtx(myrank*nxy+ij, k,  2) = snpom(ij, k)
            sdtx(myrank*nxy+ij, k,  3) = sncar(ij, k)
            sdtx(myrank*nxy+ij, k,  4) = snopa(ij, k)
            sdtx(myrank*nxy+ij, k,  5) = snage(ij, k)
            sdtx(myrank*nxy+ij, k,  6) = cnnit(ij, k)
            sdtx(myrank*nxy+ij, k,  7) = cndic(ij, k)
            sdtx(myrank*nxy+ij, k,  8) = cnalk(ij, k)
            sdtx(myrank*nxy+ij, k,  9) = cnsil(ij, k)
            sdtx(myrank*nxy+ij, k, 10) = cndo2(ij, k)
         end do
      end do
      if ( calc13 .and. sedtdm == 16 ) then
      do k = skst, sken
         do ij = 1, nxy
            sdtx(myrank*nxy+ij, k, 11) = sn13pom(ij, k)
            sdtx(myrank*nxy+ij, k, 12) = sn13car(ij, k)
            sdtx(myrank*nxy+ij, k, 13) = cn13dic(ij, k)
         end do
      end do
      end if
      if ( calc14 .and. sedtdm == 16 ) then
      do k = skst, sken
         do ij = 1, nxy
            sdtx(myrank*nxy+ij, k, 14) = sn14pom(ij, k)
            sdtx(myrank*nxy+ij, k, 15) = sn14car(ij, k)
            sdtx(myrank*nxy+ij, k, 16) = cn14dic(ij, k)
         end do
      end do
      end if
else
      do k = skst, sken
         do ij = 1, nxy
            sdtx(ij, k,  1) = sncla(ij, k)
            sdtx(ij, k,  2) = snpom(ij, k)
            sdtx(ij, k,  3) = sncar(ij, k)
            sdtx(ij, k,  4) = snopa(ij, k)
            sdtx(ij, k,  5) = snage(ij, k)
            sdtx(ij, k,  6) = cnnit(ij, k)
            sdtx(ij, k,  7) = cndic(ij, k)
            sdtx(ij, k,  8) = cnalk(ij, k)
            sdtx(ij, k,  9) = cnsil(ij, k)
            sdtx(ij, k, 10) = cndo2(ij, k)
         end do
      end do
      if ( calc13 .and. sedtdm == 16 ) then
      do k = skst, sken
         do ij = 1, nxy
            sdtx(ij, k, 11) = sn13pom(ij, k)
            sdtx(ij, k, 12) = sn13car(ij, k)
            sdtx(ij, k, 13) = cn13dic(ij, k)
         end do
      end do
      end if
      if ( calc14 .and. sedtdm == 16 ) then
      do k = skst, sken
         do ij = 1, nxy
            sdtx(ij, k, 14) = sn14pom(ij, k)
            sdtx(ij, k, 15) = sn14car(ij, k)
            sdtx(ij, k, 16) = cn14dic(ij, k)
         end do
      end do
      end if
end if

      if (hstavg) then

!     -- Update the time mean for variables --

        call update_time_mean(mean_data(1),  scla)
        call update_time_mean(mean_data(2),  spom)
        call update_time_mean(mean_data(3),  scar)
        call update_time_mean(mean_data(4),  sopa)
        call update_time_mean(mean_data(5),  sage)
        call update_time_mean(mean_data(6),  cnit)
        call update_time_mean(mean_data(7),  cdic)
        call update_time_mean(mean_data(8),  calk)
        call update_time_mean(mean_data(9),  csil)
        call update_time_mean(mean_data(10), cdo2)

        allocate(varavgdt(nxy,1))
        do ij = 1, nxy
           varavgdt(ij, 1) = ficla(ij)
        end do
        call update_time_mean(mean_data(11),  varavgdt)
        if (allocated(varavgdt)) deallocate(varavgdt)

        allocate(varavgdt(nxy,1))
        do ij = 1, nxy
           varavgdt(ij, 1) = fipom(ij)
        end do
        call update_time_mean(mean_data(12),  varavgdt)
        if (allocated(varavgdt)) deallocate(varavgdt)

        allocate(varavgdt(nxy,1))
        do ij = 1, nxy
           varavgdt(ij, 1) = ficar(ij)
        end do
        call update_time_mean(mean_data(13),  varavgdt)
        if (allocated(varavgdt)) deallocate(varavgdt)

        allocate(varavgdt(nxy,1))
        do ij = 1, nxy
           varavgdt(ij, 1) = fiopa(ij)
        end do
        call update_time_mean(mean_data(14),  varavgdt)
        if (allocated(varavgdt)) deallocate(varavgdt)

        allocate(varavgdt(nxy,1))
        do ij = 1, nxy
           varavgdt(ij, 1) = focla(ij)
        end do
        call update_time_mean(mean_data(15),  varavgdt)
        if (allocated(varavgdt)) deallocate(varavgdt)

        allocate(varavgdt(nxy,1))
        do ij = 1, nxy
           varavgdt(ij, 1) = fopom(ij)
        end do
        call update_time_mean(mean_data(16),  varavgdt)
        if (allocated(varavgdt)) deallocate(varavgdt)

        allocate(varavgdt(nxy,1))
        do ij = 1, nxy
           varavgdt(ij, 1) = focar(ij)
        end do
        call update_time_mean(mean_data(17),  varavgdt)
        if (allocated(varavgdt)) deallocate(varavgdt)

        allocate(varavgdt(nxy,1))
        do ij = 1, nxy
           varavgdt(ij, 1) = foopa(ij)
        end do
        call update_time_mean(mean_data(18),  varavgdt)
        if (allocated(varavgdt)) deallocate(varavgdt)

        allocate(varavgdt(nxy,1))
        do ij = 1, nxy
           varavgdt(ij, 1) = fiocla(ij)
        end do
        call update_time_mean(mean_data(19),  varavgdt)
        if (allocated(varavgdt)) deallocate(varavgdt)

        allocate(varavgdt(nxy,1))
        do ij = 1, nxy
           varavgdt(ij, 1) = fiopom(ij)
        end do
        call update_time_mean(mean_data(20),  varavgdt)
        if (allocated(varavgdt)) deallocate(varavgdt)

        allocate(varavgdt(nxy,1))
        do ij = 1, nxy
           varavgdt(ij, 1) = fiocar(ij)
        end do
        call update_time_mean(mean_data(21),  varavgdt)
        if (allocated(varavgdt)) deallocate(varavgdt)

        allocate(varavgdt(nxy,1))
        do ij = 1, nxy
           varavgdt(ij, 1) = fioopa(ij)
        end do
        call update_time_mean(mean_data(22),  varavgdt)
        if (allocated(varavgdt)) deallocate(varavgdt)

        allocate(varavgdt(nxy,1))
        do ij = 1, nxy
           varavgdt(ij, 1) = bfcla(ij)
        end do
        call update_time_mean(mean_data(23),  varavgdt)
        if (allocated(varavgdt)) deallocate(varavgdt)

        allocate(varavgdt(nxy,1))
        do ij = 1, nxy
           varavgdt(ij, 1) = bfpom(ij)
        end do
        call update_time_mean(mean_data(24),  varavgdt)
        if (allocated(varavgdt)) deallocate(varavgdt)

        allocate(varavgdt(nxy,1))
        do ij = 1, nxy
           varavgdt(ij, 1) = bfcar(ij)
        end do
        call update_time_mean(mean_data(25),  varavgdt)
        if (allocated(varavgdt)) deallocate(varavgdt)

        allocate(varavgdt(nxy,1))
        do ij = 1, nxy
           varavgdt(ij, 1) = bfopa(ij)
        end do
        call update_time_mean(mean_data(26),  varavgdt)
        if (allocated(varavgdt)) deallocate(varavgdt)

        call update_time_mean(mean_data(27), cco32)
        call update_time_mean(mean_data(28), satsil)

        allocate(varavgdt(nxy,1))
        do ij = 1, nxy
           varavgdt(ij, 1) = satco32(ij)
        end do
        call update_time_mean(mean_data(29),  varavgdt)
        if (allocated(varavgdt)) deallocate(varavgdt)

        call update_time_mean(mean_data(30), w)

        if ( calc13 .and. sedtdm == 16 ) then

        call update_time_mean(mean_data(31), s13pom)
        call update_time_mean(mean_data(32), s13car)
        call update_time_mean(mean_data(33), c13dic)
        call update_time_mean(mean_data(34), cdc13)

        allocate(varavgdt(nxy,1))
        do ij = 1, nxy
           varavgdt(ij, 1) = fi13pom(ij)
        end do
        call update_time_mean(mean_data(35),  varavgdt)
        if (allocated(varavgdt)) deallocate(varavgdt)

        allocate(varavgdt(nxy,1))
        do ij = 1, nxy
           varavgdt(ij, 1) = fi13car(ij)
        end do
        call update_time_mean(mean_data(36),  varavgdt)
        if (allocated(varavgdt)) deallocate(varavgdt)

        allocate(varavgdt(nxy,1))
        do ij = 1, nxy
           varavgdt(ij, 1) = fo13pom(ij)
        end do
        call update_time_mean(mean_data(37),  varavgdt)
        if (allocated(varavgdt)) deallocate(varavgdt)

        allocate(varavgdt(nxy,1))
        do ij = 1, nxy
           varavgdt(ij, 1) = fo13car(ij)
        end do
        call update_time_mean(mean_data(38),  varavgdt)
        if (allocated(varavgdt)) deallocate(varavgdt)

        allocate(varavgdt(nxy,1))
        do ij = 1, nxy
           varavgdt(ij, 1) = bf13pom(ij)
        end do
        call update_time_mean(mean_data(39),  varavgdt)
        if (allocated(varavgdt)) deallocate(varavgdt)

        allocate(varavgdt(nxy,1))
        do ij = 1, nxy
           varavgdt(ij, 1) = bf13car(ij)
        end do
        call update_time_mean(mean_data(40),  varavgdt)
        if (allocated(varavgdt)) deallocate(varavgdt)

        end if

        if ( calc14 .and. sedtdm == 16 ) then

        call update_time_mean(mean_data(31), s14pom)
        call update_time_mean(mean_data(32), s14car)
        call update_time_mean(mean_data(33), c14dic)
        call update_time_mean(mean_data(34), cdc14)

        allocate(varavgdt(nxy,1))
        do ij = 1, nxy
           varavgdt(ij, 1) = fi14pom(ij)
        end do
        call update_time_mean(mean_data(35),  varavgdt)
        if (allocated(varavgdt)) deallocate(varavgdt)

        allocate(varavgdt(nxy,1))
        do ij = 1, nxy
           varavgdt(ij, 1) = fi14car(ij)
        end do
        call update_time_mean(mean_data(36),  varavgdt)
        if (allocated(varavgdt)) deallocate(varavgdt)

        allocate(varavgdt(nxy,1))
        do ij = 1, nxy
           varavgdt(ij, 1) = fo14pom(ij)
        end do
        call update_time_mean(mean_data(37),  varavgdt)
        if (allocated(varavgdt)) deallocate(varavgdt)

        allocate(varavgdt(nxy,1))
        do ij = 1, nxy
           varavgdt(ij, 1) = fo14car(ij)
        end do
        call update_time_mean(mean_data(38),  varavgdt)
        if (allocated(varavgdt)) deallocate(varavgdt)

        allocate(varavgdt(nxy,1))
        do ij = 1, nxy
           varavgdt(ij, 1) = bf14pom(ij)
        end do
        call update_time_mean(mean_data(39),  varavgdt)
        if (allocated(varavgdt)) deallocate(varavgdt)

        allocate(varavgdt(nxy,1))
        do ij = 1, nxy
           varavgdt(ij, 1) = bf14car(ij)
        end do
        call update_time_mean(mean_data(40),  varavgdt)
        if (allocated(varavgdt)) deallocate(varavgdt)

        end if

      end if ! -- End of hstavg

!     -- Time step --
      if ( wotstp ) then
         write(nfsots, *), 'Time Step :', tt, 'step', yy, 'year'
         call flush(nfsots)
      end if

!     -- Write the output --
!     if ( wohist ) then
        if ( tt == tend .or. mod(tt, thst) == 0 ) then
!     else
!       if ( tt == tend ) then
!     end if

      call clcstr('PRDC_HIST')

if (use_mpi) then

     !if (sdofirst) then
      allocate(nrcvct2d(0:nprocs-1))
      allocate(ndispt2d(0:nprocs-1))
      allocate(nrcvct3d(0:nprocs-1))
      allocate(ndispt3d(0:nprocs-1))
      allocate(nrcvcw3d(0:nprocs-1))
      allocate(ndispw3d(0:nprocs-1))

      allocate(vardt2d(nxydim+nprocs))
      allocate(vardt3d(nxyzdm+nprocs))
      allocate(vardw3d(nxyzdm+nxydim+nprocs))

      nsndct2d = nxy
      do irank = 0, nprocs-1
         nrcvct2d(irank) = nxy
         ndispt2d(irank) = nxy * irank
      end do
      nsndct3d = nxyz
      do irank = 0, nprocs-1
         nrcvct3d(irank) = nxyz
         ndispt3d(irank) = nxyz * irank
      end do
      nsndcw3d = ( nxyz + nxy )
      do irank = 0, nprocs-1
         nrcvcw3d(irank) = ( nxyz + nxy )
         ndispw3d(irank) = ( nxyz + nxy ) * irank
      end do

     !sdofirst = .false.
     !end if

if (hstavg) then

!     -- Compute the time mean for variables --

      allocate(varavgdt(nxy,sdz))
      call get_time_mean(mean_data(1),  varavgdt)
      call mpi_gatherv( varavgdt, nsndct3d, mpi_real8, &
     &                  vardt3d,  nrcvct3d, ndispt3d, mpi_real8, &
     &                  iroot,    mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            ijkg = ndispt3d(irank)
            do k = skst, sken
               do ij = jjstr(irank), jjend(irank)
                  ijkg = ijkg + 1
                  sclag(ij, k) = vardt3d(ijkg)
               end do
            end do
         end do
      end if
      call reset_time_mean(mean_data(1))
      if (allocated(varavgdt)) deallocate(varavgdt)

      allocate(varavgdt(nxy,sdz))
      call get_time_mean(mean_data(2),  varavgdt)
      call mpi_gatherv( varavgdt, nsndct3d, mpi_real8, &
     &                  vardt3d,  nrcvct3d, ndispt3d, mpi_real8, &
     &                  iroot,    mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            ijkg = ndispt3d(irank)
            do k = skst, sken
               do ij = jjstr(irank), jjend(irank)
                  ijkg = ijkg + 1
                  spomg(ij, k) = vardt3d(ijkg)
               end do
            end do
         end do
      end if
      call reset_time_mean(mean_data(2))
      if (allocated(varavgdt)) deallocate(varavgdt)

      allocate(varavgdt(nxy,sdz))
      call get_time_mean(mean_data(3),  varavgdt)
      call mpi_gatherv( varavgdt, nsndct3d, mpi_real8, &
     &                  vardt3d,  nrcvct3d, ndispt3d, mpi_real8, &
     &                  iroot,    mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            ijkg = ndispt3d(irank)
            do k = skst, sken
               do ij = jjstr(irank), jjend(irank)
                  ijkg = ijkg + 1
                  scarg(ij, k) = vardt3d(ijkg)
               end do
            end do
         end do
      end if
      call reset_time_mean(mean_data(3))
      if (allocated(varavgdt)) deallocate(varavgdt)

      allocate(varavgdt(nxy,sdz))
      call get_time_mean(mean_data(4),  varavgdt)
      call mpi_gatherv( varavgdt, nsndct3d, mpi_real8, &
     &                  vardt3d,  nrcvct3d, ndispt3d, mpi_real8, &
     &                  iroot,    mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            ijkg = ndispt3d(irank)
            do k = skst, sken
               do ij = jjstr(irank), jjend(irank)
                  ijkg = ijkg + 1
                  sopag(ij, k) = vardt3d(ijkg)
               end do
            end do
         end do
      end if
      call reset_time_mean(mean_data(4))
      if (allocated(varavgdt)) deallocate(varavgdt)

      allocate(varavgdt(nxy,sdz))
      call get_time_mean(mean_data(5),  varavgdt)
      call mpi_gatherv( varavgdt, nsndct3d, mpi_real8, &
     &                  vardt3d,  nrcvct3d, ndispt3d, mpi_real8, &
     &                  iroot,    mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            ijkg = ndispt3d(irank)
            do k = skst, sken
               do ij = jjstr(irank), jjend(irank)
                  ijkg = ijkg + 1
                  sageg(ij, k) = vardt3d(ijkg)
               end do
            end do
         end do
      end if
      call reset_time_mean(mean_data(5))
      if (allocated(varavgdt)) deallocate(varavgdt)

      allocate(varavgdt(nxy,sdz))
      call get_time_mean(mean_data(6),  varavgdt)
      call mpi_gatherv( varavgdt, nsndct3d, mpi_real8, &
     &                  vardt3d,  nrcvct3d, ndispt3d, mpi_real8, &
     &                  iroot,    mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            ijkg = ndispt3d(irank)
            do k = skst, sken
               do ij = jjstr(irank), jjend(irank)
                  ijkg = ijkg + 1
                  cnitg(ij, k) = vardt3d(ijkg)
               end do
            end do
         end do
      end if
      call reset_time_mean(mean_data(6))
      if (allocated(varavgdt)) deallocate(varavgdt)

      allocate(varavgdt(nxy,sdz))
      call get_time_mean(mean_data(7),  varavgdt)
      call mpi_gatherv( varavgdt, nsndct3d, mpi_real8, &
     &                  vardt3d,  nrcvct3d, ndispt3d, mpi_real8, &
     &                  iroot,    mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            ijkg = ndispt3d(irank)
            do k = skst, sken
               do ij = jjstr(irank), jjend(irank)
                  ijkg = ijkg + 1
                  cdicg(ij, k) = vardt3d(ijkg)
               end do
            end do
         end do
      end if
      call reset_time_mean(mean_data(7))
      if (allocated(varavgdt)) deallocate(varavgdt)

      allocate(varavgdt(nxy,sdz))
      call get_time_mean(mean_data(8),  varavgdt)
      call mpi_gatherv( varavgdt, nsndct3d, mpi_real8, &
     &                  vardt3d,  nrcvct3d, ndispt3d, mpi_real8, &
     &                  iroot,    mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            ijkg = ndispt3d(irank)
            do k = skst, sken
               do ij = jjstr(irank), jjend(irank)
                  ijkg = ijkg + 1
                  calkg(ij, k) = vardt3d(ijkg)
               end do
            end do
         end do
      end if
      call reset_time_mean(mean_data(8))
      if (allocated(varavgdt)) deallocate(varavgdt)

      allocate(varavgdt(nxy,sdz))
      call get_time_mean(mean_data(9),  varavgdt)
      call mpi_gatherv( varavgdt, nsndct3d, mpi_real8, &
     &                  vardt3d,  nrcvct3d, ndispt3d, mpi_real8, &
     &                  iroot,    mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            ijkg = ndispt3d(irank)
            do k = skst, sken
               do ij = jjstr(irank), jjend(irank)
                  ijkg = ijkg + 1
                  csilg(ij, k) = vardt3d(ijkg)
               end do
            end do
         end do
      end if
      call reset_time_mean(mean_data(9))
      if (allocated(varavgdt)) deallocate(varavgdt)

      allocate(varavgdt(nxy,sdz))
      call get_time_mean(mean_data(10),  varavgdt)
      call mpi_gatherv( varavgdt, nsndct3d, mpi_real8, &
     &                  vardt3d,  nrcvct3d, ndispt3d, mpi_real8, &
     &                  iroot,    mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            ijkg = ndispt3d(irank)
            do k = skst, sken
               do ij = jjstr(irank), jjend(irank)
                  ijkg = ijkg + 1
                  cdo2g(ij, k) = vardt3d(ijkg)
               end do
            end do
         end do
      end if
      call reset_time_mean(mean_data(10))
      if (allocated(varavgdt)) deallocate(varavgdt)

      allocate(varavgdt(nxy,sdz))
      call get_time_mean(mean_data(27),  varavgdt)
      call mpi_gatherv( varavgdt, nsndct3d, mpi_real8, &
     &                  vardt3d,  nrcvct3d, ndispt3d, mpi_real8, &
     &                  iroot,    mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            ijkg = ndispt3d(irank)
            do k = skst, sken
               do ij = jjstr(irank), jjend(irank)
                  ijkg = ijkg + 1
                  cco32g(ij, k) = vardt3d(ijkg)
               end do
            end do
         end do
      end if
      call reset_time_mean(mean_data(27))
      if (allocated(varavgdt)) deallocate(varavgdt)

      allocate(varavgdt(nxy,sdz))
      call get_time_mean(mean_data(28),  varavgdt)
      call mpi_gatherv( varavgdt, nsndct3d, mpi_real8, &
     &                  vardt3d,  nrcvct3d, ndispt3d, mpi_real8, &
     &                  iroot,    mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            ijkg = ndispt3d(irank)
            do k = skst, sken
               do ij = jjstr(irank), jjend(irank)
                  ijkg = ijkg + 1
                  satsilg(ij, k) = vardt3d(ijkg)
               end do
            end do
         end do
      end if
      call reset_time_mean(mean_data(28))
      if (allocated(varavgdt)) deallocate(varavgdt)

      allocate(varavgdt(nxy,1))
      call get_time_mean(mean_data(29),  varavgdt)
      call mpi_gatherv( varavgdt, nsndct2d, mpi_real8, &
     &                  vardt2d,  nrcvct2d, ndispt3d, mpi_real8, &
     &                  iroot,    mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            do ij = jjstr(irank), jjend(irank)
               satco32g(ij) = vardt2d(ij)
            end do
         end do
      end if
      call reset_time_mean(mean_data(29))
      if (allocated(varavgdt)) deallocate(varavgdt)

      allocate(varavgdt(nxy,sdz+1))
      call get_time_mean(mean_data(30),  varavgdt)
      call mpi_gatherv( varavgdt, nsndcw3d, mpi_real8, &
     &                  vardw3d,  nrcvcw3d, ndispw3d, mpi_real8, &
     &                  iroot,    mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            ijkg = ndispw3d(irank)
            do k = skst, sken+1
               do ij = jjstr(irank), jjend(irank)
                  ijkg = ijkg + 1
                  wg(ij, k) = vardw3d(ijkg)
               end do
            end do
         end do
      end if
      call reset_time_mean(mean_data(30))
      if (allocated(varavgdt)) deallocate(varavgdt)

      allocate(varavgdt(nxy,1))
      call get_time_mean(mean_data(11),  varavgdt)
      call mpi_gatherv( ficla,   nsndct2d, mpi_real8, &
     &                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            do ij = jjstr(irank), jjend(irank)
               ficlag(ij) = vardt2d(ij)
            end do
         end do
      end if
      call reset_time_mean(mean_data(11))
      if (allocated(varavgdt)) deallocate(varavgdt)

      allocate(varavgdt(nxy,1))
      call get_time_mean(mean_data(12),  varavgdt)
      call mpi_gatherv( fipom,   nsndct2d, mpi_real8, &
     &                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            do ij = jjstr(irank), jjend(irank)
               fipomg(ij) = vardt2d(ij)
            end do
         end do
      end if
      call reset_time_mean(mean_data(12))
      if (allocated(varavgdt)) deallocate(varavgdt)

      allocate(varavgdt(nxy,1))
      call get_time_mean(mean_data(13),  varavgdt)
      call mpi_gatherv( ficar,   nsndct2d, mpi_real8, &
     &                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            do ij = jjstr(irank), jjend(irank)
               ficarg(ij) = vardt2d(ij)
            end do
         end do
      end if
      call reset_time_mean(mean_data(13))
      if (allocated(varavgdt)) deallocate(varavgdt)

      allocate(varavgdt(nxy,1))
      call get_time_mean(mean_data(14),  varavgdt)
      call mpi_gatherv( fiopa,   nsndct2d, mpi_real8, &
     &                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            do ij = jjstr(irank), jjend(irank)
               fiopag(ij) = vardt2d(ij)
            end do
         end do
      end if
      call reset_time_mean(mean_data(14))
      if (allocated(varavgdt)) deallocate(varavgdt)

 !!!! allocate(varavgdt(nxy,1))
 !!!! call get_time_mean(mean_data(15),  varavgdt)
 !!!! call mpi_gatherv( focla,   nsndct2d, mpi_real8, &
 !!!!&                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
 !!!!&                  iroot,   mpi_comm_world, ierr )
 !!!! if ( myrank .eq. iroot ) then
 !!!!    do irank = 0, nprocs-1
 !!!!       do ij = jjstr(irank), jjend(irank)
 !!!!          foclag(ij) = vardt2d(ij)
 !!!!       end do
 !!!!    end do
 !!!! end if
 !!!! call reset_time_mean(mean_data(15))
 !!!! if (allocated(varavgdt)) deallocate(varavgdt)

      allocate(varavgdt(nxy,1))
      call get_time_mean(mean_data(16),  varavgdt)
      call mpi_gatherv( fopom,   nsndct2d, mpi_real8, &
     &                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            do ij = jjstr(irank), jjend(irank)
               fopomg(ij) = vardt2d(ij)
            end do
         end do
      end if
      call reset_time_mean(mean_data(16))
      if (allocated(varavgdt)) deallocate(varavgdt)

      allocate(varavgdt(nxy,1))
      call get_time_mean(mean_data(17),  varavgdt)
      call mpi_gatherv( focar,   nsndct2d, mpi_real8, &
     &                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            do ij = jjstr(irank), jjend(irank)
               focarg(ij) = vardt2d(ij)
            end do
         end do
      end if
      call reset_time_mean(mean_data(17))
      if (allocated(varavgdt)) deallocate(varavgdt)

      allocate(varavgdt(nxy,1))
      call get_time_mean(mean_data(18),  varavgdt)
      call mpi_gatherv( foopa,   nsndct2d, mpi_real8, &
     &                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            do ij = jjstr(irank), jjend(irank)
               foopag(ij) = vardt2d(ij)
            end do
         end do
      end if
      call reset_time_mean(mean_data(18))
      if (allocated(varavgdt)) deallocate(varavgdt)

 !!!! allocate(varavgdt(nxy,1))
 !!!! call get_time_mean(mean_data(19),  varavgdt)
 !!!! call mpi_gatherv( fiocla,  nsndct2d, mpi_real8, &
 !!!!&                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
 !!!!&                  iroot,   mpi_comm_world, ierr )
 !!!! if ( myrank .eq. iroot ) then
 !!!!    do irank = 0, nprocs-1
 !!!!       do ij = jjstr(irank), jjend(irank)
 !!!!          fioclag(ij) = vardt2d(ij)
 !!!!       end do
 !!!!    end do
 !!!! end if
 !!!! call reset_time_mean(mean_data(19))
 !!!! if (allocated(varavgdt)) deallocate(varavgdt)

      allocate(varavgdt(nxy,1))
      call get_time_mean(mean_data(20),  varavgdt)
      call mpi_gatherv( fiopom,  nsndct2d, mpi_real8, &
     &                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            do ij = jjstr(irank), jjend(irank)
               fiopomg(ij) = vardt2d(ij)
            end do
         end do
      end if
      call reset_time_mean(mean_data(20))
      if (allocated(varavgdt)) deallocate(varavgdt)

      allocate(varavgdt(nxy,1))
      call get_time_mean(mean_data(21),  varavgdt)
      call mpi_gatherv( fiocar,  nsndct2d, mpi_real8, &
     &                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            do ij = jjstr(irank), jjend(irank)
               fiocarg(ij) = vardt2d(ij)
            end do
         end do
      end if
      call reset_time_mean(mean_data(21))
      if (allocated(varavgdt)) deallocate(varavgdt)

      allocate(varavgdt(nxy,1))
      call get_time_mean(mean_data(22),  varavgdt)
      call mpi_gatherv( fioopa,  nsndct2d, mpi_real8, &
     &                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            do ij = jjstr(irank), jjend(irank)
               fioopag(ij) = vardt2d(ij)
            end do
         end do
      end if
      call reset_time_mean(mean_data(22))
      if (allocated(varavgdt)) deallocate(varavgdt)

      allocate(varavgdt(nxy,1))
      call get_time_mean(mean_data(23),  varavgdt)
      call mpi_gatherv( bfcla,   nsndct2d, mpi_real8, &
     &                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            do ij = jjstr(irank), jjend(irank)
               bfclag(ij) = vardt2d(ij)
            end do
         end do
      end if
      call reset_time_mean(mean_data(23))
      if (allocated(varavgdt)) deallocate(varavgdt)

      allocate(varavgdt(nxy,1))
      call get_time_mean(mean_data(24),  varavgdt)
      call mpi_gatherv( bfpom,   nsndct2d, mpi_real8, &
     &                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            do ij = jjstr(irank), jjend(irank)
               bfpomg(ij) = vardt2d(ij)
            end do
         end do
      end if
      call reset_time_mean(mean_data(24))
      if (allocated(varavgdt)) deallocate(varavgdt)

      allocate(varavgdt(nxy,1))
      call get_time_mean(mean_data(25),  varavgdt)
      call mpi_gatherv( bfcar,   nsndct2d, mpi_real8, &
     &                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            do ij = jjstr(irank), jjend(irank)
               bfcarg(ij) = vardt2d(ij)
            end do
         end do
      end if
      call reset_time_mean(mean_data(25))
      if (allocated(varavgdt)) deallocate(varavgdt)

      allocate(varavgdt(nxy,1))
      call get_time_mean(mean_data(26),  varavgdt)
      call mpi_gatherv( bfopa,   nsndct2d, mpi_real8, &
     &                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            do ij = jjstr(irank), jjend(irank)
               bfopag(ij) = vardt2d(ij)
            end do
         end do
      end if
      call reset_time_mean(mean_data(26))
      if (allocated(varavgdt)) deallocate(varavgdt)

      if ( calc13 .and. sedtdm == 16 ) then

      allocate(varavgdt(nxy,sdz))
      call get_time_mean(mean_data(31),  varavgdt)
      call mpi_gatherv( s13pom,  nsndct3d, mpi_real8, &
     &                  vardt3d, nrcvct3d, ndispt3d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            ijkg = ndispt3d(irank)
            do k = skst, sken
               do ij = jjstr(irank), jjend(irank)
                  ijkg = ijkg + 1
                  s13pomg(ij, k) = vardt3d(ijkg)
               end do
            end do
         end do
      end if
      call reset_time_mean(mean_data(31))
      if (allocated(varavgdt)) deallocate(varavgdt)

      allocate(varavgdt(nxy,sdz))
      call get_time_mean(mean_data(32),  varavgdt)
      call mpi_gatherv( s13car,  nsndct3d, mpi_real8, &
     &                  vardt3d, nrcvct3d, ndispt3d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            ijkg = ndispt3d(irank)
            do k = skst, sken
               do ij = jjstr(irank), jjend(irank)
                  ijkg = ijkg + 1
                  s13carg(ij, k) = vardt3d(ijkg)
               end do
            end do
         end do
      end if
      call reset_time_mean(mean_data(32))
      if (allocated(varavgdt)) deallocate(varavgdt)

      allocate(varavgdt(nxy,sdz))
      call get_time_mean(mean_data(33),  varavgdt)
      call mpi_gatherv( c13dic,  nsndct3d, mpi_real8, &
     &                  vardt3d, nrcvct3d, ndispt3d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            ijkg = ndispt3d(irank)
            do k = skst, sken
               do ij = jjstr(irank), jjend(irank)
                  ijkg = ijkg + 1
                  c13dicg(ij, k) = vardt3d(ijkg)
               end do
            end do
         end do
      end if
      call reset_time_mean(mean_data(33))
      if (allocated(varavgdt)) deallocate(varavgdt)

      allocate(varavgdt(nxy,sdz))
      call get_time_mean(mean_data(34),  varavgdt)
      call mpi_gatherv( cdc13,   nsndct3d, mpi_real8, &
     &                  vardt3d, nrcvct3d, ndispt3d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            ijkg = ndispt3d(irank)
            do k = skst, sken
               do ij = jjstr(irank), jjend(irank)
                  ijkg = ijkg + 1
                  cdc13g(ij, k) = vardt3d(ijkg)
               end do
            end do
         end do
      end if
      call reset_time_mean(mean_data(34))
      if (allocated(varavgdt)) deallocate(varavgdt)

 !!!! allocate(varavgdt(nxy,1))
 !!!! call get_time_mean(mean_data(35),  varavgdt)
 !!!! call mpi_gatherv( fi13pom, nsndct2d, mpi_real8, &
 !!!!&                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
 !!!!&                  iroot,   mpi_comm_world, ierr )
 !!!! if ( myrank .eq. iroot ) then
 !!!!    do irank = 0, nprocs-1
 !!!!       do ij = jjstr(irank), jjend(irank)
 !!!!          fi13pomg(ij) = vardt2d(ij)
 !!!!       end do
 !!!!    end do
 !!!! end if
 !!!! call reset_time_mean(mean_data(35))
 !!!! if (allocated(varavgdt)) deallocate(varavgdt)

 !!!! allocate(varavgdt(nxy,1))
 !!!! call get_time_mean(mean_data(36),  varavgdt)
 !!!! call mpi_gatherv( fi13car, nsndct2d, mpi_real8, &
 !!!!&                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
 !!!!&                  iroot,   mpi_comm_world, ierr )
 !!!! if ( myrank .eq. iroot ) then
 !!!!    do irank = 0, nprocs-1
 !!!!       do ij = jjstr(irank), jjend(irank)
 !!!!          fi13carg(ij) = vardt2d(ij)
 !!!!       end do
 !!!!    end do
 !!!! end if
 !!!! call reset_time_mean(mean_data(36))
 !!!! if (allocated(varavgdt)) deallocate(varavgdt)

 !!!! allocate(varavgdt(nxy,1))
 !!!! call get_time_mean(mean_data(37),  varavgdt)
 !!!! call mpi_gatherv( fo13pom, nsndct2d, mpi_real8, &
 !!!!&                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
 !!!!&                  iroot,   mpi_comm_world, ierr )
 !!!! if ( myrank .eq. iroot ) then
 !!!!    do irank = 0, nprocs-1
 !!!!       do ij = jjstr(irank), jjend(irank)
 !!!!          fo13pomg(ij) = vardt2d(ij)
 !!!!       end do
 !!!!    end do
 !!!! end if
 !!!! call reset_time_mean(mean_data(37))
 !!!! if (allocated(varavgdt)) deallocate(varavgdt)

 !!!! allocate(varavgdt(nxy,1))
 !!!! call get_time_mean(mean_data(38),  varavgdt)
 !!!! call mpi_gatherv( fo13car, nsndct2d, mpi_real8, &
 !!!!&                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
 !!!!&                  iroot,   mpi_comm_world, ierr )
 !!!! if ( myrank .eq. iroot ) then
 !!!!    do irank = 0, nprocs-1
 !!!!       do ij = jjstr(irank), jjend(irank)
 !!!!          fo13carg(ij) = vardt2d(ij)
 !!!!       end do
 !!!!    end do
 !!!! end if
 !!!! call reset_time_mean(mean_data(38))
 !!!! if (allocated(varavgdt)) deallocate(varavgdt)

 !!!! allocate(varavgdt(nxy,1))
 !!!! call get_time_mean(mean_data(39),  varavgdt)
 !!!! call mpi_gatherv( bf13pom, nsndct2d, mpi_real8, &
 !!!!&                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
 !!!!&                  iroot,   mpi_comm_world, ierr )
 !!!! if ( myrank .eq. iroot ) then
 !!!!    do irank = 0, nprocs-1
 !!!!       do ij = jjstr(irank), jjend(irank)
 !!!!          bf13pomg(ij) = vardt2d(ij)
 !!!!       end do
 !!!!    end do
 !!!! end if
 !!!! call reset_time_mean(mean_data(39))
 !!!! if (allocated(varavgdt)) deallocate(varavgdt)

 !!!! allocate(varavgdt(nxy,1))
 !!!! call get_time_mean(mean_data(40),  varavgdt)
 !!!! call mpi_gatherv( bf13car, nsndct2d, mpi_real8, &
 !!!!&                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
 !!!!&                  iroot,   mpi_comm_world, ierr )
 !!!! if ( myrank .eq. iroot ) then
 !!!!    do irank = 0, nprocs-1
 !!!!       do ij = jjstr(irank), jjend(irank)
 !!!!          bf13carg(ij) = vardt2d(ij)
 !!!!       end do
 !!!!    end do
 !!!! end if
 !!!! call reset_time_mean(mean_data(40))
 !!!! if (allocated(varavgdt)) deallocate(varavgdt)

      end if

      if ( calc14 .and. sedtdm == 16 ) then

      allocate(varavgdt(nxy,sdz))
      call get_time_mean(mean_data(41),  varavgdt)
      call mpi_gatherv( s14pom,  nsndct3d, mpi_real8, &
     &                  vardt3d, nrcvct3d, ndispt3d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            ijkg = ndispt3d(irank)
            do k = skst, sken
               do ij = jjstr(irank), jjend(irank)
                  ijkg = ijkg + 1
                  s14pomg(ij, k) = vardt3d(ijkg)
               end do
            end do
         end do
      end if
      call reset_time_mean(mean_data(41))
      if (allocated(varavgdt)) deallocate(varavgdt)

      allocate(varavgdt(nxy,sdz))
      call get_time_mean(mean_data(42),  varavgdt)
      call mpi_gatherv( s14car,  nsndct3d, mpi_real8, &
     &                  vardt3d, nrcvct3d, ndispt3d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            ijkg = ndispt3d(irank)
            do k = skst, sken
               do ij = jjstr(irank), jjend(irank)
                  ijkg = ijkg + 1
                  s14carg(ij, k) = vardt3d(ijkg)
               end do
            end do
         end do
      end if
      call reset_time_mean(mean_data(42))
      if (allocated(varavgdt)) deallocate(varavgdt)

      allocate(varavgdt(nxy,sdz))
      call get_time_mean(mean_data(43),  varavgdt)
      call mpi_gatherv( c14dic,  nsndct3d, mpi_real8, &
     &                  vardt3d, nrcvct3d, ndispt3d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            ijkg = ndispt3d(irank)
            do k = skst, sken
               do ij = jjstr(irank), jjend(irank)
                  ijkg = ijkg + 1
                  c14dicg(ij, k) = vardt3d(ijkg)
               end do
            end do
         end do
      end if
      call reset_time_mean(mean_data(43))
      if (allocated(varavgdt)) deallocate(varavgdt)

      allocate(varavgdt(nxy,sdz))
      call get_time_mean(mean_data(44),  varavgdt)
      call mpi_gatherv( cdc14,   nsndct3d, mpi_real8, &
     &                  vardt3d, nrcvct3d, ndispt3d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            ijkg = ndispt3d(irank)
            do k = skst, sken
               do ij = jjstr(irank), jjend(irank)
                  ijkg = ijkg + 1
                  cdc14g(ij, k) = vardt3d(ijkg)
               end do
            end do
         end do
      end if
      call reset_time_mean(mean_data(44))
      if (allocated(varavgdt)) deallocate(varavgdt)

 !!!! allocate(varavgdt(nxy,1))
 !!!! call get_time_mean(mean_data(45),  varavgdt)
 !!!! call mpi_gatherv( fi14pom, nsndct2d, mpi_real8, &
 !!!!&                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
 !!!!&                  iroot,   mpi_comm_world, ierr )
 !!!! if ( myrank .eq. iroot ) then
 !!!!    do irank = 0, nprocs-1
 !!!!       do ij = jjstr(irank), jjend(irank)
 !!!!          fi14pomg(ij) = vardt2d(ij)
 !!!!       end do
 !!!!    end do
 !!!! end if
 !!!! call reset_time_mean(mean_data(45))
 !!!! if (allocated(varavgdt)) deallocate(varavgdt)

 !!!! allocate(varavgdt(nxy,1))
 !!!! call get_time_mean(mean_data(46),  varavgdt)
 !!!! call mpi_gatherv( fi14car, nsndct2d, mpi_real8, &
 !!!!&                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
 !!!!&                  iroot,   mpi_comm_world, ierr )
 !!!! if ( myrank .eq. iroot ) then
 !!!!    do irank = 0, nprocs-1
 !!!!       do ij = jjstr(irank), jjend(irank)
 !!!!          fi14carg(ij) = vardt2d(ij)
 !!!!       end do
 !!!!    end do
 !!!! end if
 !!!! call reset_time_mean(mean_data(46))
 !!!! if (allocated(varavgdt)) deallocate(varavgdt)

 !!!! allocate(varavgdt(nxy,1))
 !!!! call get_time_mean(mean_data(47),  varavgdt)
 !!!! call mpi_gatherv( fo14pom, nsndct2d, mpi_real8, &
 !!!!&                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
 !!!!&                  iroot,   mpi_comm_world, ierr )
 !!!! if ( myrank .eq. iroot ) then
 !!!!    do irank = 0, nprocs-1
 !!!!       do ij = jjstr(irank), jjend(irank)
 !!!!          fo14pomg(ij) = vardt2d(ij)
 !!!!       end do
 !!!!    end do
 !!!! end if
 !!!! call reset_time_mean(mean_data(47))
 !!!! if (allocated(varavgdt)) deallocate(varavgdt)

 !!!! allocate(varavgdt(nxy,1))
 !!!! call get_time_mean(mean_data(48),  varavgdt)
 !!!! call mpi_gatherv( fo14car, nsndct2d, mpi_real8, &
 !!!!&                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
 !!!!&                  iroot,   mpi_comm_world, ierr )
 !!!! if ( myrank .eq. iroot ) then
 !!!!    do irank = 0, nprocs-1
 !!!!       do ij = jjstr(irank), jjend(irank)
 !!!!          fo14carg(ij) = vardt2d(ij)
 !!!!       end do
 !!!!    end do
 !!!! end if
 !!!! call reset_time_mean(mean_data(48))
 !!!! if (allocated(varavgdt)) deallocate(varavgdt)

 !!!! allocate(varavgdt(nxy,1))
 !!!! call get_time_mean(mean_data(49),  varavgdt)
 !!!! call mpi_gatherv( bf14pom, nsndct2d, mpi_real8, &
 !!!!&                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
 !!!!&                  iroot,   mpi_comm_world, ierr )
 !!!! if ( myrank .eq. iroot ) then
 !!!!    do irank = 0, nprocs-1
 !!!!       do ij = jjstr(irank), jjend(irank)
 !!!!          bf14pomg(ij) = vardt2d(ij)
 !!!!       end do
 !!!!    end do
 !!!! end if
 !!!! call reset_time_mean(mean_data(49))
 !!!! if (allocated(varavgdt)) deallocate(varavgdt)

 !!!! allocate(varavgdt(nxy,1))
 !!!! call get_time_mean(mean_data(50),  varavgdt)
 !!!! call mpi_gatherv( bf14car, nsndct2d, mpi_real8, &
 !!!!&                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
 !!!!&                  iroot,   mpi_comm_world, ierr )
 !!!! if ( myrank .eq. iroot ) then
 !!!!    do irank = 0, nprocs-1
 !!!!       do ij = jjstr(irank), jjend(irank)
 !!!!          bf14carg(ij) = vardt2d(ij)
 !!!!       end do
 !!!!    end do
 !!!! end if
 !!!! call reset_time_mean(mean_data(50))
 !!!! if (allocated(varavgdt)) deallocate(varavgdt)

      end if

 !!!! if ( myrank .eq. iroot ) then
 !!!!    do k = skst, sken
 !!!!       do ij = 1, nxydim
 !!!!          wtg(ij, k)   = sclag(ij, k) &
 !!!!&                      + ( spomg(ij, k) * molpom * 1.0d-3 ) &
 !!!!&                      + ( scarg(ij, k) * molcar * 1.0d-3 ) &
 !!!!&                      + ( sopag(ij, k) * molopa * 1.0d-3 )
 !!!!          pclag(ij, k) = sclag(ij, k) &
 !!!!&                      / wtg(ij, k) * 1.0d+2
 !!!!          ppomg(ij, k) = spomg(ij, k) * molpom * 1.0d-3 &
 !!!!&                      / wtg(ij, k) * 1.0d+2
 !!!!          pcarg(ij, k) = scarg(ij, k) * molcar * 1.0d-3 &
 !!!!&                      / wtg(ij, k) * 1.0d+2
 !!!!          popag(ij, k) = sopag(ij, k) * molopa * 1.0d-3 &
 !!!!&                      / wtg(ij, k) * 1.0d+2
 !!!!       end do
 !!!!    end do

 !!!! end if


else

      call mpi_gatherv( scla,    nsndct3d, mpi_real8, &
     &                  vardt3d, nrcvct3d, ndispt3d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            ijkg = ndispt3d(irank)
            do k = skst, sken
               do ij = jjstr(irank), jjend(irank)
                  ijkg = ijkg + 1
                  sclag(ij, k) = vardt3d(ijkg)
               end do
            end do
         end do
      end if

      call mpi_gatherv( spom,    nsndct3d, mpi_real8, &
     &                  vardt3d, nrcvct3d, ndispt3d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            ijkg = ndispt3d(irank)
            do k = skst, sken
               do ij = jjstr(irank), jjend(irank)
                  ijkg = ijkg + 1
                  spomg(ij, k) = vardt3d(ijkg)
               end do
            end do
         end do
      end if

      call mpi_gatherv( scar,    nsndct3d, mpi_real8, &
     &                  vardt3d, nrcvct3d, ndispt3d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            ijkg = ndispt3d(irank)
            do k = skst, sken
               do ij = jjstr(irank), jjend(irank)
                  ijkg = ijkg + 1
                  scarg(ij, k) = vardt3d(ijkg)
               end do
            end do
         end do
      end if

      call mpi_gatherv( sopa,    nsndct3d, mpi_real8, &
     &                  vardt3d, nrcvct3d, ndispt3d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            ijkg = ndispt3d(irank)
            do k = skst, sken
               do ij = jjstr(irank), jjend(irank)
                  ijkg = ijkg + 1
                  sopag(ij, k) = vardt3d(ijkg)
               end do
            end do
         end do
      end if

      call mpi_gatherv( sage,    nsndct3d, mpi_real8, &
     &                  vardt3d, nrcvct3d, ndispt3d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            ijkg = ndispt3d(irank)
            do k = skst, sken
               do ij = jjstr(irank), jjend(irank)
                  ijkg = ijkg + 1
                  sageg(ij, k) = vardt3d(ijkg)
               end do
            end do
         end do
      end if

      call mpi_gatherv( cnit,    nsndct3d, mpi_real8, &
     &                  vardt3d, nrcvct3d, ndispt3d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            ijkg = ndispt3d(irank)
            do k = skst, sken
               do ij = jjstr(irank), jjend(irank)
                  ijkg = ijkg + 1
                  cnitg(ij, k) = vardt3d(ijkg)
               end do
            end do
         end do
      end if

      call mpi_gatherv( cdic,    nsndct3d, mpi_real8, &
     &                  vardt3d, nrcvct3d, ndispt3d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            ijkg = ndispt3d(irank)
            do k = skst, sken
               do ij = jjstr(irank), jjend(irank)
                  ijkg = ijkg + 1
                  cdicg(ij, k) = vardt3d(ijkg)
               end do
            end do
         end do
      end if

      call mpi_gatherv( calk,    nsndct3d, mpi_real8, &
     &                  vardt3d, nrcvct3d, ndispt3d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            ijkg = ndispt3d(irank)
            do k = skst, sken
               do ij = jjstr(irank), jjend(irank)
                  ijkg = ijkg + 1
                  calkg(ij, k) = vardt3d(ijkg)
               end do
            end do
         end do
      end if

      call mpi_gatherv( csil,    nsndct3d, mpi_real8, &
     &                  vardt3d, nrcvct3d, ndispt3d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            ijkg = ndispt3d(irank)
            do k = skst, sken
               do ij = jjstr(irank), jjend(irank)
                  ijkg = ijkg + 1
                  csilg(ij, k) = vardt3d(ijkg)
               end do
            end do
         end do
      end if

      call mpi_gatherv( cdo2,    nsndct3d, mpi_real8, &
     &                  vardt3d, nrcvct3d, ndispt3d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            ijkg = ndispt3d(irank)
            do k = skst, sken
               do ij = jjstr(irank), jjend(irank)
                  ijkg = ijkg + 1
                  cdo2g(ij, k) = vardt3d(ijkg)
               end do
            end do
         end do
      end if

      call mpi_gatherv( cco32,   nsndct3d, mpi_real8, &
     &                  vardt3d, nrcvct3d, ndispt3d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            ijkg = ndispt3d(irank)
            do k = skst, sken
               do ij = jjstr(irank), jjend(irank)
                  ijkg = ijkg + 1
                  cco32g(ij, k) = vardt3d(ijkg)
               end do
            end do
         end do
      end if

      call mpi_gatherv( satsil,  nsndct3d, mpi_real8, &
     &                  vardt3d, nrcvct3d, ndispt3d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            ijkg = ndispt3d(irank)
            do k = skst, sken
               do ij = jjstr(irank), jjend(irank)
                  ijkg = ijkg + 1
                  satsilg(ij, k) = vardt3d(ijkg)
               end do
            end do
         end do
      end if

      call mpi_gatherv( satco32, nsndct2d, mpi_real8, &
     &                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            do ij = jjstr(irank), jjend(irank)
               satco32g(ij) = vardt2d(ij)
            end do
         end do
      end if

      call mpi_gatherv( w,       nsndcw3d, mpi_real8, &
     &                  vardw3d, nrcvcw3d, ndispw3d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            ijkg = ndispw3d(irank)
            do k = skst, sken+1
               do ij = jjstr(irank), jjend(irank)
                  ijkg = ijkg + 1
                  wg(ij, k) = vardw3d(ijkg)
               end do
            end do
         end do
      end if

      call mpi_gatherv( ficla,   nsndct2d, mpi_real8, &
     &                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            do ij = jjstr(irank), jjend(irank)
               ficlag(ij) = vardt2d(ij)
            end do
         end do
      end if

      call mpi_gatherv( fipom,   nsndct2d, mpi_real8, &
     &                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            do ij = jjstr(irank), jjend(irank)
               fipomg(ij) = vardt2d(ij)
            end do
         end do
      end if

      call mpi_gatherv( ficar,   nsndct2d, mpi_real8, &
     &                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            do ij = jjstr(irank), jjend(irank)
               ficarg(ij) = vardt2d(ij)
            end do
         end do
      end if

      call mpi_gatherv( fiopa,   nsndct2d, mpi_real8, &
     &                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            do ij = jjstr(irank), jjend(irank)
               fiopag(ij) = vardt2d(ij)
            end do
         end do
      end if

 !!!! call mpi_gatherv( focla,   nsndct2d, mpi_real8, &
 !!!!&                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
 !!!!&                  iroot,   mpi_comm_world, ierr )
 !!!! if ( myrank .eq. iroot ) then
 !!!!    do irank = 0, nprocs-1
 !!!!       do ij = jjstr(irank), jjend(irank)
 !!!!          foclag(ij) = vardt2d(ij)
 !!!!       end do
 !!!!    end do
 !!!! end if

      call mpi_gatherv( fopom,   nsndct2d, mpi_real8, &
     &                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            do ij = jjstr(irank), jjend(irank)
               fopomg(ij) = vardt2d(ij)
            end do
         end do
      end if

      call mpi_gatherv( focar,   nsndct2d, mpi_real8, &
     &                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            do ij = jjstr(irank), jjend(irank)
               focarg(ij) = vardt2d(ij)
            end do
         end do
      end if

      call mpi_gatherv( foopa,   nsndct2d, mpi_real8, &
     &                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            do ij = jjstr(irank), jjend(irank)
               foopag(ij) = vardt2d(ij)
            end do
         end do
      end if

 !!!! call mpi_gatherv( fiocla,  nsndct2d, mpi_real8, &
 !!!!&                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
 !!!!&                  iroot,   mpi_comm_world, ierr )
 !!!! if ( myrank .eq. iroot ) then
 !!!!    do irank = 0, nprocs-1
 !!!!       do ij = jjstr(irank), jjend(irank)
 !!!!          fioclag(ij) = vardt2d(ij)
 !!!!       end do
 !!!!    end do
 !!!! end if

      call mpi_gatherv( fiopom,  nsndct2d, mpi_real8, &
     &                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            do ij = jjstr(irank), jjend(irank)
               fiopomg(ij) = vardt2d(ij)
            end do
         end do
      end if

      call mpi_gatherv( fiocar,  nsndct2d, mpi_real8, &
     &                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            do ij = jjstr(irank), jjend(irank)
               fiocarg(ij) = vardt2d(ij)
            end do
         end do
      end if

      call mpi_gatherv( fioopa,  nsndct2d, mpi_real8, &
     &                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            do ij = jjstr(irank), jjend(irank)
               fioopag(ij) = vardt2d(ij)
            end do
         end do
      end if

      call mpi_gatherv( bfcla,   nsndct2d, mpi_real8, &
     &                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            do ij = jjstr(irank), jjend(irank)
               bfclag(ij) = vardt2d(ij)
            end do
         end do
      end if

      call mpi_gatherv( bfpom,   nsndct2d, mpi_real8, &
     &                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            do ij = jjstr(irank), jjend(irank)
               bfpomg(ij) = vardt2d(ij)
            end do
         end do
      end if

      call mpi_gatherv( bfcar,   nsndct2d, mpi_real8, &
     &                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            do ij = jjstr(irank), jjend(irank)
               bfcarg(ij) = vardt2d(ij)
            end do
         end do
      end if

      call mpi_gatherv( bfopa,   nsndct2d, mpi_real8, &
     &                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            do ij = jjstr(irank), jjend(irank)
               bfopag(ij) = vardt2d(ij)
            end do
         end do
      end if

      if ( calc13 .and. sedtdm == 16 ) then

      call mpi_gatherv( s13pom,  nsndct3d, mpi_real8, &
     &                  vardt3d, nrcvct3d, ndispt3d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            ijkg = ndispt3d(irank)
            do k = skst, sken
               do ij = jjstr(irank), jjend(irank)
                  ijkg = ijkg + 1
                  s13pomg(ij, k) = vardt3d(ijkg)
               end do
            end do
         end do
      end if

      call mpi_gatherv( s13car,  nsndct3d, mpi_real8, &
     &                  vardt3d, nrcvct3d, ndispt3d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            ijkg = ndispt3d(irank)
            do k = skst, sken
               do ij = jjstr(irank), jjend(irank)
                  ijkg = ijkg + 1
                  s13carg(ij, k) = vardt3d(ijkg)
               end do
            end do
         end do
      end if

      call mpi_gatherv( c13dic,  nsndct3d, mpi_real8, &
     &                  vardt3d, nrcvct3d, ndispt3d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            ijkg = ndispt3d(irank)
            do k = skst, sken
               do ij = jjstr(irank), jjend(irank)
                  ijkg = ijkg + 1
                  c13dicg(ij, k) = vardt3d(ijkg)
               end do
            end do
         end do
      end if

      call mpi_gatherv( cdc13,   nsndct3d, mpi_real8, &
     &                  vardt3d, nrcvct3d, ndispt3d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            ijkg = ndispt3d(irank)
            do k = skst, sken
               do ij = jjstr(irank), jjend(irank)
                  ijkg = ijkg + 1
                  cdc13g(ij, k) = vardt3d(ijkg)
               end do
            end do
         end do
      end if

 !!!! call mpi_gatherv( fi13pom, nsndct2d, mpi_real8, &
 !!!!&                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
 !!!!&                  iroot,   mpi_comm_world, ierr )
 !!!! if ( myrank .eq. iroot ) then
 !!!!    do irank = 0, nprocs-1
 !!!!       do ij = jjstr(irank), jjend(irank)
 !!!!          fi13pomg(ij) = vardt2d(ij)
 !!!!       end do
 !!!!    end do
 !!!! end if

 !!!! call mpi_gatherv( fi13car, nsndct2d, mpi_real8, &
 !!!!&                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
 !!!!&                  iroot,   mpi_comm_world, ierr )
 !!!! if ( myrank .eq. iroot ) then
 !!!!    do irank = 0, nprocs-1
 !!!!       do ij = jjstr(irank), jjend(irank)
 !!!!          fi13carg(ij) = vardt2d(ij)
 !!!!       end do
 !!!!    end do
 !!!! end if

 !!!! call mpi_gatherv( fo13pom, nsndct2d, mpi_real8, &
 !!!!&                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
 !!!!&                  iroot,   mpi_comm_world, ierr )
 !!!! if ( myrank .eq. iroot ) then
 !!!!    do irank = 0, nprocs-1
 !!!!       do ij = jjstr(irank), jjend(irank)
 !!!!          fo13pomg(ij) = vardt2d(ij)
 !!!!       end do
 !!!!    end do
 !!!! end if

 !!!! call mpi_gatherv( fo13car, nsndct2d, mpi_real8, &
 !!!!&                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
 !!!!&                  iroot,   mpi_comm_world, ierr )
 !!!! if ( myrank .eq. iroot ) then
 !!!!    do irank = 0, nprocs-1
 !!!!       do ij = jjstr(irank), jjend(irank)
 !!!!          fo13carg(ij) = vardt2d(ij)
 !!!!       end do
 !!!!    end do
 !!!! end if

 !!!! call mpi_gatherv( bf13pom, nsndct2d, mpi_real8, &
 !!!!&                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
 !!!!&                  iroot,   mpi_comm_world, ierr )
 !!!! if ( myrank .eq. iroot ) then
 !!!!    do irank = 0, nprocs-1
 !!!!       do ij = jjstr(irank), jjend(irank)
 !!!!          bf13pomg(ij) = vardt2d(ij)
 !!!!       end do
 !!!!    end do
 !!!! end if

 !!!! call mpi_gatherv( bf13car, nsndct2d, mpi_real8, &
 !!!!&                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
 !!!!&                  iroot,   mpi_comm_world, ierr )
 !!!! if ( myrank .eq. iroot ) then
 !!!!    do irank = 0, nprocs-1
 !!!!       do ij = jjstr(irank), jjend(irank)
 !!!!          bf13carg(ij) = vardt2d(ij)
 !!!!       end do
 !!!!    end do
 !!!! end if

      end if

      if ( calc14 .and. sedtdm == 16 ) then

      call mpi_gatherv( s14pom,  nsndct3d, mpi_real8, &
     &                  vardt3d, nrcvct3d, ndispt3d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            ijkg = ndispt3d(irank)
            do k = skst, sken
               do ij = jjstr(irank), jjend(irank)
                  ijkg = ijkg + 1
                  s14pomg(ij, k) = vardt3d(ijkg)
               end do
            end do
         end do
      end if

      call mpi_gatherv( s14car,  nsndct3d, mpi_real8, &
     &                  vardt3d, nrcvct3d, ndispt3d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            ijkg = ndispt3d(irank)
            do k = skst, sken
               do ij = jjstr(irank), jjend(irank)
                  ijkg = ijkg + 1
                  s14carg(ij, k) = vardt3d(ijkg)
               end do
            end do
         end do
      end if

      call mpi_gatherv( c14dic,  nsndct3d, mpi_real8, &
     &                  vardt3d, nrcvct3d, ndispt3d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            ijkg = ndispt3d(irank)
            do k = skst, sken
               do ij = jjstr(irank), jjend(irank)
                  ijkg = ijkg + 1
                  c14dicg(ij, k) = vardt3d(ijkg)
               end do
            end do
         end do
      end if

      call mpi_gatherv( cdc14,   nsndct3d, mpi_real8, &
     &                  vardt3d, nrcvct3d, ndispt3d, mpi_real8, &
     &                  iroot,   mpi_comm_world, ierr )
      if ( myrank .eq. iroot ) then
         do irank = 0, nprocs-1
            ijkg = ndispt3d(irank)
            do k = skst, sken
               do ij = jjstr(irank), jjend(irank)
                  ijkg = ijkg + 1
                  cdc14g(ij, k) = vardt3d(ijkg)
               end do
            end do
         end do
      end if

 !!!! call mpi_gatherv( fi14pom, nsndct2d, mpi_real8, &
 !!!!&                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
 !!!!&                  iroot,   mpi_comm_world, ierr )
 !!!! if ( myrank .eq. iroot ) then
 !!!!    do irank = 0, nprocs-1
 !!!!       do ij = jjstr(irank), jjend(irank)
 !!!!          fi14pomg(ij) = vardt2d(ij)
 !!!!       end do
 !!!!    end do
 !!!! end if

 !!!! call mpi_gatherv( fi14car, nsndct2d, mpi_real8, &
 !!!!&                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
 !!!!&                  iroot,   mpi_comm_world, ierr )
 !!!! if ( myrank .eq. iroot ) then
 !!!!    do irank = 0, nprocs-1
 !!!!       do ij = jjstr(irank), jjend(irank)
 !!!!          fi14carg(ij) = vardt2d(ij)
 !!!!       end do
 !!!!    end do
 !!!! end if

 !!!! call mpi_gatherv( fo14pom, nsndct2d, mpi_real8, &
 !!!!&                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
 !!!!&                  iroot,   mpi_comm_world, ierr )
 !!!! if ( myrank .eq. iroot ) then
 !!!!    do irank = 0, nprocs-1
 !!!!       do ij = jjstr(irank), jjend(irank)
 !!!!          fo14pomg(ij) = vardt2d(ij)
 !!!!       end do
 !!!!    end do
 !!!! end if

 !!!! call mpi_gatherv( fo14car, nsndct2d, mpi_real8, &
 !!!!&                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
 !!!!&                  iroot,   mpi_comm_world, ierr )
 !!!! if ( myrank .eq. iroot ) then
 !!!!    do irank = 0, nprocs-1
 !!!!       do ij = jjstr(irank), jjend(irank)
 !!!!          fo14carg(ij) = vardt2d(ij)
 !!!!       end do
 !!!!    end do
 !!!! end if

 !!!! call mpi_gatherv( bf14pom, nsndct2d, mpi_real8, &
 !!!!&                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
 !!!!&                  iroot,   mpi_comm_world, ierr )
 !!!! if ( myrank .eq. iroot ) then
 !!!!    do irank = 0, nprocs-1
 !!!!       do ij = jjstr(irank), jjend(irank)
 !!!!          bf14pomg(ij) = vardt2d(ij)
 !!!!       end do
 !!!!    end do
 !!!! end if

 !!!! call mpi_gatherv( bf14car, nsndct2d, mpi_real8, &
 !!!!&                  vardt2d, nrcvct2d, ndispt2d, mpi_real8, &
 !!!!&                  iroot,   mpi_comm_world, ierr )
 !!!! if ( myrank .eq. iroot ) then
 !!!!    do irank = 0, nprocs-1
 !!!!       do ij = jjstr(irank), jjend(irank)
 !!!!          bf14carg(ij) = vardt2d(ij)
 !!!!       end do
 !!!!    end do
 !!!! end if

      end if

 !!!! if ( myrank .eq. iroot ) then
 !!!!    do k = skst, sken
 !!!!       do ij = 1, nxydim
 !!!!          wtg(ij, k)   = sclag(ij, k) &
 !!!!&                      + ( spomg(ij, k) * molpom * 1.0d-3 ) &
 !!!!&                      + ( scarg(ij, k) * molcar * 1.0d-3 ) &
 !!!!&                      + ( sopag(ij, k) * molopa * 1.0d-3 )
 !!!!          pclag(ij, k) = sclag(ij, k) &
 !!!!&                      / wtg(ij, k) * 1.0d+2
 !!!!          ppomg(ij, k) = spomg(ij, k) * molpom * 1.0d-3 &
 !!!!&                      / wtg(ij, k) * 1.0d+2
 !!!!          pcarg(ij, k) = scarg(ij, k) * molcar * 1.0d-3 &
 !!!!&                      / wtg(ij, k) * 1.0d+2
 !!!!          popag(ij, k) = sopag(ij, k) * molopa * 1.0d-3 &
 !!!!&                      / wtg(ij, k) * 1.0d+2
 !!!!       end do
 !!!!    end do

 !!!! end if

end if ! -- End of hstavg

      call mpi_barrier(mpi_comm_world, ierr)

      if (allocated(nrcvct2d)) deallocate(nrcvct2d)
      if (allocated(ndispt2d)) deallocate(ndispt2d)
      if (allocated(nrcvct3d)) deallocate(nrcvct3d)
      if (allocated(ndispt3d)) deallocate(ndispt3d)
      if (allocated(nrcvcw3d)) deallocate(nrcvcw3d)
      if (allocated(ndispw3d)) deallocate(ndispw3d)

      if (allocated(vardt2d))  deallocate(vardt2d)
      if (allocated(vardt3d))  deallocate(vardt3d)
      if (allocated(vardw3d))  deallocate(vardw3d)

else

         do k = skst, sken
            do ij = 1, nxydim
               sclag(ij, k) = scla(ij, k)
               spomg(ij, k) = spom(ij, k)
               scarg(ij, k) = scar(ij, k)
               sopag(ij, k) = sopa(ij, k)
               sageg(ij, k) = sage(ij, k)
            end do
         end do

!!       do k = skst, sken
!!          do ij = 1, nxydim
!!             pclag(ij, k) = pcla(ij, k)
!!             ppomg(ij, k) = ppom(ij, k)
!!             pcarg(ij, k) = pcar(ij, k)
!!             popag(ij, k) = popa(ij, k)
!!          end do
!!       end do

         do k = skst, sken
            do ij = 1, nxydim
                cnitg(ij, k)   = cnit(ij, k)
                cdicg(ij, k)   = cdic(ij, k)
                calkg(ij, k)   = calk(ij, k)
                cdo2g(ij, k)   = cdo2(ij, k)
                csilg(ij, k)   = csil(ij, k)
                cco32g(ij, k)  = cco32(ij, k)
            end do
         end do

         do k = skst, sken
            do ij = 1, nxydim
                satsilg(ij, k) = satsil(ij, k)
            end do
         end do

         do ij = 1, nxydim
            satco32g(ij) = satco32(ij)
         end do

         do k = skst, sken
            do ij = 1, nxydim
                wg(ij, k) = w(ij, k)
            end do
         end do

         if ( calc13 .and. sedtdm == 16 ) then
         do k = skst, sken
            do ij = 1, nxydim
                s13pomg(ij, k) = s13pom(ij, k)
                s13carg(ij, k) = s13car(ij, k)
                c13dicg(ij, k) = c13dic(ij, k)
                cdc13g(ij, k)  = cdc13(ij, k)
            end do
         end do
         end if

         if ( calc14 .and. sedtdm == 16 ) then
         do k = skst, sken
            do ij = 1, nxydim
                s14pomg(ij, k) = s14pom(ij, k)
                s14carg(ij, k) = s14car(ij, k)
                c14dicg(ij, k) = c14dic(ij, k)
                cdc14g(ij, k)  = cdc14(ij, k)
            end do
         end do
         end if

         do ij = 1, nxydim
            ficlag(ij)     = ficla(ij)
            fipomg(ij)     = fipom(ij)
            ficarg(ij)     = ficar(ij)
            fiopag(ij)     = fiopa(ij)
            foclag(ij)     = focla(ij)
            fopomg(ij)     = fopom(ij)
            focarg(ij)     = focar(ij)
            foopag(ij)     = foopa(ij)
            fioclag(ij)    = fiocla(ij)
            fiopomg(ij)    = fiopom(ij)
            fiocarg(ij)    = fiocar(ij)
            fioopag(ij)    = fioopa(ij)
            bfclag(ij)     = bfcla(ij)
            bfpomg(ij)     = bfpom(ij)
            bfcarg(ij)     = bfcar(ij)
            bfopag(ij)     = bfopa(ij)
         end do

         if ( calc13 .and. sedtdm == 16 ) then
         do ij = 1, nxydim
            fi13pomg(ij)   = fi13pom(ij)
            fi13carg(ij)   = fi13car(ij)
            fo13pomg(ij)   = fo13pom(ij)
            fo13carg(ij)   = fo13car(ij)
            bf13pomg(ij)   = bf13pom(ij)
            bf13carg(ij)   = bf13car(ij)
         end do
         end if

         if ( calc14 .and. sedtdm == 16 ) then
         do ij = 1, nxydim
            fi14pomg(ij)   = fi14pom(ij)
            fi14carg(ij)   = fi14car(ij)
            fo14pomg(ij)   = fo14pom(ij)
            fo14carg(ij)   = fo14car(ij)
            bf14pomg(ij)   = bf14pom(ij)
            bf14carg(ij)   = bf14car(ij)
         end do
         end if

end if ! -- End of use_mpi

if (use_mpi) then
      if ( myrank .eq. iroot ) then

         !-- NetCDF output --
         if (allocated(yy_array)) deallocate(yy_array)
         allocate(yy_array(1))
         yy_array(1) = yy

         allocate(var3d_values(nxdim,nydim,sdzdim,1))
         filename = "scla.nc"
         varname = "Clay"
         varunit = "g m-3"
         do k = 1, sdzdim
           do j = 1, nydim
             do i = 1, nxdim
               ij = i + (j - 1) * nxdim
               var3d_values(i, j, k, 1) = sclag(ij, k)
             end do
           end do
         end do
         call write_netcdf_3dwt(filename, varname, varunit, var3d_values, lon_values, lat_values, lz, yy_array)
         call calc_integral_and_mean_3d(sclag, lon_values, lat_values, lz, integral_value, mean_value)
         write(nfstdo, *) "'--------------------------------------------------'"
         write(nfstdo,'(A, A)') " Filename: ", trim(filename)
         write(nfstdo,'(A, A, A, I6)') " Variable: ", trim(varname), ", Year: ", YY
         write(nfstdo,'(A, E12.5, A)') " Global Averaged: ", mean_value, " g m-3"
         write(nfstdo,'(A, E12.5, A)') " Global Integrated: ", integral_value, " g"
         write(nfstdo, *) "'--------------------------------------------------'"
         if (allocated(var3d_values)) deallocate(var3d_values)

         allocate(var3d_values(nxdim,nydim,sdzdim,1))
         filename = "spom.nc"
         varname = "OrgM"
         varunit = "mmol m-3"
         do k = 1, sdzdim
           do j = 1, nydim
             do i = 1, nxdim
               ij = i + (j - 1) * nxdim
               var3d_values(i, j, k, 1) = spomg(ij, k)
             end do
           end do
         end do
         call write_netcdf_3dwt(filename, varname, varunit, var3d_values, lon_values, lat_values, lz, yy_array)
         call calc_integral_and_mean_3d(spomg, lon_values, lat_values, lz, integral_value, mean_value)
         write(nfstdo, *) "'--------------------------------------------------'"
         write(nfstdo,'(A, A)') " Filename: ", trim(filename)
         write(nfstdo,'(A, A, A, I6)') " Variable: ", trim(varname), ", Year: ", YY
         write(nfstdo,'(A, E12.5, A)') " Global Averaged: ", mean_value, " mmol m-3"
         write(nfstdo,'(A, E12.5, A)') " Global Integrated: ", integral_value * 1.d-3, " mol"
         write(nfstdo, *) "'--------------------------------------------------'"
         if (allocated(var3d_values)) deallocate(var3d_values)

         allocate(var3d_values(nxdim,nydim,sdzdim,1))
         filename = "scar.nc"
         varname = "CaCO3"
         varunit = "mmol m-3"
         do k = 1, sdzdim
           do j = 1, nydim
             do i = 1, nxdim
               ij = i + (j - 1) * nxdim
               var3d_values(i, j, k, 1) = scarg(ij, k)
             end do
           end do
         end do
         call write_netcdf_3dwt(filename, varname, varunit, var3d_values, lon_values, lat_values, lz, yy_array)
         call calc_integral_and_mean_3d(scarg, lon_values, lat_values, lz, integral_value, mean_value)
         write(nfstdo, *) "'--------------------------------------------------'"
         write(nfstdo,'(A, A)') " Filename: ", trim(filename)
         write(nfstdo,'(A, A, A, I6)') " Variable: ", trim(varname), ", Year: ", YY
         write(nfstdo,'(A, E12.5, A)') " Global Averaged: ", mean_value, " mmol m-3"
         write(nfstdo,'(A, E12.5, A)') " Global Integrated: ", integral_value * 1.d-3, " mol"
         write(nfstdo, *) "'--------------------------------------------------'"
         if (allocated(var3d_values)) deallocate(var3d_values)

         allocate(var3d_values(nxdim,nydim,sdzdim,1))
         filename = "sopa.nc"
         varname = "SiO2"
         varunit = "mmol m-3"
         do k = 1, sdzdim
           do j = 1, nydim
             do i = 1, nxdim
               ij = i + (j - 1) * nxdim
               var3d_values(i, j, k, 1) = sopag(ij, k)
             end do
           end do
         end do
         call write_netcdf_3dwt(filename, varname, varunit, var3d_values, lon_values, lat_values, lz, yy_array)
         call calc_integral_and_mean_3d(sopag, lon_values, lat_values, lz, integral_value, mean_value)
         write(nfstdo, *) "'--------------------------------------------------'"
         write(nfstdo,'(A, A)') " Filename: ", trim(filename)
         write(nfstdo,'(A, A, A, I6)') " Variable: ", trim(varname), ", Year: ", YY
         write(nfstdo,'(A, E12.5, A)') " Global Averaged: ", mean_value, " mmol m-3"
         write(nfstdo,'(A, E12.5, A)') " Global Integrated: ", integral_value * 1.d-3, " mol"
         write(nfstdo, *) "'--------------------------------------------------'"
         if (allocated(var3d_values)) deallocate(var3d_values)

         allocate(var3d_values(nxdim,nydim,sdzdim,1))
         filename = "sage.nc"
         varname = "Agetracer"
         varunit = "year"
         do k = 1, sdzdim
           do j = 1, nydim
             do i = 1, nxdim
               ij = i + (j - 1) * nxdim
               var3d_values(i, j, k, 1) = sageg(ij, k)
             end do
           end do
         end do
         call write_netcdf_3dwt(filename, varname, varunit, var3d_values, lon_values, lat_values, lz, yy_array)
         call calc_integral_and_mean_3d(sageg, lon_values, lat_values, lz, integral_value, mean_value)
         write(nfstdo, *) "'--------------------------------------------------'"
         write(nfstdo,'(A, A)') " Filename: ", trim(filename)
         write(nfstdo,'(A, A, A, I6)') " Variable: ", trim(varname), ", Year: ", YY
         write(nfstdo,'(A, E12.5, A)') " Global Averaged: ", mean_value, " yr"
         write(nfstdo, *) "'--------------------------------------------------'"
         if (allocated(var3d_values)) deallocate(var3d_values)

         allocate(var3d_values(nxdim,nydim,sdzdim,1))
         filename = "cnit.nc"
         varname = "Phosphate"
         varunit = "mmol m-3"
         do k = 1, sdzdim
           do j = 1, nydim
             do i = 1, nxdim
               ij = i + (j - 1) * nxdim
               var3d_values(i, j, k, 1) = cnitg(ij, k)
             end do
           end do
         end do
         call write_netcdf_3dwt(filename, varname, varunit, var3d_values, lon_values, lat_values, lz, yy_array)
         call calc_integral_and_mean_3d(cnitg, lon_values, lat_values, lz, integral_value, mean_value)
         write(nfstdo, *) "'--------------------------------------------------'"
         write(nfstdo,'(A, A)') " Filename: ", trim(filename)
         write(nfstdo,'(A, A, A, I6)') " Variable: ", trim(varname), ", Year: ", YY
         write(nfstdo,'(A, E12.5, A)') " Global Averaged: ", mean_value, " mmol m-3"
         write(nfstdo,'(A, E12.5, A)') " Global Integrated: ", integral_value * 1.d-3, " mol"
         write(nfstdo, *) "'--------------------------------------------------'"
         if (allocated(var3d_values)) deallocate(var3d_values)

         allocate(var3d_values(nxdim,nydim,sdzdim,1))
         filename = "cdic.nc"
         varname = "DIC"
         varunit = "mmol m-3"
         do k = 1, sdzdim
           do j = 1, nydim
             do i = 1, nxdim
               ij = i + (j - 1) * nxdim
               var3d_values(i, j, k, 1) = cdicg(ij, k)
             end do
           end do
         end do
         call write_netcdf_3dwt(filename, varname, varunit, var3d_values, lon_values, lat_values, lz, yy_array)
         call calc_integral_and_mean_3d(cdicg, lon_values, lat_values, lz, integral_value, mean_value)
         write(nfstdo, *) "'--------------------------------------------------'"
         write(nfstdo,'(A, A)') " Filename: ", trim(filename)
         write(nfstdo,'(A, A, A, I6)') " Variable: ", trim(varname), ", Year: ", YY
         write(nfstdo,'(A, E12.5, A)') " Global Averaged: ", mean_value, " mmol m-3"
         write(nfstdo,'(A, E12.5, A)') " Global Integrated: ", integral_value * 1.d-3, " mol"
         write(nfstdo, *) "'--------------------------------------------------'"
         if (allocated(var3d_values)) deallocate(var3d_values)

         allocate(var3d_values(nxdim,nydim,sdzdim,1))
         filename = "calk.nc"
         varname = "Alkalinity"
         varunit = "mmol m-3"
         do k = 1, sdzdim
           do j = 1, nydim
             do i = 1, nxdim
               ij = i + (j - 1) * nxdim
               var3d_values(i, j, k, 1) = calkg(ij, k)
             end do
           end do
         end do
         call write_netcdf_3dwt(filename, varname, varunit, var3d_values, lon_values, lat_values, lz, yy_array)
         call calc_integral_and_mean_3d(calkg, lon_values, lat_values, lz, integral_value, mean_value)
         write(nfstdo, *) "'--------------------------------------------------'"
         write(nfstdo,'(A, A)') " Filename: ", trim(filename)
         write(nfstdo,'(A, A, A, I6)') " Variable: ", trim(varname), ", Year: ", YY
         write(nfstdo,'(A, E12.5, A)') " Global Averaged: ", mean_value, " mmol m-3"
         write(nfstdo,'(A, E12.5, A)') " Global Integrated: ", integral_value * 1.d-3, " mol"
         write(nfstdo, *) "'--------------------------------------------------'"
         if (allocated(var3d_values)) deallocate(var3d_values)

         allocate(var3d_values(nxdim,nydim,sdzdim,1))
         filename = "cdo2.nc"
         varname = "Oxygen"
         varunit = "mmol m-3"
         do k = 1, sdzdim
           do j = 1, nydim
             do i = 1, nxdim
               ij = i + (j - 1) * nxdim
               var3d_values(i, j, k, 1) = cdo2g(ij, k)
             end do
           end do
         end do
         call write_netcdf_3dwt(filename, varname, varunit, var3d_values, lon_values, lat_values, lz, yy_array)
         call calc_integral_and_mean_3d(cdo2g, lon_values, lat_values, lz, integral_value, mean_value)
         write(nfstdo, *) "'--------------------------------------------------'"
         write(nfstdo,'(A, A)') " Filename: ", trim(filename)
         write(nfstdo,'(A, A, A, I6)') " Variable: ", trim(varname), ", Year: ", YY
         write(nfstdo,'(A, E12.5, A)') " Global Averaged: ", mean_value, " mmol m-3"
         write(nfstdo,'(A, E12.5, A)') " Global Integrated: ", integral_value * 1.d-3, " mol"
         write(nfstdo, *) "'--------------------------------------------------'"
         if (allocated(var3d_values)) deallocate(var3d_values)

         allocate(var3d_values(nxdim,nydim,sdzdim,1))
         filename = "csil.nc"
         varname = "Silicate"
         varunit = "mmol m-3"
         do k = 1, sdzdim
           do j = 1, nydim
             do i = 1, nxdim
               ij = i + (j - 1) * nxdim
               var3d_values(i, j, k, 1) = csilg(ij, k)
             end do
           end do
         end do
         call write_netcdf_3dwt(filename, varname, varunit, var3d_values, lon_values, lat_values, lz, yy_array)
         call calc_integral_and_mean_3d(csilg, lon_values, lat_values, lz, integral_value, mean_value)
         write(nfstdo, *) "'--------------------------------------------------'"
         write(nfstdo,'(A, A)') " Filename: ", trim(filename)
         write(nfstdo,'(A, A, A, I6)') " Variable: ", trim(varname), ", Year: ", YY
         write(nfstdo,'(A, E12.5, A)') " Global Averaged: ", mean_value, " mmol m-3"
         write(nfstdo,'(A, E12.5, A)') " Global Integrated: ", integral_value * 1.d-3, " mol"
         write(nfstdo, *) "'--------------------------------------------------'"
         if (allocated(var3d_values)) deallocate(var3d_values)

         allocate(var2d_values(nxdim,nydim,1))
         filename = "ficla.nc"
         varname = "Clay_bottom_flux"
         varunit = "g m-3 * cm s-1"
         do j = 1, nydim
           do i = 1, nxdim
             ij = i + (j - 1) * nxdim
             var2d_values(i, j, 1) = ficlag(ij)
           end do
         end do
        !call write_netcdf_2dwt(filename, varname, varunit, var2d_values, lon_values, lat_values, yy_array)
         call calc_integral_and_mean_2d(ficlag, lon_values, lat_values, integral_value, mean_value)
         write(nfstdo, *) "'--------------------------------------------------'"
         write(nfstdo,'(A, A)') " Filename: ", trim(filename)
         write(nfstdo,'(A, A, A, I6)') " Variable: ", trim(varname), ", Year: ", YY
         write(nfstdo,'(A, E12.5, A)') " Global Averaged: ", mean_value, " g m-3 * cm s-1"
         write(nfstdo,'(A, E12.5, A)') " Global Integrated: ", integral_value * 3.1536d+5, " g * yr-1"
         write(nfstdo, *) "'--------------------------------------------------'"
         if (allocated(var2d_values)) deallocate(var2d_values)

         allocate(var2d_values(nxdim,nydim,1))
         filename = "fipom.nc"
         varname = "OrgM_bottom_flux"
         varunit = "mMol m-3 * cm s-1"
         do j = 1, nydim
           do i = 1, nxdim
             ij = i + (j - 1) * nxdim
             var2d_values(i, j, 1) = fipomg(ij)
           end do
         end do
        !call write_netcdf_2dwt(filename, varname, varunit, var2d_values, lon_values, lat_values, yy_array)
         call calc_integral_and_mean_2d(fipomg, lon_values, lat_values, integral_value, mean_value)
         write(nfstdo, *) "'--------------------------------------------------'"
         write(nfstdo,'(A, A)') " Filename: ", trim(filename)
         write(nfstdo,'(A, A, A, I6)') " Variable: ", trim(varname), ", Year: ", YY
         write(nfstdo,'(A, E12.5, A)') " Global Averaged: ", mean_value, " mmol m-3 * cm s-1"
         write(nfstdo,'(A, E12.5, A)') " Global Integrated: ", integral_value * 3.1536d+2, " mol * yr-1"
         write(nfstdo, *) "'--------------------------------------------------'"
         if (allocated(var2d_values)) deallocate(var2d_values)

         allocate(var2d_values(nxdim,nydim,1))
         filename = "ficar.nc"
         varname = "CaCO3_bottom_flux"
         varunit = "mMol m-3 * cm s-1"
         do j = 1, nydim
           do i = 1, nxdim
             ij = i + (j - 1) * nxdim
             var2d_values(i, j, 1) = ficarg(ij)
           end do
         end do
        !call write_netcdf_2dwt(filename, varname, varunit, var2d_values, lon_values, lat_values, yy_array)
         call calc_integral_and_mean_2d(ficarg, lon_values, lat_values, integral_value, mean_value)
         write(nfstdo, *) "'--------------------------------------------------'"
         write(nfstdo,'(A, A)') " Filename: ", trim(filename)
         write(nfstdo,'(A, A, A, I6)') " Variable: ", trim(varname), ", Year: ", YY
         write(nfstdo,'(A, E12.5, A)') " Global Averaged: ", mean_value, " mmol m-3 * cm s-1"
         write(nfstdo,'(A, E12.5, A)') " Global Integrated: ", integral_value * 3.1536d+2, " mol * yr-1"
         write(nfstdo, *) "'--------------------------------------------------'"
         if (allocated(var2d_values)) deallocate(var2d_values)

         allocate(var2d_values(nxdim,nydim,1))
         filename = "fiopa.nc"
         varname = "SiO2_bottom_flux"
         varunit = "mMol m-3 * cm s-1"
         do j = 1, nydim
           do i = 1, nxdim
             ij = i + (j - 1) * nxdim
             var2d_values(i, j, 1) = fiopag(ij)
           end do
         end do
        !call write_netcdf_2dwt(filename, varname, varunit, var2d_values, lon_values, lat_values, yy_array)
         call calc_integral_and_mean_2d(fiopag, lon_values, lat_values, integral_value, mean_value)
         write(nfstdo, *) "'--------------------------------------------------'"
         write(nfstdo,'(A, A)') " Filename: ", trim(filename)
         write(nfstdo,'(A, A, A, I6)') " Variable: ", trim(varname), ", Year: ", YY
         write(nfstdo,'(A, E12.5, A)') " Global Averaged: ", mean_value, " mmol m-3 * cm s-1"
         write(nfstdo,'(A, E12.5, A)') " Global Integrated: ", integral_value * 3.1536d+2, " mol * yr-1"
         write(nfstdo, *) "'--------------------------------------------------'"
         if (allocated(var2d_values)) deallocate(var2d_values)

         allocate(var2d_values(nxdim,nydim,1))
         filename = "focla.nc"
         varname = "Clay_diffusion"
         varunit = "g m-3 * cm s-1"
         do j = 1, nydim
           do i = 1, nxdim
             ij = i + (j - 1) * nxdim
             var2d_values(i, j, 1) = foclag(ij)
           end do
         end do
        !call write_netcdf_2dwt(filename, varname, varunit, var2d_values, lon_values, lat_values, yy_array)
         call calc_integral_and_mean_2d(foclag, lon_values, lat_values, integral_value, mean_value)
         write(nfstdo, *) "'--------------------------------------------------'"
         write(nfstdo,'(A, A)') " Filename: ", trim(filename)
         write(nfstdo,'(A, A, A, I6)') " Variable: ", trim(varname), ", Year: ", YY
         write(nfstdo,'(A, E12.5, A)') " Global Averaged: ", mean_value, " mmol m-3 * cm s-1"
         write(nfstdo,'(A, E12.5, A)') " Global Integrated: ", integral_value * 3.1536d+2, " mol * yr-1"
         write(nfstdo, *) "'--------------------------------------------------'"
         if (allocated(var2d_values)) deallocate(var2d_values)

         allocate(var2d_values(nxdim,nydim,1))
         filename = "fopom.nc"
         varname = "OrgM_diffusion"
         varunit = "mMol m-3 * cm s-1"
         do j = 1, nydim
           do i = 1, nxdim
             ij = i + (j - 1) * nxdim
             var2d_values(i, j, 1) = fopomg(ij)
           end do
         end do
        !call write_netcdf_2dwt(filename, varname, varunit, var2d_values, lon_values, lat_values, yy_array)
         call calc_integral_and_mean_2d(fopomg, lon_values, lat_values, integral_value, mean_value)
         write(nfstdo, *) "'--------------------------------------------------'"
         write(nfstdo,'(A, A)') " Filename: ", trim(filename)
         write(nfstdo,'(A, A, A, I6)') " Variable: ", trim(varname), ", Year: ", YY
         write(nfstdo,'(A, E12.5, A)') " Global Averaged: ", mean_value, " mmol m-3 * cm s-1"
         write(nfstdo,'(A, E12.5, A)') " Global Integrated: ", integral_value * 3.1536d+2, " mol * yr-1"
         write(nfstdo, *) "'--------------------------------------------------'"
         if (allocated(var2d_values)) deallocate(var2d_values)

         allocate(var2d_values(nxdim,nydim,1))
         filename = "focar.nc"
         varname = "CaCO3_diffusion"
         varunit = "mMol m-3 * cm s-1"
         do j = 1, nydim
           do i = 1, nxdim
             ij = i + (j - 1) * nxdim
             var2d_values(i, j, 1) = focarg(ij)
           end do
         end do
        !call write_netcdf_2dwt(filename, varname, varunit, var2d_values, lon_values, lat_values, yy_array)
         call calc_integral_and_mean_2d(focarg, lon_values, lat_values, integral_value, mean_value)
         write(nfstdo, *) "'--------------------------------------------------'"
         write(nfstdo,'(A, A)') " Filename: ", trim(filename)
         write(nfstdo,'(A, A, A, I6)') " Variable: ", trim(varname), ", Year: ", YY
         write(nfstdo,'(A, E12.5, A)') " Global Averaged: ", mean_value, " mmol m-3 * cm s-1"
         write(nfstdo,'(A, E12.5, A)') " Global Integrated: ", integral_value * 3.1536d+2, " mol * yr-1"
         write(nfstdo, *) "'--------------------------------------------------'"
         if (allocated(var2d_values)) deallocate(var2d_values)

         allocate(var2d_values(nxdim,nydim,1))
         filename = "foopa.nc"
         varname = "SiO2_diffusion"
         varunit = "mMol m-3 * cm s-1"
         do j = 1, nydim
           do i = 1, nxdim
             ij = i + (j - 1) * nxdim
             var2d_values(i, j, 1) = foopag(ij)
           end do
         end do
        !call write_netcdf_2dwt(filename, varname, varunit, var2d_values, lon_values, lat_values, yy_array)
         call calc_integral_and_mean_2d(foopag, lon_values, lat_values, integral_value, mean_value)
         write(nfstdo, *) "'--------------------------------------------------'"
         write(nfstdo,'(A, A)') " Filename: ", trim(filename)
         write(nfstdo,'(A, A, A, I6)') " Variable: ", trim(varname), ", Year: ", YY
         write(nfstdo,'(A, E12.5, A)') " Global Averaged: ", mean_value, " mmol m-3 * cm s-1"
         write(nfstdo,'(A, E12.5, A)') " Global Integrated: ", integral_value * 3.1536d+2, " mol * yr-1"
         write(nfstdo, *) "'--------------------------------------------------'"
         if (allocated(var2d_values)) deallocate(var2d_values)

         allocate(var2d_values(nxdim,nydim,1))
         filename = "fiocla.nc"
         varname = "Clay_burial"
         varunit = "g m-3 * cm s-1"
         do j = 1, nydim
           do i = 1, nxdim
             ij = i + (j - 1) * nxdim
             var2d_values(i, j, 1) = fioclag(ij)
           end do
         end do
        !call write_netcdf_2dwt(filename, varname, varunit, var2d_values, lon_values, lat_values, yy_array)
         call calc_integral_and_mean_2d(fioclag, lon_values, lat_values, integral_value, mean_value)
         write(nfstdo, *) "'--------------------------------------------------'"
         write(nfstdo,'(A, A)') " Filename: ", trim(filename)
         write(nfstdo,'(A, A, A, I6)') " Variable: ", trim(varname), ", Year: ", YY
         write(nfstdo,'(A, E12.5, A)') " Global Averaged: ", mean_value, " g m-3 * cm s-1"
         write(nfstdo,'(A, E12.5, A)') " Global Integrated: ", integral_value * 3.1536d+5, " g * yr-1"
         write(nfstdo, *) "'--------------------------------------------------'"
         if (allocated(var2d_values)) deallocate(var2d_values)

         allocate(var2d_values(nxdim,nydim,1))
         filename = "fiopom.nc"
         varname = "OrgM_burial"
         varunit = "mMol m-3 * cm s-1"
         do j = 1, nydim
           do i = 1, nxdim
             ij = i + (j - 1) * nxdim
             var2d_values(i, j, 1) = fiopomg(ij)
           end do
         end do
         call write_netcdf_2dwt(filename, varname, varunit, var2d_values, lon_values, lat_values, yy_array)
         call calc_integral_and_mean_2d(fiopomg, lon_values, lat_values, integral_value, mean_value)
         write(nfstdo, *) "'--------------------------------------------------'"
         write(nfstdo,'(A, A)') " Filename: ", trim(filename)
         write(nfstdo,'(A, A, A, I6)') " Variable: ", trim(varname), ", Year: ", YY
         write(nfstdo,'(A, E12.5, A)') " Global Averaged: ", mean_value, " mmol m-3 * cm s-1"
         write(nfstdo,'(A, E12.5, A)') " Global Integrated: ", integral_value * 3.1536d+2, " mol * yr-1"
         write(nfstdo, *) "'--------------------------------------------------'"
         if (allocated(var2d_values)) deallocate(var2d_values)

         allocate(var2d_values(nxdim,nydim,1))
         filename = "fiocar.nc"
         varname = "CaCO3_burial"
         varunit = "mMol m-3 * cm s-1"
         do j = 1, nydim
           do i = 1, nxdim
             ij = i + (j - 1) * nxdim
             var2d_values(i, j, 1) = fiocarg(ij)
           end do
         end do
         call write_netcdf_2dwt(filename, varname, varunit, var2d_values, lon_values, lat_values, yy_array)
         call calc_integral_and_mean_2d(fiocarg, lon_values, lat_values, integral_value, mean_value)
         write(nfstdo, *) "'--------------------------------------------------'"
         write(nfstdo,'(A, A)') " Filename: ", trim(filename)
         write(nfstdo,'(A, A, A, I6)') " Variable: ", trim(varname), ", Year: ", YY
         write(nfstdo,'(A, E12.5, A)') " Global Averaged: ", mean_value, " mmol m-3 * cm s-1"
         write(nfstdo,'(A, E12.5, A)') " Global Integrated: ", integral_value * 3.1536d+2, " mol * yr-1"
         write(nfstdo, *) "'--------------------------------------------------'"
         if (allocated(var2d_values)) deallocate(var2d_values)

         allocate(var2d_values(nxdim,nydim,1))
         filename = "fioopa.nc"
         varname = "SiO2_burial"
         varunit = "mMol m-3 * cm s-1"
         do j = 1, nydim
           do i = 1, nxdim
             ij = i + (j - 1) * nxdim
             var2d_values(i, j, 1) = fioopag(ij)
           end do
         end do
         call write_netcdf_2dwt(filename, varname, varunit, var2d_values, lon_values, lat_values, yy_array)
         call calc_integral_and_mean_2d(fioopag, lon_values, lat_values, integral_value, mean_value)
         write(nfstdo, *) "'--------------------------------------------------'"
         write(nfstdo,'(A, A)') " Filename: ", trim(filename)
         write(nfstdo,'(A, A, A, I6)') " Variable: ", trim(varname), ", Year: ", YY
         write(nfstdo,'(A, E12.5, A)') " Global Averaged: ", mean_value, " mmol m-3 * cm s-1"
         write(nfstdo,'(A, E12.5, A)') " Global Integrated: ", integral_value * 3.1536d+2, " mol * yr-1"
         write(nfstdo, *) "'--------------------------------------------------'"
         if (allocated(var2d_values)) deallocate(var2d_values)

         allocate(var2d_values(nxdim,nydim,1))
         filename = "bfcla.nc"
         varname = "Clay_burial"
         varunit = "g m-3 * cm s-1"
         do j = 1, nydim
           do i = 1, nxdim
             ij = i + (j - 1) * nxdim
             var2d_values(i, j, 1) = bfclag(ij)
           end do
         end do
        !call write_netcdf_2dwt(filename, varname, varunit, var2d_values, lon_values, lat_values, yy_array)
         call calc_integral_and_mean_2d(bfclag, lon_values, lat_values, integral_value, mean_value)
         write(nfstdo, *) "'--------------------------------------------------'"
         write(nfstdo,'(A, A)') " Filename: ", trim(filename)
         write(nfstdo,'(A, A, A, I6)') " Variable: ", trim(varname), ", Year: ", YY
         write(nfstdo,'(A, E12.5, A)') " Global Averaged: ", mean_value, " mmol m-3 * cm s-1"
         write(nfstdo,'(A, E12.5, A)') " Global Integrated: ", integral_value * 3.1536d+2, " mol * yr-1"
         write(nfstdo, *) "'--------------------------------------------------'"
         if (allocated(var2d_values)) deallocate(var2d_values)

         allocate(var2d_values(nxdim,nydim,1))
         filename = "bfpom.nc"
         varname = "OrgM_burial"
         varunit = "mMol m-3 * cm s-1"
         do j = 1, nydim
           do i = 1, nxdim
             ij = i + (j - 1) * nxdim
             var2d_values(i, j, 1) = bfpomg(ij)
           end do
         end do
        !call write_netcdf_2dwt(filename, varname, varunit, var2d_values, lon_values, lat_values, yy_array)
         call calc_integral_and_mean_2d(bfpomg, lon_values, lat_values, integral_value, mean_value)
         write(nfstdo, *) "'--------------------------------------------------'"
         write(nfstdo,'(A, A)') " Filename: ", trim(filename)
         write(nfstdo,'(A, A, A, I6)') " Variable: ", trim(varname), ", Year: ", YY
         write(nfstdo,'(A, E12.5, A)') " Global Averaged: ", mean_value, " mmol m-3 * cm s-1"
         write(nfstdo,'(A, E12.5, A)') " Global Integrated: ", integral_value * 3.1536d+2, " mol * yr-1"
         write(nfstdo, *) "'--------------------------------------------------'"
         if (allocated(var2d_values)) deallocate(var2d_values)

         allocate(var2d_values(nxdim,nydim,1))
         filename = "bfcar.nc"
         varname = "CaCO3_burial"
         varunit = "mMol m-3 * cm s-1"
         do j = 1, nydim
           do i = 1, nxdim
             ij = i + (j - 1) * nxdim
             var2d_values(i, j, 1) = bfcarg(ij)
           end do
         end do
        !call write_netcdf_2dwt(filename, varname, varunit, var2d_values, lon_values, lat_values, yy_array)
         call calc_integral_and_mean_2d(bfcarg, lon_values, lat_values, integral_value, mean_value)
         write(nfstdo, *) "'--------------------------------------------------'"
         write(nfstdo,'(A, A)') " Filename: ", trim(filename)
         write(nfstdo,'(A, A, A, I6)') " Variable: ", trim(varname), ", Year: ", YY
         write(nfstdo,'(A, E12.5, A)') " Global Averaged: ", mean_value, " mmol m-3 * cm s-1"
         write(nfstdo,'(A, E12.5, A)') " Global Integrated: ", integral_value * 3.1536d+2, " mol * yr-1"
         write(nfstdo, *) "'--------------------------------------------------'"
         if (allocated(var2d_values)) deallocate(var2d_values)

         allocate(var2d_values(nxdim,nydim,1))
         filename = "bfopa.nc"
         varname = "SiO2_burial"
         varunit = "mMol m-3 * cm s-1"
         do j = 1, nydim
           do i = 1, nxdim
             ij = i + (j - 1) * nxdim
             var2d_values(i, j, 1) = bfopag(ij)
           end do
         end do
        !call write_netcdf_2dwt(filename, varname, varunit, var2d_values, lon_values, lat_values, yy_array)
         call calc_integral_and_mean_2d(bfopag, lon_values, lat_values, integral_value, mean_value)
         write(nfstdo, *) "'--------------------------------------------------'"
         write(nfstdo,'(A, A)') " Filename: ", trim(filename)
         write(nfstdo,'(A, A, A, I6)') " Variable: ", trim(varname), ", Year: ", YY
         write(nfstdo,'(A, E12.5, A)') " Global Averaged: ", mean_value, " mmol m-3 * cm s-1"
         write(nfstdo,'(A, E12.5, A)') " Global Integrated: ", integral_value * 3.1536d+2, " mol * yr-1"
         write(nfstdo, *) "'--------------------------------------------------'"
         if (allocated(var2d_values)) deallocate(var2d_values)

         allocate(var3d_values(nxdim,nydim,sdzdim+1,1))
         filename = "w.nc"
         varname = "Sediment_rate"
         varunit = "cm s-1"
         do k = 1, sdzdim+1
           do j = 1, nydim
             do i = 1, nxdim
               ij = i + (j - 1) * nxdim
               var3d_values(i, j, k, 1) = wg(ij, k)
             end do
           end do
         end do
         call write_netcdf_3dwt(filename, varname, varunit, var3d_values, lon_values, lat_values, lz+1, yy_array)
         call calc_integral_and_mean_3d(wg, lon_values, lat_values, lz+1, integral_value, mean_value)
         write(nfstdo, *) "'--------------------------------------------------'"
         write(nfstdo,'(A, A)') " Filename: ", trim(filename)
         write(nfstdo,'(A, A, A, I6)') " Variable: ", trim(varname), ", Year: ", YY
         write(nfstdo,'(A, E12.5, A)') " Global Averaged: ", mean_value, " cm s-1"
         write(nfstdo, *) "'--------------------------------------------------'"
         if (allocated(var3d_values)) deallocate(var3d_values)

        !allocate(var3d_values(nxdim,nydim,sdzdim,1))
        !filename = "cdc13.nc"
        !varname = "delta_carbon-13"
        !varunit = "permil"
        !do k = 1, sdzdim
        !  do j = 1, nydim
        !    do i = 1, nxdim
        !      ij = i + (j - 1) * nxdim
        !      var3d_values(i, j, k, 1) = cdc13g(ij, k)
        !    end do
        !  end do
        !end do
        !call write_netcdf_3dwt(filename, varname, varunit, var3d_values, lon_values, lat_values, lz, yy_array)
        !call calc_integral_and_mean_3d(cdc13g, lon_values, lat_values, lz, integral_value, mean_value)
        !write(nfstdo, *) "'--------------------------------------------------'"
        !write(nfstdo,'(A, A)') " Filename: ", trim(filename)
        !write(nfstdo,'(A, A, A, I6)') " Variable: ", trim(varname), ", Year: ", YY
        !write(nfstdo,'(A, E12.5, A)') " Global Averaged: ", mean_value, " mmol m-3"
        !write(nfstdo,'(A, E12.5, A)') " Global Integrated: ", integral_value * 1.d-3, " mol"
        !write(nfstdo, *) "'--------------------------------------------------'"
        !if (allocated(var3d_values)) deallocate(var3d_values)

        !allocate(var3d_values(nxdim,nydim,sdzdim,1))
        !filename = "cdc14.nc"
        !varname = "Delta_carbon-14"
        !varunit = "permil"
        !do k = 1, sdzdim+1
        !  do j = 1, nydim
        !    do i = 1, nxdim
        !      ij = i + (j - 1) * nxdim
        !      var3d_values(i, j, k, 1) = cdc14g(ij, k)
        !    end do
        !  end do
        !end do
        !call write_netcdf_3dwt(filename, varname, varunit, var3d_values, lon_values, lat_values, lz, yy_array)
        !call calc_integral_and_mean_3d(cdc14g, lon_values, lat_values, lz, integral_value, mean_value)
        !write(nfstdo, *) "'--------------------------------------------------'"
        !write(nfstdo,'(A, A)') " Filename: ", trim(filename)
        !write(nfstdo,'(A, A, A, I6)') " Variable: ", trim(varname), ", Year: ", YY
        !write(nfstdo,'(A, E12.5, A)') " Global Averaged: ", mean_value, " mmol m-3"
        !write(nfstdo,'(A, E12.5, A)') " Global Integrated: ", integral_value * 1.d-3, " mol"
        !write(nfstdo, *) "'--------------------------------------------------'"
        !if (allocated(var3d_values)) deallocate(var3d_values)

        !! -- Gtool output --
        !chead  = ''
        !cfile  = 'scla'
        !citem  = 'scla'
        !ctitle = 'Clay'
        !cunit  = 'g/m^3'
        !call write_gtool_beta( cfile, chead, citem, ctitle, cunit, nxdim, nydim, sdzdim, 1, sclag)
        !cfile  = 'spom'
        !citem  = 'spom'
        !ctitle = 'OrgM'
        !cunit  = 'mmol/m^3'
        !call write_gtool_beta( cfile, chead, citem, ctitle, cunit, nxdim, nydim, sdzdim, 1, spomg)
        !cfile  = 'scar'
        !citem  = 'scar'
        !ctitle = 'CaCO3'
        !cunit  = 'mmol/m^3'
        !call write_gtool_beta( cfile, chead, citem, ctitle, cunit, nxdim, nydim, sdzdim, 1, scarg)
        !cfile  = 'sopa'
        !citem  = 'sopa'
        !ctitle = 'SiO2'
        !cunit  = 'mmol/m^3'
        !call write_gtool_beta( cfile, chead, citem, ctitle, cunit, nxdim, nydim, sdzdim, 1, sopag)
        !cfile  = 'sage'
        !citem  = 'sage'
        !ctitle = 'Age tracer'
        !cunit  = 'year'
        !call write_gtool_beta( cfile, chead, citem, ctitle, cunit, nxdim, nydim, sdzdim, 1, sageg)

        !cfile  = 'cnit'
        !citem  = 'cnit'
        !ctitle = 'Phosphate'
        !cunit  = 'mmol/m^3'
        !call write_gtool_beta( cfile, chead, citem, ctitle, cunit, nxdim, nydim, sdzdim, 1, cnitg)
        !cfile  = 'cdic'
        !citem  = 'cdic'
        !ctitle = 'dic'
        !cunit  = 'mmol/m^3'
        !call write_gtool_beta( cfile, chead, citem, ctitle, cunit, nxdim, nydim, sdzdim, 1, cdicg)
        !cfile  = 'calk'
        !citem  = 'calk'
        !ctitle = 'Alkalinity'
        !cunit  = 'mmol/m^3'
        !call write_gtool_beta( cfile, chead, citem, ctitle, cunit, nxdim, nydim, sdzdim, 1, calkg)
        !cfile  = 'cdo2'
        !citem  = 'cdo2'
        !ctitle = 'Oxygen'
        !cunit  = 'mmol/m^3'
        !call write_gtool_beta( cfile, chead, citem, ctitle, cunit, nxdim, nydim, sdzdim, 1, cdo2g)
        !cfile  = 'csil'
        !citem  = 'csil'
        !ctitle = 'Silicate'
        !cunit  = 'mmol/m^3'
        !call write_gtool_beta( cfile, chead, citem, ctitle, cunit, nxdim, nydim, sdzdim, 1, csilg)

        !cfile  = 'ficla'
        !citem  = 'ficla'
        !ctitle = 'Clay bottom flux'
        !cunit  = 'g/m^3*cm/s'
        !call write_gtool_beta( cfile, chead, citem, ctitle, cunit, nxdim, nydim, 1, 1, ficlag)
        !cfile  = 'fipom'
        !citem  = 'fipom'
        !ctitle = 'OrgM bottom flux'
        !cunit  = 'mMol/m^3*cm/s'
        !call write_gtool_beta( cfile, chead, citem, ctitle, cunit, nxdim, nydim, 1, 1, fipomg)
        !cfile  = 'ficar'
        !citem  = 'ficar'
        !ctitle = 'CaCO3 bottom flux'
        !cunit  = 'mMol/m^3*cm/s'
        !call write_gtool_beta( cfile, chead, citem, ctitle, cunit, nxdim, nydim, 1, 1, ficarg)
        !cfile  = 'fiopa'
        !citem  = 'fiopa'
        !ctitle = 'SiO2 bottom flux'
        !cunit  = 'mMol/m^3*cm/s'
        !call write_gtool_beta( cfile, chead, citem, ctitle, cunit, nxdim, nydim, 1, 1, fiopag)

        !cfile  = 'fopom'
        !citem  = 'fopom'
        !ctitle = 'OrgM diffusion'
        !cunit  = 'mMol/m^3*cm/s'
        !call write_gtool_beta( cfile, chead, citem, ctitle, cunit, nxdim, nydim, 1, 1, fopomg)
        !cfile  = 'focar'
        !citem  = 'focar'
        !ctitle = 'CaCO3 diffusion'
        !cunit  = 'mMol/m^3*cm/s'
        !call write_gtool_beta( cfile, chead, citem, ctitle, cunit, nxdim, nydim, 1, 1, focarg)
        !cfile  = 'foopa'
        !citem  = 'foopa'
        !ctitle = 'SiO2 diffusion'
        !cunit  = 'mMol/m^3*cm/s'
        !call write_gtool_beta( cfile, chead, citem, ctitle, cunit, nxdim, nydim, 1, 1, foopag)

        !cfile  = 'bfcla'
        !citem  = 'bfcla'
        !ctitle = 'Clay burial'
        !cunit  = 'g/m^3*cm/s'
        !call write_gtool_beta( cfile, chead, citem, ctitle, cunit, nxdim, nydim, 1, 1, bfclag)
        !cfile  = 'bfpom'
        !citem  = 'bfpom'
        !ctitle = 'OrgM burial'
        !cunit  = 'mMol/m^3*cm/s'
        !call write_gtool_beta( cfile, chead, citem, ctitle, cunit, nxdim, nydim, 1, 1, bfpomg)
        !cfile  = 'bfcar'
        !citem  = 'bfcar'
        !ctitle = 'CaCO3 burial'
        !cunit  = 'mMol/m^3*cm/s'
        !call write_gtool_beta( cfile, chead, citem, ctitle, cunit, nxdim, nydim, 1, 1, bfcarg)
        !cfile  = 'bfopa'
        !citem  = 'bfopa'
        !ctitle = 'SiO2 burial'
        !cunit  = 'mMol/m^3*cm/s'
        !call write_gtool_beta( cfile, chead, citem, ctitle, cunit, nxdim, nydim, 1, 1, bfopag)

        !cfile  = 'w'
        !citem  = 'w'
        !ctitle = 'Sediment rate'
        !cunit  = 'cm/s'
        !call write_gtool_beta( cfile, chead, citem, ctitle, cunit, nxdim, nydim, sdzdim+1, 1, wg)

        !cfile  = 'cdc13'
        !citem  = 'cdc13'
        !ctitle = 'delta carbon-13'
        !cunit  = 'permil'
        !call write_gtool_beta( cfile, chead, citem, ctitle, cunit, nxdim, nydim, sdzdim, 1, cdc13g)
        !cfile  = 'cdc14'
        !citem  = 'cdc14'
        !ctitle = 'Delta carbon-14'
        !cunit  = 'permil'
        !call write_gtool_beta( cfile, chead, citem, ctitle, cunit, nxdim, nydim, sdzdim, 1, cdc14g)

         ! -- Pass the tracer value to the output program --
         do k = skst, sken
            do ij = 1, nxydim
               sdtx(ij, k, 1) = sclag(ij, k)
               sdtx(ij, k, 2) = spomg(ij, k)
               sdtx(ij, k, 3) = scarg(ij, k)
               sdtx(ij, k, 4) = sopag(ij, k)
               sdtx(ij, k, 5) = sageg(ij, k)
               sdtx(ij, k, 6) = cnitg(ij, k)
               sdtx(ij, k, 7) = cdicg(ij, k)
               sdtx(ij, k, 8) = calkg(ij, k)
               sdtx(ij, k, 9) = csilg(ij, k)
               sdtx(ij, k, 10) = cdo2g(ij, k)
            if ( calc13 .and. sedtdm == 16 ) then
               sdtx(ij, k, 11) = s13pomg(ij, k)
               sdtx(ij, k, 12) = s13carg(ij, k)
               sdtx(ij, k, 13) = c13dicg(ij, k)
            end if
            if ( calc14 .and. sedtdm == 16 ) then
               sdtx(ij, k, 14) = s14pomg(ij, k)
               sdtx(ij, k, 15) = s14carg(ij, k)
               sdtx(ij, k, 16) = c14dicg(ij, k)
            end if
            end do
         end do

      end if

      call mpi_bcast(sdtx, nxydim*sdzdim*sedtdm, mpi_real8, iroot, mpi_comm_world, ierr)
      call mpi_barrier(mpi_comm_world, ierr)

end if ! -- End of use_mpi

      call clcend('PRDC_HIST')

      end if ! -- End of Output

      end do ! -- End of Time integration

      call clcend('PRDC_STEP')

      write(nfstdo, *) "'--------------------------------------------------'"
      write(nfstdo, *) "'-- DIGDA Prediction:                            --'"
      write(nfstdo, *) "'-- Finished time integration                    --'"
      write(nfstdo, *) "'--------------------------------------------------'"
      call flush(nfstdo)

      if (hstavg) then
         if (allocated(mean_data)) deallocate(mean_data)
      end if

      end subroutine digda_cprdc

      subroutine calc_integral_and_mean_2d(var_values, lon_values, lat_values, integral_value, mean_value)

        double precision, dimension(:), intent(in) :: var_values
        double precision, dimension(:), intent(in) :: lon_values, lat_values
        double precision, intent(out) :: integral_value, mean_value
        double precision :: grid_volume
        integer(kind=4) :: nx, ny, nz
        double precision, dimension(:), allocatable :: lon_grid_size, lat_grid_size
        double precision, dimension(:), allocatable :: lon_grid_size_m, lat_grid_size_m
        integer :: i, j, ij

        nx = size(lon_values)
        ny = size(lat_values)
        allocate(lon_grid_size(nx))
        allocate(lat_grid_size(ny))
        allocate(lon_grid_size_m(nx*ny))
        allocate(lat_grid_size_m(nx*ny))

        ! Calculate grid sizes
          lon_grid_size(1) = (lon_values(1) - 0.d0) * 2.d0
        do i = 2, nx
          lon_grid_size(i) = ( abs(lon_values(i) - lon_values(i-1)) &
                             - 0.5d0 * lon_grid_size(i-1) ) &
                             * 2.d0
        end do
          lat_grid_size(1) = (lat_values(1) + 90.d0) * 2.d0
        do j = 2, ny
          lat_grid_size(j) = ( abs(lat_values(j) - lat_values(j-1)) &
                             - 0.5d0 * lat_grid_size(j-1) ) &
                             * 2.d0
        end do

        ! Convert grid sizes to meters
        do j = 1, ny
          do i = 1, nx
            ij = i + (j-1)*nx
            lon_grid_size_m(ij) = lon_grid_size(i) * (pi/180.d00) * er * cos(lat_values(j) * (pi/180.d00))
            lat_grid_size_m(ij) = lat_grid_size(j) * (pi/180.d00) * er
          end do
        end do

        grid_volume = 0.d0
        integral_value = 0.d0
        mean_value = 0.d0
        do j = 1, ny
          do i = 1, nx
            ij = i + (j-1)*nx
            grid_volume = grid_volume &
              + lon_grid_size_m(ij) &
              * lat_grid_size_m(ij)
            integral_value = integral_value &
              + var_values(ij) &
              * lon_grid_size_m(ij) &
              * lat_grid_size_m(ij)
          end do
        end do
        mean_value = integral_value / grid_volume

        if (allocated(lon_grid_size))   deallocate(lon_grid_size)
        if (allocated(lat_grid_size))   deallocate(lat_grid_size)
        if (allocated(lon_grid_size_m)) deallocate(lon_grid_size_m)
        if (allocated(lat_grid_size_m)) deallocate(lat_grid_size_m)

      end subroutine calc_integral_and_mean_2d

      subroutine calc_integral_and_mean_3d(var_values, lon_values, lat_values, dep_values, integral_value, mean_value)

        double precision, dimension(:,:), intent(in) :: var_values
        double precision, dimension(:), intent(in) :: lon_values, lat_values, dep_values
        double precision, intent(out) :: integral_value, mean_value
        double precision :: grid_volume
        integer(kind=4) :: nx, ny, nz
        double precision, dimension(:), allocatable :: lon_grid_size, lat_grid_size, dep_grid_size
        double precision, dimension(:), allocatable :: lon_grid_size_m, lat_grid_size_m, dep_grid_size_m
        integer :: i, j, k, ij

        nx = size(lon_values)
        ny = size(lat_values)
        nz = size(dep_values)
        allocate(lon_grid_size(nx))
        allocate(lat_grid_size(ny))
        allocate(dep_grid_size(nz))
        allocate(lon_grid_size_m(nx*ny))
        allocate(lat_grid_size_m(nx*ny))
        allocate(dep_grid_size_m(nz))

        ! Calculate grid sizes
          lon_grid_size(1) = (lon_values(1) - 0.d0) * 2.d0
        do i = 2, nx
          lon_grid_size(i) = ( abs(lon_values(i) - lon_values(i-1)) &
                             - 0.5d0 * lon_grid_size(i-1) ) &
                             * 2.d0
        end do
          lat_grid_size(1) = (lat_values(1) + 90.d0) * 2.d0
        do j = 2, ny
          lat_grid_size(j) = ( abs(lat_values(j) - lat_values(j-1)) &
                             - 0.5d0 * lat_grid_size(j-1) ) &
                             * 2.d0
        end do
          dep_grid_size(1) = (dep_values(1) - 0.d0) * 2.d0
        do k = 2, nz
          dep_grid_size(k) = ( abs(dep_values(k) - dep_values(k-1)) &
                             - 0.5d0 * dep_grid_size(k-1) ) &
                             * 2.d0
        end do

        ! Convert grid sizes to meters
        do j = 1, ny
          do i = 1, nx
            ij = i + (j-1)*nx
            lon_grid_size_m(ij) = lon_grid_size(i) * (pi/180.d00) * er * cos(lat_values(j) * (pi/180.d00))
            lat_grid_size_m(ij) = lat_grid_size(j) * (pi/180.d00) * er
          end do
        end do
        do k = 1, nz
          dep_grid_size_m(k) = dep_grid_size(k) * 1.d-2
        end do

        grid_volume = 0.d0
        integral_value = 0.d0
        mean_value = 0.d0
        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              ij = i + (j-1)*nx
              grid_volume = grid_volume &
                + lon_grid_size_m(ij) &
                * lat_grid_size_m(ij) &
                * dep_grid_size_m(k)
              integral_value = integral_value &
                + var_values(ij, k) &
                * lon_grid_size_m(ij) &
                * lat_grid_size_m(ij) &
                * dep_grid_size_m(k)
            end do
          end do
        end do
        mean_value = integral_value / grid_volume

        if (allocated(lon_grid_size))   deallocate(lon_grid_size)
        if (allocated(lat_grid_size))   deallocate(lat_grid_size)
        if (allocated(dep_grid_size))   deallocate(dep_grid_size)
        if (allocated(lon_grid_size_m)) deallocate(lon_grid_size_m)
        if (allocated(lat_grid_size_m)) deallocate(lat_grid_size_m)
        if (allocated(dep_grid_size_m)) deallocate(dep_grid_size_m)

      end subroutine calc_integral_and_mean_3d

end module module_digda_cprdc
