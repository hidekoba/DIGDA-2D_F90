!######################################################################
module module_digda_grid
  use module_digda_mpi
  implicit none

  ! Declaration of parameters and variables
  integer, parameter :: benfdm = 14
  ! cla,pom,car,opa,age,tem,sal,nit,dic,alk,sil,do2,13dic,14dic
  integer, parameter :: sedtdm = 10
 !integer, parameter :: sedtdm = 16
  ! scla,spom,scar,sopa,sage,cnit,cdic,calk,csil,cdo2,s13pom,s13car,c13dic,s14pom,s14car,c14dic

  ! Vertical grid settings
  !-- ZGRID_CONSTANT --
 !integer, parameter :: sdzdim = 20
 !integer, parameter :: skst = 1
 !integer, parameter :: sken = 20
 !integer :: sztdim = sdzdim * sedtdm
 !double precision :: lz(skst:sken) = &
 !  (/ 2.500d-3, 7.500d-3, 1.250d-2, 1.750d-2, 2.250d-2, &
 !     2.750d-2, 3.250d-2, 3.750d-2, 4.250d-2, 4.750d-2, &
 !     5.250d-2, 5.750d-2, 6.250d-2, 6.750d-2, 7.250d-2, &
 !     7.750d-2, 8.250d-2, 8.750d-2, 9.250d-2, 9.750d-2 /)
 !double precision :: lzi(skst:sken+1) = &
 !  (/ 0.0000d0, 5.000d-3, 1.000d-2, 1.500d-2, 2.000d-2, &
 !     2.500d-2, 3.000d-2, 3.500d-2, 4.000d-2, 4.500d-2, &
 !     5.000d-2, 5.500d-2, 6.000d-2, 6.500d-2, 7.000d-2, &
 !     7.500d-2, 8.000d-2, 8.500d-2, 9.000d-2, 9.500d-2, &
 !     1.000d-1 /)
 !double precision :: dlz(skst:sken) = &
 !  (/ 5.000d-3, 5.000d-3, 5.000d-3, 5.000d-3, 5.000d-3, &
 !     5.000d-3, 5.000d-3, 5.000d-3, 5.000d-3, 5.000d-3, &
 !     5.000d-3, 5.000d-3, 5.000d-3, 5.000d-3, 5.000d-3, &
 !     5.000d-3, 5.000d-3, 5.000d-3, 5.000d-3, 5.000d-3 /)
 !double precision :: dlzi(skst:sken+1) = &
 !  (/ 5.000d-3, 5.000d-3, 5.000d-3, 5.000d-3, 5.000d-3, &
 !     5.000d-3, 5.000d-3, 5.000d-3, 5.000d-3, 5.000d-3, &
 !     5.000d-3, 5.000d-3, 5.000d-3, 5.000d-3, 5.000d-3, &
 !     5.000d-3, 5.000d-3, 5.000d-3, 5.000d-3, 5.000d-3, &
 !     5.000d-3 /)

  !-- ZGRID_ARCHER --
 !integer, parameter :: sdzdim = 16
 !integer, parameter :: skst = 1
 !integer, parameter :: sken = 16
 !integer :: sztdim = sdzdim * sedtdm
 !double precision :: lz(skst:sken) = &
 !  (/ 2.500d-4, 1.750d-3, 5.000d-3, 1.100d-2, 2.000d-2, &
 !     3.250d-2, 5.000d-2, 7.250d-2, 1.025d-1, 1.400d-1, &
 !     1.850d-1, 2.300d-1, 2.750d-1, 4.000d-1, 6.250d-1, &
 !     8.750d-1 /)
 !double precision :: lzi(skst:sken+1) = &
 !  (/ 0.000d0,  5.000d-4, 3.000d-3, 7.000d-3, 1.500d-2, &
 !     2.500d-2, 4.000d-2, 6.000d-2, 8.500d-2, 1.200d-1, &
 !     1.600d-1, 2.100d-1, 2.500d-1, 3.000d-1, 5.000d-1, &
 !     7.500d-1, 1.000d0 /)
 !double precision :: dlz(skst:sken) = &
 !  (/ 5.000d-4, 2.500d-3, 4.000d-3, 8.000d-3, 1.000d-2, &
 !     1.500d-2, 2.000d-2, 2.500d-2, 3.500d-2, 4.000d-2, &
 !     5.000d-2, 4.000d-2, 5.000d-2, 2.000d-1, 2.500d-1, &
 !     2.500d-1 /)
 !double precision :: dlzi(skst:sken+1) = &
 !  (/ 5.000d-4, 1.500d-3, 3.250d-3, 6.000d-3, 9.000d-3, &
 !     1.250d-2, 1.750d-2, 2.250d-2, 3.000d-2, 3.750d-2, &
 !     4.500d-2, 4.500d-2, 4.500d-2, 1.250d-1, 2.250d-1, &
 !     2.500d-1, 2.500d-1 /)

  !-- ZGRID_HEInzE --
  integer, parameter :: sdzdim = 10
  integer, parameter :: skst = 1
  integer, parameter :: sken = 10
  integer, parameter :: sztdim = sdzdim * sedtdm
  double precision :: lz(skst:sken) = &
    (/ 1.500d-3, 4.500d-3, 8.500d-3, 1.350d-2, 1.850d-2, &
       2.600d-2, 3.600d-2, 4.600d-2, 6.325d-2, 8.775d-2 /)
  double precision :: lzi(skst:sken+1) = &
    (/ 0.0000d0, 3.000d-3, 6.000d-3, 1.100d-2, 1.600d-2, &
       2.100d-2, 3.100d-2, 4.100d-2, 5.100d-2, 7.550d-2, &
       1.000d-1 /)
  double precision :: dlz(skst:sken) = &
    (/ 3.000d-3, 3.000d-3, 5.000d-3, 5.000d-3, 5.000d-3, &
       1.000d-2, 1.000d-2, 1.000d-2, 2.450d-2, 2.450d-2 /)
  double precision :: dlzi(skst:sken+1) = &
    (/ 3.000d-3, 3.000d-3, 4.000d-3, 5.000d-3, 5.000d-3, &
       7.500d-3, 1.000d-2, 1.000d-2, 1.725d-2, 2.450d-2, &
       2.450d-2 /)

  ! Horizontal grid settings
  !-- ONE_Dimension --
 !integer, parameter :: nxdim = 1
 !integer, parameter :: nydim = 1
 !integer, parameter :: nzdim = 44

  !-- Grid_128x64 --
 !integer, parameter :: nxdim = 128
 !integer, parameter :: nydim = 64
 !integer, parameter :: nzdim = 44

  !-- Grid_256x192 --
  integer, parameter :: nxdim = 256
  integer, parameter :: nydim = 192
  integer, parameter :: nzdim = 44

  !-- Grid_360x180 --
 !integer, parameter :: nxdim = 360
 !integer, parameter :: nydim = 180
 !integer, parameter :: nzdim = 44

 !!-- Grid_2160x1080 --
 !integer, parameter :: nxdim = 2160
 !integer, parameter :: nydim = 1080
 !integer, parameter :: nzdim = 44

  character(len=100) :: bound_dir = 'not-specified/'
  namelist /nmocbc/ bound_dir
  character(len=100) :: init_dir = 'not-specified/'
  character(len=100) :: rstt_dir = 'not-specified/'
  namelist /nmocrt/ init_dir, rstt_dir

  integer, parameter :: nxydim = nxdim * nydim
  integer, parameter :: nx = nxdim
  integer, parameter :: ny = nydim / ncpu
  integer, parameter :: nz = nzdim
  integer, parameter :: sdz = sdzdim
  integer, parameter :: nxy = nx * ny
  integer, parameter :: nxyz = nxy * sdz
  integer, parameter :: nxyzdm = nxydim * sdz

  integer, parameter :: istr = 1
  integer, parameter :: iend = nx
  integer, parameter :: jstr = 1
  integer, parameter :: jend = ny
  integer, parameter :: kstr = 1
  integer, parameter :: kend = nz
  integer, parameter :: ijtstr = 1
  integer, parameter :: ijtend = nxy

  ! Ocean depth settings
  double precision :: odepth(nzdim)
  integer :: nbot(nxydim)
  double precision :: dz(nxydim,nzdim)

  ! Initialize ocean depth values if nzdim is 44
  data odepth /  5.0d0,   10.0d0,   15.0d0,   20.0d0,   25.0d0, &
                30.0d0,   37.5d0,   45.0d0,   55.0d0,   65.0d0, &
                75.0d0,   85.0d0,  100.0d0,  120.0d0,  145.0d0, &
               170.0d0,  200.0d0,  250.0d0,  325.0d0,  400.0d0, &
               500.0d0,  600.0d0,  750.0d0,  900.0d0, 1050.0d0, &
              1200.0d0, 1400.0d0, 1650.0d0, 1900.0d0, 2150.0d0, &
              2400.0d0, 2650.0d0, 2900.0d0, 3150.0d0, 3400.0d0, &
              3650.0d0, 3900.0d0, 4150.0d0, 4400.0d0, 4650.0d0, &
              4900.0d0, 5150.0d0, 5400.0d0, 5500.0d0 /

  ! Missing value
  double precision, parameter :: rmiss = 999.00d0
 !double precision, parameter :: rmiss = 0.0000d0

end module module_digda_grid
