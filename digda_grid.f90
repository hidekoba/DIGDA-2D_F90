!######################################################################
module module_digda_grid
  use module_digda_config
  use mpi
  implicit none

  ! Declaration of parameters and variables
  integer, parameter :: IROOT = 0
  integer :: IERR, NPROCS, MYRANK
  integer :: IREQ, ISTATUS
 !integer, parameter :: NCPU = 12
  integer, parameter :: NCPU = 16
 !integer, parameter :: NCPU = 24
  integer :: IRANK, ISIZE1, ISIZE2, ISIZE3, IIST, IIEN, IJKG
  integer :: ISEND1, ISEND2, ISEND3, ISEND4, ISEND5
  integer :: IRECV1, IRECV2, IRECV3, IRECV4, IRECV5
  integer :: IREQ1, IREQ2, IREQ3, IREQ4, IREQ5
  integer :: IREQ11, IREQ12, IREQ13, IREQ14, IREQ15
  integer :: IREQ21, IREQ22, IREQ23, IREQ24, IREQ25
  integer :: JJSTR(0:NCPU-1), JJEND(0:NCPU-1)

  integer, parameter :: BENFDM = 14
  ! cla,pom,car,opa,age,tem,sal,nit,dic,alk,sil,do2,13dic,14dic,
  integer, parameter :: SEDTDM = 16
  ! scla,spom,scar,sopa,sage,cnit,cdic,calk,csil,cdo2,s13pom,s13car,c13dic,s14pom,s14car,c14dic,

  ! Vertical grid settings
  ! ZGRID_CONSTANT
!!integer, parameter :: SDZDIM = 20
!!integer, parameter :: SKST = 1
!!integer, parameter :: SKEN = 20
!!integer :: SZTDIM = SDZDIM * SEDTDM
!!double precision :: LZ(SKST:SKEN) = &
!!  (/ 2.500D-3, 7.500D-3, 1.250D-2, 1.750D-2, 2.250D-2, &
!!     2.750D-2, 3.250D-2, 3.750D-2, 4.250D-2, 4.750D-2, &
!!     5.250D-2, 5.750D-2, 6.250D-2, 6.750D-2, 7.250D-2, &
!!     7.750D-2, 8.250D-2, 8.750D-2, 9.250D-2, 9.750D-2 /)
!!double precision :: LZI(SKST:SKEN+1) = &
!!  (/ 0.0000D0, 5.000D-3, 1.000D-2, 1.500D-2, 2.000D-2, &
!!     2.500D-2, 3.000D-2, 3.500D-2, 4.000D-2, 4.500D-2, &
!!     5.000D-2, 5.500D-2, 6.000D-2, 6.500D-2, 7.000D-2, &
!!     7.500D-2, 8.000D-2, 8.500D-2, 9.000D-2, 9.500D-2, &
!!     1.000D-1 /)
!!double precision :: DLZ(SKST:SKEN) = &
!!  (/ 5.000D-3, 5.000D-3, 5.000D-3, 5.000D-3, 5.000D-3, &
!!     5.000D-3, 5.000D-3, 5.000D-3, 5.000D-3, 5.000D-3, &
!!     5.000D-3, 5.000D-3, 5.000D-3, 5.000D-3, 5.000D-3, &
!!     5.000D-3, 5.000D-3, 5.000D-3, 5.000D-3, 5.000D-3 /)
!!double precision :: DLZI(SKST:SKEN+1) = &
!!  (/ 5.000D-3, 5.000D-3, 5.000D-3, 5.000D-3, 5.000D-3, &
!!     5.000D-3, 5.000D-3, 5.000D-3, 5.000D-3, 5.000D-3, &
!!     5.000D-3, 5.000D-3, 5.000D-3, 5.000D-3, 5.000D-3, &
!!     5.000D-3, 5.000D-3, 5.000D-3, 5.000D-3, 5.000D-3, &
!!     5.000D-3 /)

  ! ZGRID_ARCHER
!!integer, parameter :: SDZDIM = 16
!!integer, parameter :: SKST = 1
!!integer, parameter :: SKEN = 16
!!integer :: SZTDIM = SDZDIM * SEDTDM
!!double precision :: LZ(SKST:SKEN) = &
!!  (/ 2.500D-4, 1.750D-3, 5.000D-3, 1.100D-2, 2.000D-2, &
!!     3.250D-2, 5.000D-2, 7.250D-2, 1.025D-1, 1.400D-1, &
!!     1.850D-1, 2.300D-1, 2.750D-1, 4.000D-1, 6.250D-1, &
!!     8.750D-1 /)
!!double precision :: LZI(SKST:SKEN+1) = &
!!  (/ 0.000D0,  5.000D-4, 3.000D-3, 7.000D-3, 1.500D-2, &
!!     2.500D-2, 4.000D-2, 6.000D-2, 8.500D-2, 1.200D-1, &
!!     1.600D-1, 2.100D-1, 2.500D-1, 3.000D-1, 5.000D-1, &
!!     7.500D-1, 1.000D0 /)
!!double precision :: DLZ(SKST:SKEN) = &
!!  (/ 5.000D-4, 2.500D-3, 4.000D-3, 8.000D-3, 1.000D-2, &
!!     1.500D-2, 2.000D-2, 2.500D-2, 3.500D-2, 4.000D-2, &
!!     5.000D-2, 4.000D-2, 5.000D-2, 2.000D-1, 2.500D-1, &
!!     2.500D-1 /)
!!double precision :: DLZI(SKST:SKEN+1) = &
!!  (/ 5.000D-4, 1.500D-3, 3.250D-3, 6.000D-3, 9.000D-3, &
!!     1.250D-2, 1.750D-2, 2.250D-2, 3.000D-2, 3.750D-2, &
!!     4.500D-2, 4.500D-2, 4.500D-2, 1.250D-1, 2.250D-1, &
!!     2.500D-1, 2.500D-1 /)

  ! ZGRID_HEINZE
  integer, parameter :: SDZDIM = 10
  integer, parameter :: SKST = 1
  integer, parameter :: SKEN = 10
  integer, parameter :: SZTDIM = SDZDIM * SEDTDM
  double precision :: LZ(SKST:SKEN) = &
    (/ 1.500D-3, 4.500D-3, 8.500D-3, 1.350D-2, 1.850D-2, &
       2.600D-2, 3.600D-2, 4.600D-2, 6.325D-2, 8.775D-2 /)
  double precision :: LZI(SKST:SKEN+1) = &
    (/ 0.0000D0, 3.000D-3, 6.000D-3, 1.100D-2, 1.600D-2, &
       2.100D-2, 3.100D-2, 4.100D-2, 5.100D-2, 7.550D-2, &
       1.000D-1 /)
  double precision :: DLZ(SKST:SKEN) = &
    (/ 3.000D-3, 3.000D-3, 5.000D-3, 5.000D-3, 5.000D-3, &
       1.000D-2, 1.000D-2, 1.000D-2, 2.450D-2, 2.450D-2 /)
  double precision :: DLZI(SKST:SKEN+1) = &
    (/ 3.000D-3, 3.000D-3, 4.000D-3, 5.000D-3, 5.000D-3, &
       7.500D-3, 1.000D-2, 1.000D-2, 1.725D-2, 2.450D-2, &
       2.450D-2 /)

  ! Horizontal grid settings
  ! ONE_DIM
!!integer, parameter :: NXDIM = 1
!!integer, parameter :: NYDIM = 1
!!integer, parameter :: NZDIM = 44
!!character(len=100), parameter :: BOUND_DIR = "bound1d/"

  ! Grid_128x64
!!integer, parameter :: NXDIM = 128
!!integer, parameter :: NYDIM = 64
!!integer, parameter :: NZDIM = 44
!!character(len=100), parameter :: BOUND_DIR = "bound12864/"

  ! Grid_256x192
  integer, parameter :: NXDIM = 256
  integer, parameter :: NYDIM = 192
  integer, parameter :: NZDIM = 44
  character(len=100), parameter :: BOUND_DIR = "bound256192/"

!!! Grid_360x180
!!integer, parameter :: NXDIM = 360
!!integer, parameter :: NYDIM = 180
!!integer, parameter :: NZDIM = 44
!!character(len=100), parameter :: BOUND_DIR = "bound360180/"

!!! Grid_2160x1080
!!integer, parameter :: NXDIM = 2160
!!integer, parameter :: NYDIM = 1080
!!integer, parameter :: NZDIM = 44
!!character(len=100), parameter :: BOUND_DIR = "bound21601080/reversed_files/"

  integer, parameter :: NXYDIM = NXDIM * NYDIM
  integer, parameter :: NX = NXDIM
  integer, parameter :: NY = NYDIM / NCPU
  integer, parameter :: NZ = NZDIM
  integer, parameter :: SDZ = SDZDIM
  integer, parameter :: NXY = NX * NY
  integer, parameter :: NXYZ = NXY * SDZ
  integer, parameter :: NXYZDM = NXYDIM * SDZ

  integer, parameter :: ISTR = 1
  integer, parameter :: IEND = NX
  integer, parameter :: JSTR = 1
  integer, parameter :: JEND = NY
  integer, parameter :: KSTR = 1
  integer, parameter :: KEND = NZ
  integer, parameter :: IJTSTR = 1
  integer, parameter :: IJTEND = NXY

  ! Ocean depth settings
  double precision :: DEPTH(NZDIM)
  integer :: NBOT(NXYDIM)
  double precision :: DZ(NXYDIM,NZDIM)

  ! Initialize DEPTH values if NZDIM is 44
  data DEPTH /   5.0D0,   10.0D0,   15.0D0,   20.0D0,   25.0D0, &
                30.0D0,   37.5D0,   45.0D0,   55.0D0,   65.0D0, &
                75.0D0,   85.0D0,  100.0D0,  120.0D0,  145.0D0, &
               170.0D0,  200.0D0,  250.0D0,  325.0D0,  400.0D0, &
               500.0D0,  600.0D0,  750.0D0,  900.0D0, 1050.0D0, &
              1200.0D0, 1400.0D0, 1650.0D0, 1900.0D0, 2150.0D0, &
              2400.0D0, 2650.0D0, 2900.0D0, 3150.0D0, 3400.0D0, &
              3650.0D0, 3900.0D0, 4150.0D0, 4400.0D0, 4650.0D0, &
              4900.0D0, 5150.0D0, 5400.0D0, 5500.0D0 /

  ! Missing value
  double precision, parameter :: rmiss = 999.00D0
 !double precision, parameter :: rmiss = 0.0000D0

contains

  subroutine allocate_var1d(iidim, var1d)
    integer, intent(in) :: iidim
    double precision, allocatable, intent(out) :: var1d(:)
    allocate(var1d(iidim))
  end subroutine allocate_var1d

  subroutine allocate_var2d(iidim, jjdim, var2d)
    integer, intent(in) :: iidim, jjdim
    double precision, allocatable, intent(out) :: var2d(:,:)
    allocate(var2d(iidim, jjdim))
  end subroutine allocate_var2d

  subroutine allocate_var3d(iidim, jjdim, kkdim, var3d)
    integer, intent(in) :: iidim, jjdim, kkdim
    double precision, allocatable, intent(out) :: var3d(:,:,:)
    allocate(var3d(iidim, jjdim, kkdim))
  end subroutine allocate_var3d

end module module_digda_grid
