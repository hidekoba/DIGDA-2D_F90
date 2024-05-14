program io_netcdf
    use netcdf_utilities
    use sdmain_module
    implicit none

    character(len=100) :: filename
    character(len=100) :: varname
!   integer(kind=4) :: nsize
!   integer(kind=4) :: nx = 2560
!   integer(kind=4) :: ny = 1920
!   integer(kind=4) :: nz =  44
!   character(len=16) :: citem
!   character(len=16) :: ctitle
!   character(len=16) :: cunit

!   real(kind=4), dimension(:,:), allocatable :: bdep
!   real(kind=4), dimension(:,:), allocatable :: btem, bsal, bnit, bdic, balk, bdo2, bsil
!   real(kind=4), dimension(:,:), allocatable :: bc13, bc14
!   real(kind=4), dimension(:,:), allocatable :: ficla, fipom, ficar, fiopa
!   real(kind=8), dimension(:,:,:), allocatable :: depth
!   real(kind=8), dimension(:,:), allocatable :: bdepth
!   real(kind=8), dimension(:), allocatable :: data_flat
    integer :: tmpi, tmpj, tmpk

    double precision, dimension(:,:), allocatable :: inficla   !! benthic clay flux [g/m3*cm/s]
    double precision, dimension(:,:), allocatable :: infipom   !! benthic pom flux [mMol/m3*cm/s]
    double precision, dimension(:,:), allocatable :: inficar   !! benthic caco3 flux [mMol/m3*cm/s]
    double precision, dimension(:,:), allocatable :: infiopa   !! benthic opal flux [mMol/m3*cm/s]
    double precision, dimension(:,:), allocatable :: infiage   !! benthic opal flux [mMol/m3*cm/s]

    double precision, dimension(:,:), allocatable :: inbdep    !! ocean depth [m]
    double precision, dimension(:,:), allocatable :: inbwtem   !! temperature [mMol/m3]
    double precision, dimension(:,:), allocatable :: inbwsal   !! salinity [mMol/m3]
    double precision, dimension(:,:), allocatable :: inbwnit   !! phosphate [mMol/m3]
    double precision, dimension(:,:), allocatable :: inbwdic   !! DIC [mMol/m3]
    double precision, dimension(:,:), allocatable :: inbwalk   !! alkalinity [mMol/m3]
    double precision, dimension(:,:), allocatable :: inbwdo2   !! oxygen [mMol/m3]
    double precision, dimension(:,:), allocatable :: inbwsil   !! silicate [mMol/m3]o
    double precision, dimension(:,:), allocatable :: inbw13dic !! 13C-DIC [mMol/m3]
    double precision, dimension(:,:), allocatable :: inbw14dic !! 14C-DIC [mMol/m3]

    ! Allocate the array with appropriate dimensions
    allocate(inbdep(nxdim, nydim))
    allocate(inficla(nxdim, nydim))
    allocate(infipom(nxdim, nydim))
    allocate(inficar(nxdim, nydim))
    allocate(infiopa(nxdim, nydim))
    allocate(inbwtem(nxdim, nydim))
    allocate(inbwsal(nxdim, nydim))
    allocate(inbwnit(nxdim, nydim))
    allocate(inbwdic(nxdim, nydim))
    allocate(inbwalk(nxdim, nydim))
    allocate(inbwdo2(nxdim, nydim))
    allocate(inbwsil(nxdim, nydim))
    allocate(inbw13dic(nxdim, nydim))
    allocate(inbw14dic(nxdim, nydim))

    ! Call the read_netcdf subroutine
    filename = "bound256192/bdep.nc"
    varname = "depth"
    call read_netcdf_2d(filename, varname, inbdep)
    filename = "bound256192/btem.nc"
    varname = "TEM"
    call read_netcdf_2d(filename, varname, inbwtem)
    filename = "bound256192/bsal.nc"
    varname = "SAL"
    call read_netcdf_2d(filename, varname, inbwsal)
    filename = "bound256192/bnit.nc"
    varname = "NIT"
    call read_netcdf_2d(filename, varname, inbwnit)
    filename = "bound256192/bdic.nc"
    varname = "DIC"
    call read_netcdf_2d(filename, varname, inbwdic)
    filename = "bound256192/balk.nc"
    varname = "ALK"
    call read_netcdf_2d(filename, varname, inbwalk)
    filename = "bound256192/bdo2.nc"
    varname = "DO2"
    call read_netcdf_2d(filename, varname, inbwdo2)
    filename = "bound256192/bsil.nc"
    varname = "SIL"
    call read_netcdf_2d(filename, varname, inbwsil)
    filename = "bound256192/bc13.nc"
    varname = "C13"
    call read_netcdf_2d(filename, varname, inbw13dic)
    filename = "bound256192/bc14.nc"
    varname = "C14"
    call read_netcdf_2d(filename, varname, inbw14dic)
    filename = "bound256192/ficla.nc"
    varname = "FICLA"
    call read_netcdf_2d(filename, varname, inficla)
    filename = "bound256192/fipom.nc"
    varname = "FIPOM"
    call read_netcdf_2d(filename, varname, infipom)
    filename = "bound256192/ficar.nc"
    varname = "FICAR"
    call read_netcdf_2d(filename, varname, inficar)
    filename = "bound256192/fiopa.nc"
    varname = "FIOPA"
    call read_netcdf_2d(filename, varname, infiopa)

    infiage = 0.d0

!    call sd
     call sd( &
    &     inbdep, &
    &     inficla, infipom, inficar, infiopa, infiage, &
    &     inbwtem, inbwsal, inbwnit, inbwdic, inbwalk, inbwsil, inbwdo2, &
    &     inbw13dic, inbw14dic )
!!  &     inbw13dic, inbw14dic, &
!!  &     infopom, focar, foopa )
     do tmpj = 1, nydim
     do tmpi = 1, nxdim
     write(*,*) tmpi, tmpj, inbwtem(tmpi,tmpj)
     end do
     end do

!
!!   ! Call the read_gtool subroutine
!!   filename = "bound256192/bdep.gt"
!!   nsize = nx * ny
!!   call read_gtool(filename, nsize, data_flat)
!
!!    ! Now reshape data_flat into bdepth
!!   if (allocated(data_flat)) then
!!       allocate(bdepth(nx, ny))
!!       bdepth = reshape(data_flat, shape(bdepth))
!!       deallocate(data_flat)
!!   else
!!       write(*,*) "Data allocation failed for data_flat."
!!       stop
!!   endif
!
!!   ! Process
!!   if (allocated(depth)) then
!!       write(*,*) "First element of ", varname, " = ", depth(1,1,1)
!!   else
!!       write(*,*) "Data allocation failed."
!!       stop
!!   endif
!
!!   if (allocated(bdepth)) then
!!       write(*,*) "First element of ", varname, " = ", bdepth(128,96)
!!   else
!!       write(*,*) "Data allocation failed."
!!       stop
!!   endif
!
!!   do j = 1, ny
!!       do i = 1, nx
!!           bdepth(i,j) = -5. * bdepth(i,j)
!!       end do
!!   end do
!
!    ! Call the write_netcdf subroutine
!    filename = "bound256192/output.nc"
!    varname = "ALK"
!    call write_netcdf_2d(filename, varname, balk)
!
!!   ! Call the write_netcdf subroutine
!!   filename = "bound256192/output.gt"
!!   citem = "dEpTh"
!!   ctitle = "DePtH"
!!   cunit = "m"
!!   call write_gtool(filename, citem, ctitle, cunit, nx, ny, 1, bdepth)
!
!!   ! Free memory for variables
!!   if (allocated(depth)) deallocate(depth)
!!   if (allocated(bdepth)) deallocate(bdepth)
!!   if (allocated(data_flat)) deallocate(data_flat)
!
end program io_netcdf
