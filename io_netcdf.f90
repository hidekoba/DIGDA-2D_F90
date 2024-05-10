program io_netcdf
    use netcdf_utilities
!   use gtool_utilities
    implicit none

    character(len=100) :: filename
    character(len=100) :: varname
    integer(kind=4) :: lonsize, latsize, depsize
    integer(kind=4) :: nsize
    integer(kind=4) :: nx = 2560
    integer(kind=4) :: ny = 1920
    integer(kind=4) :: nz =  44
    character(len=16) :: citem
    character(len=16) :: ctitle
    character(len=16) :: cunit
    real(kind=4), dimension(:,:), allocatable :: bdep
    real(kind=4), dimension(:,:), allocatable :: btem, bsal, bdip, bdic, balk, bdo2, bsil
    real(kind=4), dimension(:,:), allocatable :: bc13, bc14
    real(kind=4), dimension(:,:), allocatable :: ficla, fipom, ficar, fiopa
    real(kind=8), dimension(:,:,:), allocatable :: depth
    real(kind=8), dimension(:,:), allocatable :: bdepth
    real(kind=8), dimension(:), allocatable :: data_flat
    integer :: i, j, k

    ! Call the read_netcdf subroutine
    filename = "bdep.nc"
    varname = "depth"
    call read_netcdf_2d(filename, varname, bdep)
    filename = "btem.nc"
    varname = "TEM"
    call read_netcdf_2d(filename, varname, btem)
    filename = "bsal.nc"
    varname = "SAL"
    call read_netcdf_2d(filename, varname, bsal)
    filename = "bdip.nc"
    varname = "DIP"
    call read_netcdf_2d(filename, varname, bdip)
    filename = "bdic.nc"
    varname = "DIC"
    call read_netcdf_2d(filename, varname, bdic)
    filename = "balk.nc"
    varname = "ALK"
    call read_netcdf_2d(filename, varname, balk)
    filename = "bdo2.nc"
    varname = "DO2"
    call read_netcdf_2d(filename, varname, bdo2)
    filename = "bsil.nc"
    varname = "SIL"
    call read_netcdf_2d(filename, varname, bsil)
    filename = "bc13.nc"
    varname = "C13"
    call read_netcdf_2d(filename, varname, bc13)
    filename = "bc14.nc"
    varname = "C14"
    call read_netcdf_2d(filename, varname, bc14)
    filename = "ficla.nc"
    varname = "FICLA"
    call read_netcdf_2d(filename, varname, ficla)
    filename = "fipom.nc"
    varname = "FIPOM"
    call read_netcdf_2d(filename, varname, fipom)
    filename = "ficar.nc"
    varname = "FICAR"
    call read_netcdf_2d(filename, varname, ficar)
    filename = "fiopa.nc"
    varname = "FIOPA"
    call read_netcdf_2d(filename, varname, fiopa)

!   ! Call the read_gtool subroutine
!   filename = "bdep.gt"
!   nsize = nx * ny
!   call read_gtool(filename, nsize, data_flat)

!    ! Now reshape data_flat into bdepth
!   if (allocated(data_flat)) then
!       allocate(bdepth(nx, ny))
!       bdepth = reshape(data_flat, shape(bdepth))
!       deallocate(data_flat)
!   else
!       write(*,*) "Data allocation failed for data_flat."
!       stop
!   endif

!   ! Process
!   if (allocated(depth)) then
!       write(*,*) "First element of ", varname, " = ", depth(1,1,1)
!   else
!       write(*,*) "Data allocation failed."
!       stop
!   endif

!   if (allocated(bdepth)) then
!       write(*,*) "First element of ", varname, " = ", bdepth(128,96)
!   else
!       write(*,*) "Data allocation failed."
!       stop
!   endif

!   do j = 1, ny
!       do i = 1, nx
!           bdepth(i,j) = -5. * bdepth(i,j)
!       end do
!   end do

    ! Call the write_netcdf subroutine
    filename = "output.nc"
    varname = "ALK"
    call write_netcdf_2d(filename, varname, balk)

!   ! Call the write_netcdf subroutine
!   filename = "output.gt"
!   citem = "dEpTh"
!   ctitle = "DePtH"
!   cunit = "m"
!   call write_gtool(filename, citem, ctitle, cunit, nx, ny, 1, bdepth)

!   ! Free memory for variables
!   if (allocated(depth)) deallocate(depth)
!   if (allocated(bdepth)) deallocate(bdepth)
!   if (allocated(data_flat)) deallocate(data_flat)

end program io_netcdf
