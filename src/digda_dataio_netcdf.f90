program digda_dataio_netcdf
  use netcdf_utilities
  use module_digda_main
  implicit none

  character(len=100) :: filename, varname, varunit

  integer :: i, j, k, l

  double precision, dimension(:,:), allocatable :: bdep    !! ocean depth [m]
  double precision, dimension(nxdim, nydim) :: omask   !! ocean depth [m]

  double precision, dimension(nxydim, benfdm) :: bfrc
  double precision, dimension(nxydim, sdzdim, sedtdm) :: sdtx

  double precision, dimension(:), allocatable :: lon_values
  double precision, dimension(:), allocatable :: lat_values
  double precision, dimension(:), allocatable :: depth_values

  double precision, dimension(:,:), allocatable :: vard2d
  double precision, dimension(:,:,:,:), allocatable :: vard3d

  ! Initialize MPI
  if (use_mpi) then
    call initialize_digda_mpi
  end if

  ! Call the subroutine to read the namelist
  call read_namelist()

  write(*, *) "'-------------------------------------------'"
  write(*, *) "'-- Start reading the boundary conditions --'"
  write(*, *) "'-------------------------------------------'"
  write(*, *) "    Tracer no.01 = SCLA"
  write(*, *) "    Tracer no.02 = SPOM"
  write(*, *) "    Tracer no.03 = SCAR"
  write(*, *) "    Tracer no.04 = SOPA"
  write(*, *) "    Tracer no.05 = SAGE"
  write(*, *) "    Tracer no.06 = CNIT"
  write(*, *) "    Tracer no.07 = CDIC"
  write(*, *) "    Tracer no.08 = CALK"
  write(*, *) "    Tracer no.09 = CSIL"
  write(*, *) "    Tracer no.10 = CDO2"
  write(*, *) "    Tracer no.11 = S13POM"
  write(*, *) "    Tracer no.12 = S13CAR"
  write(*, *) "    Tracer no.13 = C13DIC"
  write(*, *) "    Tracer no.14 = S14POM"
  write(*, *) "    Tracer no.15 = S14CAR"
  write(*, *) "    Tracer no.16 = C14DIC"
  write(*, *) "'-------------------------------------------'"

  ! Allocate the array with appropriate dimensions
  allocate(bdep(nxdim, nydim))
  allocate(vard2d(nxdim, nydim))

  ! Call the read_netcdf subroutine
  filename = trim(bound_dir)//"/bdep.nc"
  varname = "depth"
  call read_netcdf_2d(filename, varname, bdep, lon_values, lat_values)

  do j = 1, nydim
    do i = 1, nxdim
      if (bdep(i, j) == rmiss) then
        omask(i, j) = 0.D0
      else
        omask(i, j) = 1.D0
      end if
    end do
  end do

  filename = trim(bound_dir)//"/btem.nc"
  varname = "TEM"
  call read_netcdf_2d(filename, varname, vard2d, lon_values, lat_values)
  do j = 1, nydim
    do i = 1, nxdim
      l = i + (j - 1) * nxdim
      bfrc(l, 6) = vard2d(i, j)
    end do
  end do

  filename = trim(bound_dir)//"/bsal.nc"
  varname = "SAL"
  call read_netcdf_2d(filename, varname, vard2d, lon_values, lat_values)
  do j = 1, nydim
    do i = 1, nxdim
      l = i + (j - 1) * nxdim
      bfrc(l, 7) = vard2d(i, j)
    end do
  end do

  filename = trim(bound_dir)//"/bnit.nc"
  varname = "NIT"
  call read_netcdf_2d(filename, varname, vard2d, lon_values, lat_values)
  do j = 1, nydim
    do i = 1, nxdim
      l = i + (j - 1) * nxdim
      bfrc(l, 8) = vard2d(i, j)
    end do
  end do

  filename = trim(bound_dir)//"/bdic.nc"
  varname = "DIC"
  call read_netcdf_2d(filename, varname, vard2d, lon_values, lat_values)
  do j = 1, nydim
    do i = 1, nxdim
      l = i + (j - 1) * nxdim
      bfrc(l, 9) = vard2d(i, j)
    end do
  end do

  filename = trim(bound_dir)//"/balk.nc"
  varname = "ALK"
  call read_netcdf_2d(filename, varname, vard2d, lon_values, lat_values)
  do j = 1, nydim
    do i = 1, nxdim
      l = i + (j - 1) * nxdim
      bfrc(l, 10) = vard2d(i, j)
    end do
  end do

  filename = trim(bound_dir)//"/bsil.nc"
  varname = "SIL"
  call read_netcdf_2d(filename, varname, vard2d, lon_values, lat_values)
  do j = 1, nydim
    do i = 1, nxdim
      l = i + (j - 1) * nxdim
      bfrc(l, 11) = vard2d(i, j)
    end do
  end do

  filename = trim(bound_dir)//"/bdo2.nc"
  varname = "DO2"
  call read_netcdf_2d(filename, varname, vard2d, lon_values, lat_values)
  do j = 1, nydim
    do i = 1, nxdim
      l = i + (j - 1) * nxdim
      bfrc(l, 12) = vard2d(i, j)
    end do
  end do

  filename = trim(bound_dir)//"/bc13.nc"
  varname = "C13"
  call read_netcdf_2d(filename, varname, vard2d, lon_values, lat_values)
  do j = 1, nydim
    do i = 1, nxdim
      l = i + (j - 1) * nxdim
      bfrc(l, 13) = vard2d(i, j)
    end do
  end do

  filename = trim(bound_dir)//"/bc14.nc"
  varname = "C14"
  call read_netcdf_2d(filename, varname, vard2d, lon_values, lat_values)
  do j = 1, nydim
    do i = 1, nxdim
      l = i + (j - 1) * nxdim
      bfrc(l, 14) = vard2d(i, j)
    end do
  end do

  filename = trim(bound_dir)//"/ficla.nc"
  varname = "FICLA"
  call read_netcdf_2d(filename, varname, vard2d, lon_values, lat_values)
  do j = 1, nydim
    do i = 1, nxdim
      l = i + (j - 1) * nxdim
      bfrc(l, 1) = vard2d(i, j) * omask(i, j)
    end do
  end do

  filename = trim(bound_dir)//"/fipom.nc"
  varname = "FIPOM"
  call read_netcdf_2d(filename, varname, vard2d, lon_values, lat_values)
  do j = 1, nydim
    do i = 1, nxdim
      l = i + (j - 1) * nxdim
      bfrc(l, 2) = vard2d(i, j) * omask(i, j)
    end do
  end do

  filename = trim(bound_dir)//"/ficar.nc"
  varname = "FICAR"
  call read_netcdf_2d(filename, varname, vard2d, lon_values, lat_values)
  do j = 1, nydim
    do i = 1, nxdim
      l = i + (j - 1) * nxdim
      bfrc(l, 3) = vard2d(i, j) * omask(i, j)
    end do
  end do

  filename = trim(bound_dir)//"/fiopa.nc"
  varname = "FIOPA"
  call read_netcdf_2d(filename, varname, vard2d, lon_values, lat_values)
  do j = 1, nydim
    do i = 1, nxdim
      l = i + (j - 1) * nxdim
      bfrc(l, 4) = vard2d(i, j) * omask(i, j)
    end do
  end do

  do j = 1, nydim
    do i = 1, nxdim
      l = i + (j - 1) * nxdim
      vard2d(i, j) = 0.D0
      bfrc(l, 5) = vard2d(i, j) * omask(i, j)
    end do
  end do

  ! Synchronize all processes
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  ! Free memory for variables
  if (allocated(vard2d)) then
    deallocate(vard2d)
  end if

  write(*, *) "'--------------------------------------------'"
  write(*, *) "'-- Finish reading the boundary conditions --'"
  write(*, *) "'--------------------------------------------'"

  call digda_main( &
 &     bdep, &
 &     bfrc, &
 &     sdtx )

  if ( myrank == iroot ) then

  ! Free memory for variables
  if (allocated(bdep)) then
    deallocate(bdep)
  end if

  write(*, *) "'-------------------------------'"
  write(*, *) "'-- Start writing NetCDF file --'"
  write(*, *) "'-------------------------------'"

  ! Allocate the array with appropriate dimensions
  allocate(vard3d(nxdim, nydim, sdzdim, 1))

! do k = 1, sdzdim
!   do j = 1, nydim
!     do i = 1, nxdim
!       l = i + (j - 1) * nxdim
!       vard3d(i, j, k, 1) = sdtx(l, k, 1)
!     end do
!   end do
! end do
!   filename = "scla.nc"
!   varname = "Clay"
!   varunit = "g m-3"
!   call write_netcdf_3d(filename, varname, varunit, vard3d, lon_values, lat_values, lz, 999)

! do k = 1, sdzdim
!   do j = 1, nydim
!     do i = 1, nxdim
!       l = i + (j - 1) * nxdim
!       vard3d(i, j, k, 1) = sdtx(l, k, 2)
!     end do
!   end do
! end do
!   filename = "spom.nc"
!   varname = "OrgM"
!   varunit = "mmojg m-3"
!   call write_netcdf_3d(filename, varname, varunit, vard3d, lon_values, lat_values, lz, 999)

! do k = 1, sdzdim
!   do j = 1, nydim
!     do i = 1, nxdim
!       l = i + (j - 1) * nxdim
!       vard3d(i, j, k, 1) = sdtx(l, k, 3)
!     end do
!   end do
! end do
!   filename = "scar.nc"
!   varname = "CaCO3"
!   varunit = "mmol m-3"
!   call write_netcdf_3d(filename, varname, varunit, vard3d, lon_values, lat_values, lz, 999)

! do k = 1, sdzdim
!   do j = 1, nydim
!     do i = 1, nxdim
!       l = i + (j - 1) * nxdim
!       vard3d(i, j, k, 1) = sdtx(l, k, 4)
!     end do
!   end do
! end do
!   filename = "sopa.nc"
!   varname = "SiO2"
!   varunit = "mmol m-3"
!   call write_netcdf_3d(filename, varname, varunit, vard3d, lon_values, lat_values, lz, 999)

! do k = 1, sdzdim
!   do j = 1, nydim
!     do i = 1, nxdim
!       l = i + (j - 1) * nxdim
!       vard3d(i, j, k, 1) = sdtx(l, k, 5)
!     end do
!   end do
! end do
!   filename = "sage.nc"
!   varname = "Age"
!   varunit = "yr"
!   call write_netcdf_3d(filename, varname, varunit, vard3d, lon_values, lat_values, lz, 999)

! do k = 1, sdzdim
!   do j = 1, nydim
!     do i = 1, nxdim
!       l = i + (j - 1) * nxdim
!       vard3d(i, j, k, 1) = sdtx(l, k, 6)
!     end do
!   end do
! end do
!   filename = "cnit.nc"
!   varname = "Phosphate"
!   varunit = "mmol m-3"
!   call write_netcdf_3d(filename, varname, varunit, vard3d, lon_values, lat_values, lz, 999)

! do k = 1, sdzdim
!   do j = 1, nydim
!     do i = 1, nxdim
!       l = i + (j - 1) * nxdim
!       vard3d(i, j, k, 1) = sdtx(l, k, 4)
!     end do
!   end do
! end do
!   filename = "cdic.nc"
!   varname = "DIC"
!   varunit = "mmol m-3"
!   call write_netcdf_3d(filename, varname, varunit, vard3d, lon_values, lat_values, lz, 999)

! do k = 1, sdzdim
!   do j = 1, nydim
!     do i = 1, nxdim
!       l = i + (j - 1) * nxdim
!       vard3d(i, j, k, 1) = sdtx(l, k, 4)
!     end do
!   end do
! end do
!   filename = "calk.nc"
!   varname = "Alkalinity"
!   varunit = "mmol m-3"
!   call write_netcdf_3d(filename, varname, varunit, vard3d, lon_values, lat_values, lz, 999)

! do k = 1, sdzdim
!   do j = 1, nydim
!     do i = 1, nxdim
!       l = i + (j - 1) * nxdim
!       vard3d(i, j, k, 1) = sdtx(l, k, 4)
!     end do
!   end do
! end do
!   filename = "csil.nc"
!   varname = "Silicate"
!   varunit = "mmol m-3"
!   call write_netcdf_3d(filename, varname, varunit, vard3d, lon_values, lat_values, lz, 999)

! do k = 1, sdzdim
!   do j = 1, nydim
!     do i = 1, nxdim
!       l = i + (j - 1) * nxdim
!       vard3d(i, j, k, 1) = sdtx(l, k, 4)
!     end do
!   end do
! end do
!   filename = "cdo2.nc"
!   varname = "Oxygen"
!   varunit = "mmol m-3"
!   call write_netcdf_3d(filename, varname, varunit, vard3d, lon_values, lat_values, lz, 999)


  ! Free memory for variables
  if (allocated(vard3d)) then
    deallocate(vard3d)
  end if

  write(*, *) "'--------------------------------'"
  write(*, *) "'-- Finish writing NetCDF file --'"
  write(*, *) "'--------------------------------'"

  end if

  ! Finalize MPI
  call MPI_FINALIZE(ierr)

contains

  subroutine read_namelist()
    ! Open the namelist file
    integer :: unit, iostat
    character(len=100) :: namelist_file

    namelist_file = 'PARAMET'
    open(unit=10, file=namelist_file, status='old', iostat=iostat)

    if (iostat /= 0) then
      print *, 'Error opening namelist file: ', trim(namelist_file)
    else
      read(10, nml=nmocbc, iostat=iostat)
      if (iostat /= 0) then
        print *, 'Error reading namelist file: ', trim(namelist_file)
      else
      !!print *, 'Namelist read successfully: ', trim(bound_dir)
      end if
      close(10)
    end if
  end subroutine read_namelist

end program digda_dataio_netcdf

