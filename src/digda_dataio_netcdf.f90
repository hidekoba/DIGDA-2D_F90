module module_digda_dataio_netcdf
  use netcdf_utilities
  use module_zcfile
  use module_digda_grid
  use module_digda_mpi
  implicit none

  character(len=128) :: filename
  character(len=16) :: varname, varunit
  integer :: i, j, k, t, ij

  double precision, dimension(:),       allocatable :: longitude
  double precision, dimension(:),       allocatable :: latitude
  double precision, dimension(:),       allocatable :: depth
  integer(kind=4),  dimension(:),       allocatable :: time

  double precision, dimension(:,:),     allocatable :: vard2d
  double precision, dimension(:,:,:),   allocatable :: vard3d
  double precision, dimension(:,:,:,:), allocatable :: vard4d

  logical :: get_axis

contains

  subroutine digda_read_restart(sdtx, lon_values, lat_values, depth_values, time_values)
    double precision, dimension(:), allocatable, intent(out) :: lon_values
    double precision, dimension(:), allocatable, intent(out) :: lat_values
    double precision, dimension(:), allocatable, intent(out) :: depth_values
    integer(kind=4),  dimension(:), allocatable, intent(out) :: time_values
    double precision, dimension(nxydim, sdzdim, sedtdm), intent(out) :: sdtx

    integer :: ios

    integer :: n

    write(nfstdo, *) "'--------------------------------------------------'"
    write(nfstdo, *) "'-- DIGDA DATA I/O:                              --'"
    write(nfstdo, *) "'-- Started reading Sediment RESTART file        --'"
    write(nfstdo, *) "'--------------------------------------------------'"
    rewind(nfparm)
    read(nfparm, NML=nmocrt, IOSTAT=ios)
    write(nfstdo, nmocrt)

    filename = trim(init_dir) // "/" // trim(rsfilnm)
    write(nfstdo, *) "'--------------------------------------------------'"
    write(nfstdo, *) "'-- Read initial condition file: ", filename
    write(nfstdo, *) "'--------------------------------------------------'"

    do n = 1, sedtdm
      varname = dtrcnm(n)
      call read_netcdf_auto(filename, varname, vard4d, lon_values, lat_values, depth_values, time_values)
      do k = 1, sdzdim
        do j = 1, nydim
          do i = 1, nxdim
            ij = i + (j - 1) * nxdim
            sdtx(ij, k, n) = vard4d(i, j, k, 1)
          end do
        end do
      end do
      if ( n == 1 .and. myrank == iroot ) then
        write(*, *) "'-- DIGDA DATA I/O:                              --'"
        write(*, *) "'--   Grid information (Longitude)               --'"
        do i = 1, nxdim
           write(*, *), i, lon_values(i)
        end do
        write(*, *) "'--   Grid information (Latitude)                --'"
        do j = 1, nydim
           write(*, *), j, lat_values(j)
        end do
        write(*, *) "'--   Grid information (Depth)                   --'"
        do k = 1, sdzdim
           write(*, *), k, depth_values(k)
        end do
        write(*, *) "'--   Grid information (Time)                    --'"
           write(*, *), 1, time_values
      end if
      write(nfstdo, '(A,A)') "                  variable name: ", trim(varname)
    end do

   !!-- Read: Clay
   !varname = "Clay"
   !call read_netcdf_auto(filename, varname, vard4d, lon_values, lat_values, depth_values, time_values)
   !if ( myrank == iroot ) then
   !write(nfstdo, *) "'-- DIGDA DATA I/O:                              --'"
   !write(nfstdo, *) "'--   Grid information (Longitude)               --'"
   !do i = 1, nxdim
   !   write(nfstdo, *), i, lon_values(i)
   !end do
   !write(nfstdo, *) "'--   Grid information (Latitude)                --'"
   !do j = 1, nydim
   !   write(nfstdo, *), j, lat_values(j)
   !end do
   !write(nfstdo, *) "'--   Grid information (Depth)                   --'"
   !do k = 1, sdzdim
   !   write(nfstdo, *), k, depth_values(k)
   !end do
   !write(nfstdo, *) "'--   Grid information (Time)                    --'"
   !   write(nfstdo, *), 1, time_values
   !do k = 1, sdzdim
   !  do j = 1, nydim
   !    do i = 1, nxdim
   !      ij = i + (j - 1) * nxdim
   !      sdtx(ij, k, 1) = vard4d(i, j, k, 1)
   !    end do
   !  end do
   !end do

   !!-- Read: OrgM
   !varname = "OrgM"
   !call read_netcdf_auto(filename, varname, vard4d, lon_values, lat_values, depth_values, time_values)
   !do k = 1, sdzdim
   !  do j = 1, nydim
   !    do i = 1, nxdim
   !      ij = i + (j - 1) * nxdim
   !      sdtx(ij, k, 2) = vard4d(i, j, k, 1)
   !    end do
   !  end do
   !end do

   !!-- Read: CaCO3
   !varname = "CaCO3"
   !call read_netcdf_auto(filename, varname, vard4d, lon_values, lat_values, depth_values, time_values)
   !do k = 1, sdzdim
   !  do j = 1, nydim
   !    do i = 1, nxdim
   !      ij = i + (j - 1) * nxdim
   !      sdtx(ij, k, 3) = vard4d(i, j, k, 1)
   !    end do
   !  end do
   !end do

   !!-- Read: SiO2
   !varname = "SiO2"
   !call read_netcdf_auto(filename, varname, vard4d, lon_values, lat_values, depth_values, time_values)
   !do k = 1, sdzdim
   !  do j = 1, nydim
   !    do i = 1, nxdim
   !      ij = i + (j - 1) * nxdim
   !      sdtx(ij, k, 4) = vard4d(i, j, k, 1)
   !    end do
   !  end do
   !end do

   !!-- Read: Agetracer
   !varname = "Agetracer"
   !call read_netcdf_auto(filename, varname, vard4d, lon_values, lat_values, depth_values, time_values)
   !do k = 1, sdzdim
   !  do j = 1, nydim
   !    do i = 1, nxdim
   !      ij = i + (j - 1) * nxdim
   !      sdtx(ij, k, 5) = vard4d(i, j, k, 1)
   !    end do
   !  end do
   !end do

   !!-- Read: Phosphate
   !varname = "Phosphate"
   !call read_netcdf_auto(filename, varname, vard4d, lon_values, lat_values, depth_values, time_values)
   !do k = 1, sdzdim
   !  do j = 1, nydim
   !    do i = 1, nxdim
   !      ij = i + (j - 1) * nxdim
   !      sdtx(ij, k, 6) = vard4d(i, j, k, 1)
   !    end do
   !  end do
   !end do

   !!-- Read: DIC
   !varname = "DIC"
   !call read_netcdf_auto(filename, varname, vard4d, lon_values, lat_values, depth_values, time_values)
   !do k = 1, sdzdim
   !  do j = 1, nydim
   !    do i = 1, nxdim
   !      ij = i + (j - 1) * nxdim
   !      sdtx(ij, k, 7) = vard4d(i, j, k, 1)
   !    end do
   !  end do
   !end do

   !!-- Read: Alkalinity
   !varname = "Alkalinity"
   !call read_netcdf_auto(filename, varname, vard4d, lon_values, lat_values, depth_values, time_values)
   !do k = 1, sdzdim
   !  do j = 1, nydim
   !    do i = 1, nxdim
   !      ij = i + (j - 1) * nxdim
   !      sdtx(ij, k, 8) = vard4d(i, j, k, 1)
   !    end do
   !  end do
   !end do

   !!-- Read: Silicate
   !varname = "Silicate"
   !call read_netcdf_auto(filename, varname, vard4d, lon_values, lat_values, depth_values, time_values)
   !do k = 1, sdzdim
   !  do j = 1, nydim
   !    do i = 1, nxdim
   !      ij = i + (j - 1) * nxdim
   !      sdtx(ij, k, 9) = vard4d(i, j, k, 1)
   !    end do
   !  end do
   !end do

   !!-- Read: Oxygen
   !varname = "Oxygen"
   !call read_netcdf_auto(filename, varname, vard4d, lon_values, lat_values, depth_values, time_values)
   !do k = 1, sdzdim
   !  do j = 1, nydim
   !    do i = 1, nxdim
   !      ij = i + (j - 1) * nxdim
   !      sdtx(ij, k, 10) = vard4d(i, j, k, 1)
   !    end do
   !  end do
   !end do

    !-- Synchronize all processes --
    call mpi_barrier(mpi_comm_world, ierr)

    !-- Free memory for variables --
    if (allocated(vard3d)) deallocate(vard3d)
    if (allocated(vard4d)) deallocate(vard4d)

    !-- Save the grid information as module variables --
    allocate(longitude(size(lon_values)))
    allocate(latitude(size(lat_values)))
    allocate(depth(size(depth_values)))
    allocate(time(size(time_values)))
    longitude = lon_values
    latitude  = lat_values
    depth     = depth_values
    time      = time_values

    write(nfstdo, *) "'--------------------------------------------------'"
    write(nfstdo, *) "'-- DIGDA DATA I/O:                              --'"
    write(nfstdo, *) "'-- Finised reading Sediment RESTART file        --'"
    write(nfstdo, *) "'--------------------------------------------------'"
    call flush(nfstdo)

  end subroutine digda_read_restart

  subroutine digda_write_restart(sdtx, lon_values, lat_values, depth_values, time_values)

    double precision, dimension(nxydim, sdzdim, sedtdm), intent(in) :: sdtx
    double precision, dimension(:), allocatable, intent(in) :: lon_values
    double precision, dimension(:), allocatable, intent(in) :: lat_values
    double precision, dimension(:), allocatable, intent(in) :: depth_values
    integer(kind=4), dimension(:), allocatable, intent(in) :: time_values

    character(len=100)            :: filenames(sedtdm)

    integer :: ios

    integer :: n
    character(len=2) :: nstr

    !-- Allocate the array with appropriate dimensions --
    allocate(vard4d(nxdim, nydim, sdzdim, 1))

    write(nfstdo, *) "'--------------------------------------------------'"
    write(nfstdo, *) "'-- DIGDA DATA I/O:                              --'"
    write(nfstdo, *) "'-- Started writing Sediment RESTART file        --'"
    write(nfstdo, *) "'--------------------------------------------------'"

    filename = trim(rstt_dir)//"/SEDTRR.nc"
    write(nfstdo, *) "'--------------------------------------------------'"
    write(nfstdo, *) "'-- DIGDA DATA I/O:                              --'"
    write(nfstdo, *) "'-- Write restart file: ", filename
    write(nfstdo, *) "'--------------------------------------------------'"

    if ( myrank == iroot ) then

   !do n = 1, sedtdm
   !  write(nstr, '(I2.2)') n
   !  filenames(n) = trim(rstt_dir) // "/SEDTRR" // nstr // ".nc"

   !  write(nfstdo, *) "'--------------------------------------------------'"
   !  write(nfstdo, *) "'-- Write restart file: ", filenames(n)
   !  write(nfstdo, *) "'--------------------------------------------------'"
   !  varname = dtrcnm(n)
   !  varunit = dtrcun(n)
   !  do k = 1, sdzdim
   !    do j = 1, nydim
   !      do i = 1, nxdim
   !        ij = i + (j - 1) * nxdim
   !        vard4d(i, j, k, 1) = sdtx(ij, k, n)
   !      end do
   !    end do
   !  end do
   !  call write_netcdf_3dwt(filenames(n), varname, varunit, vard4d, longitude, latitude, depth, time)
   !end do

    !-- Write: Clay
    filenames(1) = trim(rstt_dir)//"/SEDTRR01.nc"
    write(nfstdo, *) "'--------------------------------------------------'"
    write(nfstdo, *) "  write restart file: ", filenames(1)
    write(nfstdo, *) "'--------------------------------------------------'"
    varname = "Clay"
    varunit = "mmol m-3"
    do k = 1, sdzdim
      do j = 1, nydim
        do i = 1, nxdim
          ij = i + (j - 1) * nxdim
          vard4d(i, j, k, 1) = sdtx(ij, k, 1)
        end do
      end do
    end do
    call write_netcdf_3dwt(filenames(1), varname, varunit, vard4d, longitude, latitude, depth, time)

    !-- Write: OrgM
    filenames(2) = trim(rstt_dir)//"/SEDTRR02.nc"
    write(nfstdo, *) "'--------------------------------------------------'"
    write(nfstdo, *) "  write restart file: ", filenames(2)
    write(nfstdo, *) "'--------------------------------------------------'"
    varname = "OrgM"
    varunit = "mmol m-3"
    do k = 1, sdzdim
      do j = 1, nydim
        do i = 1, nxdim
          ij = i + (j - 1) * nxdim
          vard4d(i, j, k, 1) = sdtx(ij, k, 2)
        end do
      end do
    end do
    call write_netcdf_3dwt(filenames(2), varname, varunit, vard4d, longitude, latitude, depth, time)

    filenames(3) = trim(rstt_dir)//"/SEDTRR03.nc"
    write(nfstdo, *) "'--------------------------------------------------'"
    write(nfstdo, *) "  write restart file: ", filenames(3)
    write(nfstdo, *) "'--------------------------------------------------'"
    varname = "CaCO3"
    varunit = "mmol m-3"
    do k = 1, sdzdim
      do j = 1, nydim
        do i = 1, nxdim
          ij = i + (j - 1) * nxdim
          vard4d(i, j, k, 1) = sdtx(ij, k, 3)
        end do
      end do
    end do
    call write_netcdf_3dwt(filenames(3), varname, varunit, vard4d, longitude, latitude, depth, time)

    filenames(4) = trim(rstt_dir)//"/SEDTRR04.nc"
    write(nfstdo, *) "'--------------------------------------------------'"
    write(nfstdo, *) "  write restart file: ", filenames(4)
    write(nfstdo, *) "'--------------------------------------------------'"
    varname = "SiO2"
    varunit = "mmol m-3"
    do k = 1, sdzdim
      do j = 1, nydim
        do i = 1, nxdim
          ij = i + (j - 1) * nxdim
          vard4d(i, j, k, 1) = sdtx(ij, k, 4)
        end do
      end do
    end do
    call write_netcdf_3dwt(filenames(4), varname, varunit, vard4d, longitude, latitude, depth, time)

    filenames(5) = trim(rstt_dir)//"/SEDTRR05.nc"
    write(nfstdo, *) "'--------------------------------------------------'"
    write(nfstdo, *) "  write restart file: ", filenames(5)
    write(nfstdo, *) "'--------------------------------------------------'"
    varname = "Agetracer"
    varunit = "mmol m-3"
    do k = 1, sdzdim
      do j = 1, nydim
        do i = 1, nxdim
          ij = i + (j - 1) * nxdim
          vard4d(i, j, k, 1) = sdtx(ij, k, 5)
        end do
      end do
    end do
    call write_netcdf_3dwt(filenames(5), varname, varunit, vard4d, longitude, latitude, depth, time)

    filenames(6) = trim(rstt_dir)//"/SEDTRR06.nc"
    write(nfstdo, *) "'--------------------------------------------------'"
    write(nfstdo, *) "  write restart file: ", filenames(6)
    write(nfstdo, *) "'--------------------------------------------------'"
    varname = "Phosphate"
    varunit = "mmol m-3"
    do k = 1, sdzdim
      do j = 1, nydim
        do i = 1, nxdim
          ij = i + (j - 1) * nxdim
          vard4d(i, j, k, 1) = sdtx(ij, k, 6)
        end do
      end do
    end do
    call write_netcdf_3dwt(filenames(6), varname, varunit, vard4d, longitude, latitude, depth, time)

    filenames(7) = trim(rstt_dir)//"/SEDTRR07.nc"
    write(nfstdo, *) "'--------------------------------------------------'"
    write(nfstdo, *) "  write restart file: ", filenames(7)
    write(nfstdo, *) "'--------------------------------------------------'"
    varname = "DIC"
    varunit = "mmol m-3"
    do k = 1, sdzdim
      do j = 1, nydim
        do i = 1, nxdim
          ij = i + (j - 1) * nxdim
          vard4d(i, j, k, 1) = sdtx(ij, k, 7)
        end do
      end do
    end do
    call write_netcdf_3dwt(filenames(7), varname, varunit, vard4d, longitude, latitude, depth, time)

    filenames(8) = trim(rstt_dir)//"/SEDTRR08.nc"
    write(nfstdo, *) "'--------------------------------------------------'"
    write(nfstdo, *) "  write restart file: ", filenames(8)
    write(nfstdo, *) "'--------------------------------------------------'"
    varname = "Alkalinity"
    varunit = "mmol m-3"
    do k = 1, sdzdim
      do j = 1, nydim
        do i = 1, nxdim
          ij = i + (j - 1) * nxdim
          vard4d(i, j, k, 1) = sdtx(ij, k, 8)
        end do
      end do
    end do
    call write_netcdf_3dwt(filenames(8), varname, varunit, vard4d, longitude, latitude, depth, time)

    filenames(9) = trim(rstt_dir)//"/SEDTRR09.nc"
    write(nfstdo, *) "'--------------------------------------------------'"
    write(nfstdo, *) "  write restart file: ", filenames(9)
    write(nfstdo, *) "'--------------------------------------------------'"
    varname = "Silicate"
    varunit = "mmol m-3"
    do k = 1, sdzdim
      do j = 1, nydim
        do i = 1, nxdim
          ij = i + (j - 1) * nxdim
          vard4d(i, j, k, 1) = sdtx(ij, k, 9)
        end do
      end do
    end do
    call write_netcdf_3dwt(filenames(9), varname, varunit, vard4d, longitude, latitude, depth, time)

    filenames(10) = trim(rstt_dir)//"/SEDTRR10.nc"
    write(nfstdo, *) "'--------------------------------------------------'"
    write(nfstdo, *) "  write restart file: ", filenames(10)
    write(nfstdo, *) "'--------------------------------------------------'"
    varname = "Oxygen"
    varunit = "mmol m-3"
    do k = 1, sdzdim
      do j = 1, nydim
        do i = 1, nxdim
          ij = i + (j - 1) * nxdim
          vard4d(i, j, k, 1) = sdtx(ij, k, 10)
        end do
      end do
    end do
    call write_netcdf_3dwt(filenames(10), varname, varunit, vard4d, longitude, latitude, depth, time)

    end if

    !-- Synchronize all processes --
    call mpi_barrier(mpi_comm_world, ierr)

    !-- Free memory for variables --
    if (allocated(vard3d)) deallocate(vard3d)
    if (allocated(vard4d)) deallocate(vard4d)

    write(nfstdo, *) "'--------------------------------------------------'"
    write(nfstdo, *) "'-- DIGDA DATA I/O:                              --'"
    write(nfstdo, *) "'-- Finised writing Sediment RESTART file        --'"
    write(nfstdo, *) "'--------------------------------------------------'"

  end subroutine digda_write_restart

  subroutine digda_read_bottom_forcing(bdep, bfrc)
!!subroutine digda_read_bottom_forcing(bdep, bfrc, omask)

    double precision, dimension(nxydim), intent(out) :: bdep
    double precision, dimension(nxydim, benfdm), intent(out) :: bfrc
!!  double precision, dimension(nxydim), intent(out) :: omask

    double precision, dimension(nxydim) :: omask

    double precision, dimension(:), allocatable :: lon_values
    double precision, dimension(:), allocatable :: lat_values
    double precision, dimension(:), allocatable :: depth_values
    integer(kind=4),  dimension(:), allocatable :: time_values

    integer :: ios

    integer :: n

    write(nfstdo, *) "'--------------------------------------------------'"
    write(nfstdo, *) "'-- DIGDA DATA I/O:                              --'"
    write(nfstdo, *) "'-- Started reading Sed-Water boundary forcing   --'"
    write(nfstdo, *) "'--------------------------------------------------'"
    rewind(nfparm)
    read(nfparm, NML=nmocbc, IOSTAT=ios)
    write(nfstdo, nmocbc)

    !-- Allocate the array with appropriate dimensions --
    allocate(vard2d(nxdim, nydim))
    allocate(vard3d(nxdim, nydim, 1))

    !-- Read forcing: Ocean depth
    n = 0
      filename = trim(bound_dir) // "/" // trim(bcfilnm(n))
      varname = bfrcnm(n)
      write(nfstdo, *) "'--------------------------------------------------'"
      write(nfstdo, '(A,A)') "'-- Read forcing file: ", trim(filename)
      write(nfstdo, '(A,A)') "        variable name: ", trim(varname)
      write(nfstdo, *) "'--------------------------------------------------'"
      call read_netcdf_auto(filename, varname, vard4d, lon_values, lat_values, depth_values, time_values)
      do j = 1, nydim
        do i = 1, nxdim
          ij = i + (j - 1) * nxdim
          bdep(ij) = vard4d(i, j, 1, 1)
          if (vard4d(i, j, 1, 1) == rmiss .or. vard4d(i, j, 1, 1) == 0.d0) then
            omask(ij) = 0.d0
          else
            omask(ij) = 1.d0
          end if
        end do
      end do

    !-- Read forcing: Benthic particle flux
    do n = 1, 4
      filename = trim(bound_dir) // "/" // trim(bcfilnm(n))
      varname = bfrcnm(n)
      write(nfstdo, *) "'--------------------------------------------------'"
      write(nfstdo, '(A,A)') "'-- Read forcing file: ", trim(filename)
      write(nfstdo, '(A,A)') "        variable name: ", trim(varname)
      write(nfstdo, *) "'--------------------------------------------------'"
      call read_netcdf_auto(filename, varname, vard4d, lon_values, lat_values, depth_values, time_values)
      do j = 1, nydim
        do i = 1, nxdim
          ij = i + (j - 1) * nxdim
          bfrc(ij, n) = vard4d(i, j, 1, 1) * omask(ij)
        end do
      end do
    end do

    !-- Read forcing: Age tracer
    n = 5
      filename = "not-specified"
      varname = bfrcnm(n)
      write(nfstdo, *) "'--------------------------------------------------'"
      write(nfstdo, '(A,A)') "'-- Read forcing file: ", trim(filename)
      write(nfstdo, '(A,A)') "        variable name: ", trim(varname)
      write(nfstdo, *) "'--------------------------------------------------'"
      do j = 1, nydim
        do i = 1, nxdim
          ij = i + (j - 1) * nxdim
          vard4d(i, j, 1, 1) = 0.d0
          bfrc(ij, n) = vard4d(i, j, 1, 1) * omask(ij)
        end do
      end do

    !-- Read forcing: Bottom-water properties
    do n = 6, benfdm
      filename = trim(bound_dir) // "/" // trim(bcfilnm(n))
      varname = bfrcnm(n)
      write(nfstdo, *) "'--------------------------------------------------'"
      write(nfstdo, '(A,A)') "'-- Read forcing file: ", trim(filename)
      write(nfstdo, '(A,A)') "        variable name: ", trim(varname)
      write(nfstdo, *) "'--------------------------------------------------'"
      call read_netcdf_auto(filename, varname, vard4d, lon_values, lat_values, depth_values, time_values)
      do j = 1, nydim
        do i = 1, nxdim
          ij = i + (j - 1) * nxdim
          bfrc(ij, n) = vard4d(i, j, 1, 1)
        end do
      end do
    end do

   !!-- Ocean Depth
   !filename = trim(bound_dir)//"/bdep.nc"
   !varname = "DEPTH"
   !write(nfstdo, *) "'--------------------------------------------------'"
   !write(nfstdo, *) "'-- Read forcing file: ", filename
   !write(nfstdo, *) "                       ", varname
   !write(nfstdo, *) "'--------------------------------------------------'"
   !call read_netcdf_auto(filename, varname, vard4d, lon_values, lat_values, depth_values, time_values)
   !do j = 1, nydim
   !  do i = 1, nxdim
   !    ij = i + (j - 1) * nxdim
   !    bdep(ij) = vard4d(i, j, 1, 1)
   !    if (vard4d(i, j, 1, 1) == rmiss .or. vard4d(i, j, 1, 1) == 0.0) then
   !      omask(ij) = 0.d0
   !    else
   !      omask(ij) = 1.d0
   !    end if
   !  end do
   !end do

   !!-- Temperature
   !filename = trim(bound_dir)//"/btem.nc"
   !varname = "TEM"
   !write(nfstdo, *) "'--------------------------------------------------'"
   !write(nfstdo, *) "'-- Read forcing file: ", filename
   !write(nfstdo, *) "                       ", varname
   !write(nfstdo, *) "'--------------------------------------------------'"
   !call read_netcdf_auto(filename, varname, vard4d, lon_values, lat_values, depth_values, time_values)
   !do j = 1, nydim
   !  do i = 1, nxdim
   !    ij = i + (j - 1) * nxdim
   !    bfrc(ij, 6) = vard4d(i, j, 1, 1)
   !  end do
   !end do

   !!-- Salinity
   !filename = trim(bound_dir)//"/bsal.nc"
   !varname = "SAL"
   !write(nfstdo, *) "'--------------------------------------------------'"
   !write(nfstdo, *) "'-- Read forcing file: ", filename
   !write(nfstdo, *) "                       ", varname
   !write(nfstdo, *) "'--------------------------------------------------'"
   !call read_netcdf_auto(filename, varname, vard4d, lon_values, lat_values, depth_values, time_values)
   !do j = 1, nydim
   !  do i = 1, nxdim
   !    ij = i + (j - 1) * nxdim
   !    bfrc(ij, 7) = vard4d(i, j, 1, 1)
   !  end do
   !end do

   !!-- Phosphate
   !filename = trim(bound_dir)//"/bnit.nc"
   !varname = "NIT"
   !write(nfstdo, *) "'--------------------------------------------------'"
   !write(nfstdo, *) "'-- Read forcing file: ", filename
   !write(nfstdo, *) "                       ", varname
   !write(nfstdo, *) "'--------------------------------------------------'"
   !call read_netcdf_auto(filename, varname, vard4d, lon_values, lat_values, depth_values, time_values)
   !do j = 1, nydim
   !  do i = 1, nxdim
   !    ij = i + (j - 1) * nxdim
   !    bfrc(ij, 8) = vard4d(i, j, 1, 1)
   !  end do
   !end do

   !!-- DIC
   !filename = trim(bound_dir)//"/bdic.nc"
   !varname = "DIC"
   !write(nfstdo, *) "'--------------------------------------------------'"
   !write(nfstdo, *) "'-- Read forcing file: ", filename
   !write(nfstdo, *) "                       ", varname
   !write(nfstdo, *) "'--------------------------------------------------'"
   !call read_netcdf_auto(filename, varname, vard4d, lon_values, lat_values, depth_values, time_values)
   !do j = 1, nydim
   !  do i = 1, nxdim
   !    ij = i + (j - 1) * nxdim
   !    bfrc(ij, 9) = vard4d(i, j, 1, 1)
   !  end do
   !end do

   !!-- Alkalinity
   !filename = trim(bound_dir)//"/balk.nc"
   !varname = "ALK"
   !write(nfstdo, *) "'--------------------------------------------------'"
   !write(nfstdo, *) "'-- Read forcing file: ", filename
   !write(nfstdo, *) "                       ", varname
   !write(nfstdo, *) "'--------------------------------------------------'"
   !call read_netcdf_auto(filename, varname, vard4d, lon_values, lat_values, depth_values, time_values)
   !do j = 1, nydim
   !  do i = 1, nxdim
   !    ij = i + (j - 1) * nxdim
   !    bfrc(ij, 10) = vard4d(i, j, 1, 1)
   !  end do
   !end do

   !!-- Silicate
   !filename = trim(bound_dir)//"/bsil.nc"
   !varname = "SIL"
   !write(nfstdo, *) "'--------------------------------------------------'"
   !write(nfstdo, *) "'-- Read forcing file: ", filename
   !write(nfstdo, *) "                       ", varname
   !write(nfstdo, *) "'--------------------------------------------------'"
   !call read_netcdf_auto(filename, varname, vard4d, lon_values, lat_values, depth_values, time_values)
   !do j = 1, nydim
   !  do i = 1, nxdim
   !    ij = i + (j - 1) * nxdim
   !    bfrc(ij, 11) = vard4d(i, j, 1, 1)
   !  end do
   !end do

   !!-- Oxygen
   !filename = trim(bound_dir)//"/bdo2.nc"
   !varname = "DO2"
   !write(nfstdo, *) "'--------------------------------------------------'"
   !write(nfstdo, *) "'-- Read forcing file: ", filename
   !write(nfstdo, *) "                       ", varname
   !write(nfstdo, *) "'--------------------------------------------------'"
   !call read_netcdf_auto(filename, varname, vard4d, lon_values, lat_values, depth_values, time_values)
   !do j = 1, nydim
   !  do i = 1, nxdim
   !    ij = i + (j - 1) * nxdim
   !    bfrc(ij, 12) = vard4d(i, j, 1, 1)
   !  end do
   !end do

   !!-- 13C-DIC
   !filename = trim(bound_dir)//"/bc13.nc"
   !varname = "C13"
   !write(nfstdo, *) "'--------------------------------------------------'"
   !write(nfstdo, *) "'-- Read forcing file: ", filename
   !write(nfstdo, *) "                       ", varname
   !write(nfstdo, *) "'--------------------------------------------------'"
   !call read_netcdf_auto(filename, varname, vard4d, lon_values, lat_values, depth_values, time_values)
   !do j = 1, nydim
   !  do i = 1, nxdim
   !    ij = i + (j - 1) * nxdim
   !    bfrc(ij, 13) = vard4d(i, j, 1, 1)
   !  end do
   !end do

   !!-- 14C-DIC
   !filename = trim(bound_dir)//"/bc14.nc"
   !varname = "C14"
   !write(nfstdo, *) "'--------------------------------------------------'"
   !write(nfstdo, *) "'-- Read forcing file: ", filename
   !write(nfstdo, *) "                       ", varname
   !write(nfstdo, *) "'--------------------------------------------------'"
   !call read_netcdf_auto(filename, varname, vard4d, lon_values, lat_values, depth_values, time_values)
   !do j = 1, nydim
   !  do i = 1, nxdim
   !    ij = i + (j - 1) * nxdim
   !    bfrc(ij, 14) = vard4d(i, j, 1, 1)
   !  end do
   !end do

   !!-- Clay
   !filename = trim(bound_dir)//"/ficla.nc"
   !varname = "FICLA"
   !write(nfstdo, *) "'--------------------------------------------------'"
   !write(nfstdo, *) "'-- Read forcing file: ", filename
   !write(nfstdo, *) "                       ", varname
   !write(nfstdo, *) "'--------------------------------------------------'"
   !call read_netcdf_auto(filename, varname, vard4d, lon_values, lat_values, depth_values, time_values)
   !do j = 1, nydim
   !  do i = 1, nxdim
   !    ij = i + (j - 1) * nxdim
   !    bfrc(ij, 1) = vard4d(i, j, 1, 1) * omask(ij)
   !  end do
   !end do

   !!-- POM
   !filename = trim(bound_dir)//"/fipom.nc"
   !varname = "FIPOM"
   !write(nfstdo, *) "'--------------------------------------------------'"
   !write(nfstdo, *) "'-- Read forcing file: ", filename
   !write(nfstdo, *) "                       ", varname
   !write(nfstdo, *) "'--------------------------------------------------'"
   !call read_netcdf_auto(filename, varname, vard4d, lon_values, lat_values, depth_values, time_values)
   !do j = 1, nydim
   !  do i = 1, nxdim
   !    ij = i + (j - 1) * nxdim
   !    bfrc(ij, 2) = vard4d(i, j, 1, 1) * omask(ij)
   !  end do
   !end do

   !!-- CaCO3
   !filename = trim(bound_dir)//"/ficar.nc"
   !varname = "FICAR"
   !write(nfstdo, *) "'--------------------------------------------------'"
   !write(nfstdo, *) "'-- Read forcing file: ", filename
   !write(nfstdo, *) "                       ", varname
   !write(nfstdo, *) "'--------------------------------------------------'"
   !call read_netcdf_auto(filename, varname, vard4d, lon_values, lat_values, depth_values, time_values)
   !do j = 1, nydim
   !  do i = 1, nxdim
   !    ij = i + (j - 1) * nxdim
   !    bfrc(ij, 3) = vard4d(i, j, 1, 1) * omask(ij)
   !  end do
   !end do

   !!-- SiO2
   !filename = trim(bound_dir)//"/fiopa.nc"
   !varname = "FIOPA"
   !write(nfstdo, *) "'--------------------------------------------------'"
   !write(nfstdo, *) "'-- Read forcing file: ", filename
   !write(nfstdo, *) "                       ", varname
   !write(nfstdo, *) "'--------------------------------------------------'"
   !call read_netcdf_auto(filename, varname, vard4d, lon_values, lat_values, depth_values, time_values)
   !do j = 1, nydim
   !  do i = 1, nxdim
   !    ij = i + (j - 1) * nxdim
   !    bfrc(ij, 4) = vard4d(i, j, 1, 1) * omask(ij)
   !  end do
   !end do

   !!-- Age tracer
   !do j = 1, nydim
   !  do i = 1, nxdim
   !    ij = i + (j - 1) * nxdim
   !    vard4d(i, j, 1, 1) = 0.d0
   !    bfrc(ij, 5) = vard4d(i, j, 1, 1) * omask(ij)
   !  end do
   !end do

    !-- Synchronize all processes --
    call mpi_barrier(mpi_comm_world, ierr)

    !-- Free memory for variables --
    if (allocated(vard2d)) deallocate(vard2d)
    if (allocated(vard3d)) deallocate(vard3d)
    if (allocated(vard4d)) deallocate(vard4d)

    write(nfstdo, *) "'--------------------------------------------------'"
    write(nfstdo, *) "'-- DIGDA DATA I/O:                              --'"
    write(nfstdo, *) "'-- Finished reading Sed-Water boundary forcing  --'"
    write(nfstdo, *) "'--------------------------------------------------'"
    call flush(nfstdo)

  end subroutine digda_read_bottom_forcing

  subroutine digda_write_data()

    if ( myrank == iroot ) then

    write(nfstdo, *) "'--------------------------------------------------'"
    write(nfstdo, *) "'-- DIGDA DATA I/O:                              --'"
    write(nfstdo, *) "'-- Started writing NetCDF file                  --'"
    write(nfstdo, *) "'--------------------------------------------------'"

    !-- Allocate the array with appropriate dimensions --
    allocate(vard4d(nxdim, nydim, sdzdim, 1))

  ! do k = 1, sdzdim
  !   do j = 1, nydim
  !     do i = 1, nxdim
  !       ij = i + (j - 1) * nxdim
  !       vard3d(i, j, k, 1) = sdtx(l, k, 1)
  !     end do
  !   end do
  ! end do
  !   filename = "scla.nc"
  !   varname = "Clay"
  !   varunit = "g m-3"
  !   call write_netcdf_3dwt(filename, varname, varunit, vard3d, lon_values, lat_values, lz, 999)

  ! do k = 1, sdzdim
  !   do j = 1, nydim
  !     do i = 1, nxdim
  !       ij = i + (j - 1) * nxdim
  !       vard3d(i, j, k, 1) = sdtx(l, k, 2)
  !     end do
  !   end do
  ! end do
  !   filename = "spom.nc"
  !   varname = "OrgM"
  !   varunit = "mmojg m-3"
  !   call write_netcdf_3dwt(filename, varname, varunit, vard3d, lon_values, lat_values, lz, 999)

  ! do k = 1, sdzdim
  !   do j = 1, nydim
  !     do i = 1, nxdim
  !       ij = i + (j - 1) * nxdim
  !       vard3d(i, j, k, 1) = sdtx(l, k, 3)
  !     end do
  !   end do
  ! end do
  !   filename = "scar.nc"
  !   varname = "CaCO3"
  !   varunit = "mmol m-3"
  !   call write_netcdf_3dwt(filename, varname, varunit, vard3d, lon_values, lat_values, lz, 999)

  ! do k = 1, sdzdim
  !   do j = 1, nydim
  !     do i = 1, nxdim
  !       ij = i + (j - 1) * nxdim
  !       vard3d(i, j, k, 1) = sdtx(l, k, 4)
  !     end do
  !   end do
  ! end do
  !   filename = "sopa.nc"
  !   varname = "SiO2"
  !   varunit = "mmol m-3"
  !   call write_netcdf_3dwt(filename, varname, varunit, vard3d, lon_values, lat_values, lz, 999)

  ! do k = 1, sdzdim
  !   do j = 1, nydim
  !     do i = 1, nxdim
  !       ij = i + (j - 1) * nxdim
  !       vard3d(i, j, k, 1) = sdtx(l, k, 5)
  !     end do
  !   end do
  ! end do
  !   filename = "sage.nc"
  !   varname = "Age"
  !   varunit = "yr"
  !   call write_netcdf_3dwt(filename, varname, varunit, vard3d, lon_values, lat_values, lz, 999)

  ! do k = 1, sdzdim
  !   do j = 1, nydim
  !     do i = 1, nxdim
  !       ij = i + (j - 1) * nxdim
  !       vard3d(i, j, k, 1) = sdtx(l, k, 6)
  !     end do
  !   end do
  ! end do
  !   filename = "cnit.nc"
  !   varname = "Phosphate"
  !   varunit = "mmol m-3"
  !   call write_netcdf_3dwt(filename, varname, varunit, vard3d, lon_values, lat_values, lz, 999)

  ! do k = 1, sdzdim
  !   do j = 1, nydim
  !     do i = 1, nxdim
  !       ij = i + (j - 1) * nxdim
  !       vard3d(i, j, k, 1) = sdtx(l, k, 4)
  !     end do
  !   end do
  ! end do
  !   filename = "cdic.nc"
  !   varname = "DIC"
  !   varunit = "mmol m-3"
  !   call write_netcdf_3dwt(filename, varname, varunit, vard3d, lon_values, lat_values, lz, 999)

  ! do k = 1, sdzdim
  !   do j = 1, nydim
  !     do i = 1, nxdim
  !       ij = i + (j - 1) * nxdim
  !       vard3d(i, j, k, 1) = sdtx(l, k, 4)
  !     end do
  !   end do
  ! end do
  !   filename = "calk.nc"
  !   varname = "Alkalinity"
  !   varunit = "mmol m-3"
  !   call write_netcdf_3dwt(filename, varname, varunit, vard3d, lon_values, lat_values, lz, 999)

  ! do k = 1, sdzdim
  !   do j = 1, nydim
  !     do i = 1, nxdim
  !       ij = i + (j - 1) * nxdim
  !       vard3d(i, j, k, 1) = sdtx(l, k, 4)
  !     end do
  !   end do
  ! end do
  !   filename = "csil.nc"
  !   varname = "Silicate"
  !   varunit = "mmol m-3"
  !   call write_netcdf_3dwt(filename, varname, varunit, vard3d, lon_values, lat_values, lz, 999)

  ! do k = 1, sdzdim
  !   do j = 1, nydim
  !     do i = 1, nxdim
  !       ij = i + (j - 1) * nxdim
  !       vard3d(i, j, k, 1) = sdtx(l, k, 4)
  !     end do
  !   end do
  ! end do
  !   filename = "cdo2.nc"
  !   varname = "Oxygen"
  !   varunit = "mmol m-3"
  !   call write_netcdf_3dwt(filename, varname, varunit, vard3d, lon_values, lat_values, lz, 999)


    ! Free memory for variables
    if (allocated(vard3d)) deallocate(vard3d)
    if (allocated(vard4d)) deallocate(vard4d)

    write(nfstdo, *) "'--------------------------------------------------'"
    write(nfstdo, *) "'-- DIGDA DATA I/O:                              --'"
    write(nfstdo, *) "'-- Finished writing NetCDF file                 --'"
    write(nfstdo, *) "'--------------------------------------------------'"
    call flush(nfstdo)

    end if

  end subroutine digda_write_data

! subroutine read_namelist()
!   integer :: ios
!   character(len=64) :: namelist_file

!   namelist_file = "PARAMET"
!   open(UNIT=nfparm, FILE=namelist_file, STATUS="old", iostat=ios)

!   if (ios /= 0) then
!     write(*, *) "Error opening namelist file: ", trim(namelist_file)
!   else
!     read(nfparm, NML=nmocbc, IOSTAT=ios)
!     if (ios /= 0) then
!       write(*, *) "Error reading namelist file: ", trim(namelist_file)
!     else
!     !!write(*, *) "Namelist read successfully: ", trim(bound_dir)
!     end if
!   !!close(nfparm)
!   end if
! end subroutine read_namelist

end module module_digda_dataio_netcdf
