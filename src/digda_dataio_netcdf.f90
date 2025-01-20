module module_digda_dataio_netcdf
  use netcdf_utilities
  use module_zcfile
  use module_digda_grid
  use module_digda_mpi
  implicit none

  character(len=100) :: filename, varname, varunit
  integer :: i, j, k, t, ij

! double precision, dimension(nxydim, sdzdim, sedtdm) :: sdtx

  double precision, dimension(:), allocatable :: lon_values
  double precision, dimension(:), allocatable :: lat_values
  double precision, dimension(:), allocatable :: depth_values
  integer(kind=4), dimension(:), allocatable :: time_values

  double precision, dimension(:,:), allocatable :: vard2d
  double precision, dimension(:,:,:), allocatable :: vard3d
  double precision, dimension(:,:,:,:), allocatable :: vard4d

contains

  subroutine digda_read_restart(sdtx, lon_values, lat_values, depth_values, time_values)
    double precision, dimension(:), allocatable, intent(out) :: lon_values
    double precision, dimension(:), allocatable, intent(out) :: lat_values
    double precision, dimension(:), allocatable, intent(out) :: depth_values
    integer(kind=4), dimension(:), allocatable, intent(out) :: time_values
    double precision, dimension(nxydim, sdzdim, sedtdm), intent(out) :: sdtx

    integer :: ios

    write(nfstdo, *) "'--------------------------------------------------'"
    write(nfstdo, *) "'-- Started reading Sediment RESTART file        --'"
    write(nfstdo, *) "'--------------------------------------------------'"
    write(nfstdo, *) "  DIGDA TRACER No.01 = scla"
    write(nfstdo, *) "  DIGDA TRACER No.02 = spom"
    write(nfstdo, *) "  DIGDA TRACER No.03 = scar"
    write(nfstdo, *) "  DIGDA TRACER No.04 = sopa"
    write(nfstdo, *) "  DIGDA TRACER No.05 = sage"
    write(nfstdo, *) "  DIGDA TRACER No.06 = cnit"
    write(nfstdo, *) "  DIGDA TRACER No.07 = cdic"
    write(nfstdo, *) "  DIGDA TRACER No.08 = calk"
    write(nfstdo, *) "  DIGDA TRACER No.09 = csil"
    write(nfstdo, *) "  DIGDA TRACER No.10 = cdo2"
    write(nfstdo, *) "  DIGDA TRACER No.11 = s13pom"
    write(nfstdo, *) "  DIGDA TRACER No.12 = s13car"
    write(nfstdo, *) "  DIGDA TRACER No.13 = c13dic"
    write(nfstdo, *) "  DIGDA TRACER No.14 = s14pom"
    write(nfstdo, *) "  DIGDA TRACER No.15 = s14car"
    write(nfstdo, *) "  DIGDA TRACER No.16 = c14dic"
    write(nfstdo, *) "'--------------------------------------------------'"
    rewind(nfparm)
    read(nfparm, NML=nmocrt, IOSTAT=ios)
    write(nfstdo, nmocrt)

    !-- Call the read_netcdf subroutine --

    filename = trim(init_dir)//"/SEDTRR.nc"
    write(nfstdo, *) "'--------------------------------------------------'"
    write(nfstdo, *) "  Read initial condition file: ", filename
    write(nfstdo, *) "'--------------------------------------------------'"
    varname = "Clay"
    get_axis = .true.
    call read_netcdf_auto(filename, varname, vard4d, lon_values, lat_values, depth_values, time_values, get_axis)
    do k = 1, sdzdim
      do j = 1, nydim
        do i = 1, nxdim
          ij = i + (j - 1) * nxdim
          sdtx(ij, k, 1) = vard4d(i, j, k, 1)
        end do
      end do
    end do

    varname = "OrgM"
    get_axis = .false.
    call read_netcdf_auto(filename, varname, vard4d, lon_values, lat_values, depth_values, time_values, get_axis)
    do k = 1, sdzdim
      do j = 1, nydim
        do i = 1, nxdim
          ij = i + (j - 1) * nxdim
          sdtx(ij, k, 2) = vard4d(i, j, k, 1)
        end do
      end do
    end do

    varname = "CaCO3"
    get_axis = .false.
    call read_netcdf_auto(filename, varname, vard4d, lon_values, lat_values, depth_values, time_values, get_axis)
    do k = 1, sdzdim
      do j = 1, nydim
        do i = 1, nxdim
          ij = i + (j - 1) * nxdim
          sdtx(ij, k, 3) = vard4d(i, j, k, 1)
        end do
      end do
    end do

    varname = "SiO2"
    get_axis = .false.
    call read_netcdf_auto(filename, varname, vard4d, lon_values, lat_values, depth_values, time_values, get_axis)
    do k = 1, sdzdim
      do j = 1, nydim
        do i = 1, nxdim
          ij = i + (j - 1) * nxdim
          sdtx(ij, k, 4) = vard4d(i, j, k, 1)
        end do
      end do
    end do

    varname = "Agetracer"
    get_axis = .false.
    call read_netcdf_auto(filename, varname, vard4d, lon_values, lat_values, depth_values, time_values, get_axis)
    do k = 1, sdzdim
      do j = 1, nydim
        do i = 1, nxdim
          ij = i + (j - 1) * nxdim
          sdtx(ij, k, 5) = vard4d(i, j, k, 1)
        end do
      end do
    end do

    varname = "Phosphate"
    get_axis = .false.
    call read_netcdf_auto(filename, varname, vard4d, lon_values, lat_values, depth_values, time_values, get_axis)
    do k = 1, sdzdim
      do j = 1, nydim
        do i = 1, nxdim
          ij = i + (j - 1) * nxdim
          sdtx(ij, k, 6) = vard4d(i, j, k, 1)
        end do
      end do
    end do

    varname = "DIC"
    get_axis = .false.
    call read_netcdf_auto(filename, varname, vard4d, lon_values, lat_values, depth_values, time_values, get_axis)
    do k = 1, sdzdim
      do j = 1, nydim
        do i = 1, nxdim
          ij = i + (j - 1) * nxdim
          sdtx(ij, k, 7) = vard4d(i, j, k, 1)
        end do
      end do
    end do

    varname = "Alkalinity"
    get_axis = .false.
    call read_netcdf_auto(filename, varname, vard4d, lon_values, lat_values, depth_values, time_values, get_axis)
    do k = 1, sdzdim
      do j = 1, nydim
        do i = 1, nxdim
          ij = i + (j - 1) * nxdim
          sdtx(ij, k, 8) = vard4d(i, j, k, 1)
        end do
      end do
    end do

    varname = "Silicate"
    get_axis = .false.
    call read_netcdf_auto(filename, varname, vard4d, lon_values, lat_values, depth_values, time_values, get_axis)
    do k = 1, sdzdim
      do j = 1, nydim
        do i = 1, nxdim
          ij = i + (j - 1) * nxdim
          sdtx(ij, k, 9) = vard4d(i, j, k, 1)
        end do
      end do
    end do

    varname = "Oxygen"
    get_axis = .false.
    call read_netcdf_auto(filename, varname, vard4d, lon_values, lat_values, depth_values, time_values, get_axis)
    do k = 1, sdzdim
      do j = 1, nydim
        do i = 1, nxdim
          ij = i + (j - 1) * nxdim
          sdtx(ij, k, 10) = vard4d(i, j, k, 1)
        end do
      end do
    end do

    !-- Synchronize all processes --
    call mpi_barrier(mpi_comm_world, ierr)

    !-- Free memory for variables --
    if (allocated(vard3d)) deallocate(vard3d)
    if (allocated(vard4d)) deallocate(vard4d)

    write(nfstdo, *) "'--------------------------------------------------'"
    write(nfstdo, *) "'-- Finised reading Sediment RESTART file        --'"
    write(nfstdo, *) "'--------------------------------------------------'"

  end subroutine digda_read_restart


  subroutine digda_write_restart(sdtx, lon_values, lat_values, depth_values, time_values)

    double precision, dimension(nxydim, sdzdim, sedtdm), intent(in) :: sdtx
    double precision, dimension(:), allocatable, intent(in) :: lon_values
    double precision, dimension(:), allocatable, intent(in) :: lat_values
    double precision, dimension(:), allocatable, intent(in) :: depth_values
    integer(kind=4), dimension(:), allocatable, intent(in) :: time_values

    character(len=100)            :: filenames(sedtdm)

    integer :: ios

    !-- Allocate the array with appropriate dimensions --
    allocate(vard4d(nxdim, nydim, sdzdim, 1))

    write(nfstdo, *) "'--------------------------------------------------'"
    write(nfstdo, *) "'-- Started writing Sediment RESTART file        --'"
    write(nfstdo, *) "'--------------------------------------------------'"

    !-- Call the write_netcdf subroutine --

    filename = trim(rstt_dir)//"/SEDTRR.nc"
    write(nfstdo, *) "'--------------------------------------------------'"
    write(nfstdo, *) "  Write restart file: ", filename
    write(nfstdo, *) "'--------------------------------------------------'"

    if ( myrank == iroot ) then

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
    call write_netcdf_3dwt(filenames(1), varname, varunit, vard4d, lon_values, lat_values, depth_values, time_values)

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
    call write_netcdf_3dwt(filenames(2), varname, varunit, vard4d, lon_values, lat_values, depth_values, time_values)

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
    call write_netcdf_3dwt(filenames(3), varname, varunit, vard4d, lon_values, lat_values, depth_values, time_values)

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
    call write_netcdf_3dwt(filenames(4), varname, varunit, vard4d, lon_values, lat_values, depth_values, time_values)

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
    call write_netcdf_3dwt(filenames(5), varname, varunit, vard4d, lon_values, lat_values, depth_values, time_values)

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
    call write_netcdf_3dwt(filenames(6), varname, varunit, vard4d, lon_values, lat_values, depth_values, time_values)

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
    call write_netcdf_3dwt(filenames(7), varname, varunit, vard4d, lon_values, lat_values, depth_values, time_values)

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
    call write_netcdf_3dwt(filenames(8), varname, varunit, vard4d, lon_values, lat_values, depth_values, time_values)

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
    call write_netcdf_3dwt(filenames(9), varname, varunit, vard4d, lon_values, lat_values, depth_values, time_values)

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
    call write_netcdf_3dwt(filenames(10), varname, varunit, vard4d, lon_values, lat_values, depth_values, time_values)

   !call combine_multiple_netcdf_files(filenames, filename)

    end if

    !-- Synchronize all processes --
    call mpi_barrier(mpi_comm_world, ierr)

    !-- Free memory for variables --
    if (allocated(vard3d)) deallocate(vard3d)
    if (allocated(vard4d)) deallocate(vard4d)

    write(nfstdo, *) "'--------------------------------------------------'"
    write(nfstdo, *) "'-- Finised writing Sediment RESTART file        --'"
    write(nfstdo, *) "'--------------------------------------------------'"

  end subroutine digda_write_restart

  subroutine digda_read_bottom_forcing(bdep, bfrc)

    double precision, dimension(nxydim), intent(out) :: bdep
    double precision, dimension(nxydim, benfdm), intent(out) :: bfrc

    double precision, dimension(nxdim, nydim) :: omask

    integer :: ios

    write(nfstdo, *) "'--------------------------------------------------'"
    write(nfstdo, *) "'-- Started reading Sed-Water boundary forcing   --'"
    write(nfstdo, *) "'--------------------------------------------------'"
    rewind(nfparm)
    read(nfparm, NML=nmocbc, IOSTAT=ios)
    write(nfstdo, nmocbc)

    !-- Allocate the array with appropriate dimensions --
    allocate(vard2d(nxdim, nydim))
    allocate(vard3d(nxdim, nydim, 1))

    !-- Call the read_netcdf subroutine --
    !-- Ocean Depth
    filename = trim(bound_dir)//"/bdep.nc"
    varname = "DEPTH"
    get_axis = .false.
 !!!call read_netcdf_2dwt(filename, varname, vard3d, lon_values, lat_values, time_values)
   !call read_netcdf_2dwot(filename, varname, vard2d, lon_values, lat_values)
    call read_netcdf_auto(filename, varname, vard4d, lon_values, lat_values, depth_values, time_values, get_axis)
    do j = 1, nydim
      do i = 1, nxdim
        ij = i + (j - 1) * nxdim
       !bdep(ij) = vard3d(i, j, 1)
       !if (vard3d(i, j, 1) == rmiss) then
       !  omask(i, j) = 0.d0
       !else
       !  omask(i, j) = 1.d0
       !end if
       !bdep(ij) = vard2d(i, j)
        bdep(ij) = vard4d(i, j, 1, 1)
        if (vard2d(i, j) == rmiss) then
          omask(i, j) = 0.d0
        else
          omask(i, j) = 1.d0
        end if
      end do
    end do

    !-- Temperature
    filename = trim(bound_dir)//"/btem.nc"
    varname = "TEM"
    get_axis = .false.
   !call read_netcdf_2dwt(filename, varname, vard3d, lon_values, lat_values, time_values)
 !!!call read_netcdf_2dwot(filename, varname, vard2d, lon_values, lat_values)
    call read_netcdf_auto(filename, varname, vard4d, lon_values, lat_values, depth_values, time_values, get_axis)
    do j = 1, nydim
      do i = 1, nxdim
        ij = i + (j - 1) * nxdim
       !bfrc(ij, 6) = vard3d(i, j, 1)
     !!!bfrc(ij, 6) = vard2d(i, j)
        bfrc(ij, 6) = vard4d(i, j, 1, 1)
      end do
    end do

    !-- Salinity
    filename = trim(bound_dir)//"/bsal.nc"
    varname = "SAL"
    get_axis = .false.
   !call read_netcdf_2dwt(filename, varname, vard3d, lon_values, lat_values, time_values)
 !!!call read_netcdf_2dwot(filename, varname, vard2d, lon_values, lat_values)
    call read_netcdf_auto(filename, varname, vard4d, lon_values, lat_values, depth_values, time_values, get_axis)
    do j = 1, nydim
      do i = 1, nxdim
        ij = i + (j - 1) * nxdim
       !bfrc(ij, 7) = vard3d(i, j, 1)
     !!!bfrc(ij, 7) = vard2d(i, j)
        bfrc(ij, 7) = vard4d(i, j, 1, 1)
      end do
    end do

    !-- Phosphate
    filename = trim(bound_dir)//"/bnit.nc"
    varname = "NIT"
    get_axis = .false.
   !call read_netcdf_2dwt(filename, varname, vard3d, lon_values, lat_values, time_values)
 !!!call read_netcdf_2dwot(filename, varname, vard2d, lon_values, lat_values)
    call read_netcdf_auto(filename, varname, vard4d, lon_values, lat_values, depth_values, time_values, get_axis)
    do j = 1, nydim
      do i = 1, nxdim
        ij = i + (j - 1) * nxdim
       !bfrc(ij, 8) = vard3d(i, j, 1)
     !!!bfrc(ij, 8) = vard2d(i, j)
        bfrc(ij, 8) = vard4d(i, j, 1, 1)
      end do
    end do

    !-- DIC
    filename = trim(bound_dir)//"/bdic.nc"
    varname = "DIC"
    get_axis = .false.
   !call read_netcdf_2dwt(filename, varname, vard3d, lon_values, lat_values, time_values)
 !!!call read_netcdf_2dwot(filename, varname, vard2d, lon_values, lat_values)
    call read_netcdf_auto(filename, varname, vard4d, lon_values, lat_values, depth_values, time_values, get_axis)
    do j = 1, nydim
      do i = 1, nxdim
        ij = i + (j - 1) * nxdim
       !bfrc(ij, 9) = vard3d(i, j, 1)
     !!!bfrc(ij, 9) = vard2d(i, j)
        bfrc(ij, 9) = vard4d(i, j, 1, 1)
      end do
    end do

    !-- Alkalinity
    filename = trim(bound_dir)//"/balk.nc"
    varname = "ALK"
    get_axis = .false.
   !call read_netcdf_2dwt(filename, varname, vard3d, lon_values, lat_values, time_values)
 !!!call read_netcdf_2dwot(filename, varname, vard2d, lon_values, lat_values)
    call read_netcdf_auto(filename, varname, vard4d, lon_values, lat_values, depth_values, time_values, get_axis)
    do j = 1, nydim
      do i = 1, nxdim
        ij = i + (j - 1) * nxdim
       !bfrc(ij, 10) = vard3d(i, j, 1)
     !!!bfrc(ij, 10) = vard2d(i, j)
        bfrc(ij, 10) = vard4d(i, j, 1, 1)
      end do
    end do

    !-- Silicate
    filename = trim(bound_dir)//"/bsil.nc"
    varname = "SIL"
    get_axis = .false.
   !call read_netcdf_2dwt(filename, varname, vard3d, lon_values, lat_values, time_values)
 !!!call read_netcdf_2dwot(filename, varname, vard2d, lon_values, lat_values)
    call read_netcdf_auto(filename, varname, vard4d, lon_values, lat_values, depth_values, time_values, get_axis)
    do j = 1, nydim
      do i = 1, nxdim
        ij = i + (j - 1) * nxdim
       !bfrc(ij, 11) = vard3d(i, j, 1)
     !!!bfrc(ij, 11) = vard2d(i, j)
        bfrc(ij, 11) = vard4d(i, j, 1, 1)
      end do
    end do

    !-- Oxygen
    filename = trim(bound_dir)//"/bdo2.nc"
    varname = "DO2"
    get_axis = .false.
   !call read_netcdf_2dwt(filename, varname, vard3d, lon_values, lat_values, time_values)
 !!!call read_netcdf_2dwot(filename, varname, vard2d, lon_values, lat_values)
    call read_netcdf_auto(filename, varname, vard4d, lon_values, lat_values, depth_values, time_values, get_axis)
    do j = 1, nydim
      do i = 1, nxdim
        ij = i + (j - 1) * nxdim
       !bfrc(ij, 12) = vard3d(i, j, 1)
     !!!bfrc(ij, 12) = vard2d(i, j)
        bfrc(ij, 12) = vard4d(i, j, 1, 1)
      end do
    end do

    !-- 13C-DIC
    filename = trim(bound_dir)//"/bc13.nc"
    varname = "C13"
    get_axis = .false.
   !call read_netcdf_2dwt(filename, varname, vard3d, lon_values, lat_values, time_values)
 !!!call read_netcdf_2dwot(filename, varname, vard2d, lon_values, lat_values)
    call read_netcdf_auto(filename, varname, vard4d, lon_values, lat_values, depth_values, time_values, get_axis)
    do j = 1, nydim
      do i = 1, nxdim
        ij = i + (j - 1) * nxdim
       !bfrc(ij, 13) = vard3d(i, j, 1)
     !!!bfrc(ij, 13) = vard2d(i, j)
        bfrc(ij, 13) = vard4d(i, j, 1, 1)
      end do
    end do

    !-- 14C-DIC
    filename = trim(bound_dir)//"/bc14.nc"
    varname = "C14"
    get_axis = .false.
   !call read_netcdf_2dwt(filename, varname, vard3d, lon_values, lat_values, time_values)
 !!!call read_netcdf_2dwot(filename, varname, vard2d, lon_values, lat_values)
    call read_netcdf_auto(filename, varname, vard4d, lon_values, lat_values, depth_values, time_values, get_axis)
    do j = 1, nydim
      do i = 1, nxdim
        ij = i + (j - 1) * nxdim
       !bfrc(ij, 14) = vard3d(i, j, 1)
     !!!bfrc(ij, 14) = vard2d(i, j)
        bfrc(ij, 14) = vard4d(i, j, 1, 1)
      end do
    end do

    !-- Clay
    filename = trim(bound_dir)//"/ficla.nc"
    varname = "FICLA"
    get_axis = .false.
   !call read_netcdf_2dwt(filename, varname, vard3d, lon_values, lat_values, time_values)
 !!!call read_netcdf_2dwot(filename, varname, vard2d, lon_values, lat_values)
    call read_netcdf_auto(filename, varname, vard4d, lon_values, lat_values, depth_values, time_values, get_axis)
    do j = 1, nydim
      do i = 1, nxdim
        ij = i + (j - 1) * nxdim
       !bfrc(ij, 1) = vard3d(i, j, 1) * omask(i, j)
     !!!bfrc(ij, 1) = vard2d(i, j) * omask(i, j)
        bfrc(ij, 1) = vard4d(i, j, 1, 1) * omask(i, j)
      end do
    end do

    !-- POM
    filename = trim(bound_dir)//"/fipom.nc"
    varname = "FIPOM"
    get_axis = .false.
   !call read_netcdf_2dwt(filename, varname, vard3d, lon_values, lat_values, time_values)
 !!!call read_netcdf_2dwot(filename, varname, vard2d, lon_values, lat_values)
    call read_netcdf_auto(filename, varname, vard4d, lon_values, lat_values, depth_values, time_values, get_axis)
    do j = 1, nydim
      do i = 1, nxdim
        ij = i + (j - 1) * nxdim
       !bfrc(ij, 2) = vard3d(i, j, 1) * omask(i, j)
     !!!bfrc(ij, 2) = vard2d(i, j) * omask(i, j)
        bfrc(ij, 2) = vard4d(i, j, 1, 1) * omask(i, j)
      end do
    end do

    !-- CaCO3
    filename = trim(bound_dir)//"/ficar.nc"
    varname = "FICAR"
    get_axis = .false.
   !call read_netcdf_2dwt(filename, varname, vard3d, lon_values, lat_values, time_values)
 !!!call read_netcdf_2dwot(filename, varname, vard2d, lon_values, lat_values)
    call read_netcdf_auto(filename, varname, vard4d, lon_values, lat_values, depth_values, time_values, get_axis)
    do j = 1, nydim
      do i = 1, nxdim
        ij = i + (j - 1) * nxdim
       !bfrc(ij, 3) = vard3d(i, j, 1) * omask(i, j)
     !!!bfrc(ij, 3) = vard2d(i, j) * omask(i, j)
        bfrc(ij, 3) = vard4d(i, j, 1, 1) * omask(i, j)
      end do
    end do

    !-- SiO2
    filename = trim(bound_dir)//"/fiopa.nc"
    varname = "FIOPA"
    get_axis = .false.
   !call read_netcdf_2dwt(filename, varname, vard3d, lon_values, lat_values, time_values)
 !1!call read_netcdf_2dwot(filename, varname, vard2d, lon_values, lat_values)
    call read_netcdf_auto(filename, varname, vard4d, lon_values, lat_values, depth_values, time_values, get_axis)
    do j = 1, nydim
      do i = 1, nxdim
        ij = i + (j - 1) * nxdim
       !bfrc(ij, 4) = vard3d(i, j, 1) * omask(i, j)
     !!!bfrc(ij, 4) = vard2d(i, j) * omask(i, j)
        bfrc(ij, 4) = vard4d(i, j, 1, 1) * omask(i, j)
      end do
    end do

    !-- Age tracer
    do j = 1, nydim
      do i = 1, nxdim
        ij = i + (j - 1) * nxdim
        vard4d(i, j, 1, 1) = 0.d0
        bfrc(ij, 5) = vard4d(i, j, 1, 1) * omask(i, j)
      end do
    end do

    !-- Synchronize all processes --
    call mpi_barrier(mpi_comm_world, ierr)

    !-- Free memory for variables --
    if (allocated(vard2d)) deallocate(vard2d)
    if (allocated(vard3d)) deallocate(vard3d)
    if (allocated(vard4d)) deallocate(vard4d)

    write(nfstdo, *) "'--------------------------------------------------'"
    write(nfstdo, *) "'-- Finished reading Sed-Water boundary forcing  --'"
    write(nfstdo, *) "'--------------------------------------------------'"

  end subroutine digda_read_bottom_forcing

  subroutine digda_write_data()

    if ( myrank == iroot ) then

    write(nfstdo, *) "'--------------------------------------------------'"
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
    write(nfstdo, *) "'-- Finished writing NetCDF file                 --'"
    write(nfstdo, *) "'--------------------------------------------------'"

    end if

  end subroutine digda_write_data

! subroutine read_namelist()
!   integer :: ios
!   character(len=100) :: namelist_file

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
