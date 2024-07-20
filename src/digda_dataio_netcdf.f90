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

  double precision, dimension(:,:), allocatable :: vard2d
  double precision, dimension(:,:,:,:), allocatable :: vard3d

contains

  subroutine digda_read_data(bdep, bfrc)

    double precision, dimension(nxydim), intent(out) :: bdep
    double precision, dimension(nxydim, benfdm), intent(out) :: bfrc

    double precision, dimension(nxdim, nydim) :: omask

    integer :: ios

    write(nfstdo, *) "'--------------------------------------------------'"
    write(nfstdo, *) "'-- Started reading Sed-Water boundary forcing   --'"
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
    read(nfparm, NML=nmocbc, IOSTAT=ios)
    write(nfstdo, nmocbc)

    !-- Allocate the array with appropriate dimensions --
    allocate(vard2d(nxdim, nydim))

    !-- Call the read_netcdf subroutine --
    filename = trim(bound_dir)//"/bdep.nc"
    varname = "depth"
    call read_netcdf_2d(filename, varname, vard2d, lon_values, lat_values)

    do j = 1, nydim
      do i = 1, nxdim
        ij = i + (j - 1) * nxdim
        bdep(ij) = vard2d(i, j)
        if (vard2d(i, j) == rmiss) then
          omask(i, j) = 0.d0
        else
          omask(i, j) = 1.d0
        end if
      end do
    end do

    filename = trim(bound_dir)//"/btem.nc"
    varname = "TEM"
    call read_netcdf_2d(filename, varname, vard2d, lon_values, lat_values)
    do j = 1, nydim
      do i = 1, nxdim
        ij = i + (j - 1) * nxdim
        bfrc(ij, 6) = vard2d(i, j)
      end do
    end do

    filename = trim(bound_dir)//"/bsal.nc"
    varname = "SAL"
    call read_netcdf_2d(filename, varname, vard2d, lon_values, lat_values)
    do j = 1, nydim
      do i = 1, nxdim
        ij = i + (j - 1) * nxdim
        bfrc(ij, 7) = vard2d(i, j)
      end do
    end do

    filename = trim(bound_dir)//"/bnit.nc"
    varname = "NIT"
    call read_netcdf_2d(filename, varname, vard2d, lon_values, lat_values)
    do j = 1, nydim
      do i = 1, nxdim
        ij = i + (j - 1) * nxdim
        bfrc(ij, 8) = vard2d(i, j)
      end do
    end do

    filename = trim(bound_dir)//"/bdic.nc"
    varname = "DIC"
    call read_netcdf_2d(filename, varname, vard2d, lon_values, lat_values)
    do j = 1, nydim
      do i = 1, nxdim
        ij = i + (j - 1) * nxdim
        bfrc(ij, 9) = vard2d(i, j)
      end do
    end do

    filename = trim(bound_dir)//"/balk.nc"
    varname = "ALK"
    call read_netcdf_2d(filename, varname, vard2d, lon_values, lat_values)
    do j = 1, nydim
      do i = 1, nxdim
        ij = i + (j - 1) * nxdim
        bfrc(ij, 10) = vard2d(i, j)
      end do
    end do

    filename = trim(bound_dir)//"/bsil.nc"
    varname = "SIL"
    call read_netcdf_2d(filename, varname, vard2d, lon_values, lat_values)
    do j = 1, nydim
      do i = 1, nxdim
        ij = i + (j - 1) * nxdim
        bfrc(ij, 11) = vard2d(i, j)
      end do
    end do

    filename = trim(bound_dir)//"/bdo2.nc"
    varname = "DO2"
    call read_netcdf_2d(filename, varname, vard2d, lon_values, lat_values)
    do j = 1, nydim
      do i = 1, nxdim
        ij = i + (j - 1) * nxdim
        bfrc(ij, 12) = vard2d(i, j)
      end do
    end do

    filename = trim(bound_dir)//"/bc13.nc"
    varname = "C13"
    call read_netcdf_2d(filename, varname, vard2d, lon_values, lat_values)
    do j = 1, nydim
      do i = 1, nxdim
        ij = i + (j - 1) * nxdim
        bfrc(ij, 13) = vard2d(i, j)
      end do
    end do

    filename = trim(bound_dir)//"/bc14.nc"
    varname = "C14"
    call read_netcdf_2d(filename, varname, vard2d, lon_values, lat_values)
    do j = 1, nydim
      do i = 1, nxdim
        ij = i + (j - 1) * nxdim
        bfrc(ij, 14) = vard2d(i, j)
      end do
    end do

    filename = trim(bound_dir)//"/ficla.nc"
    varname = "FICLA"
    call read_netcdf_2d(filename, varname, vard2d, lon_values, lat_values)
    do j = 1, nydim
      do i = 1, nxdim
        ij = i + (j - 1) * nxdim
        bfrc(ij, 1) = vard2d(i, j) * omask(i, j)
      end do
    end do

    filename = trim(bound_dir)//"/fipom.nc"
    varname = "FIPOM"
    call read_netcdf_2d(filename, varname, vard2d, lon_values, lat_values)
    do j = 1, nydim
      do i = 1, nxdim
        ij = i + (j - 1) * nxdim
        bfrc(ij, 2) = vard2d(i, j) * omask(i, j)
      end do
    end do

    filename = trim(bound_dir)//"/ficar.nc"
    varname = "FICAR"
    call read_netcdf_2d(filename, varname, vard2d, lon_values, lat_values)
    do j = 1, nydim
      do i = 1, nxdim
        ij = i + (j - 1) * nxdim
        bfrc(ij, 3) = vard2d(i, j) * omask(i, j)
      end do
    end do

    filename = trim(bound_dir)//"/fiopa.nc"
    varname = "FIOPA"
    call read_netcdf_2d(filename, varname, vard2d, lon_values, lat_values)
    do j = 1, nydim
      do i = 1, nxdim
        ij = i + (j - 1) * nxdim
        bfrc(ij, 4) = vard2d(i, j) * omask(i, j)
      end do
    end do

    do j = 1, nydim
      do i = 1, nxdim
        ij = i + (j - 1) * nxdim
        vard2d(i, j) = 0.d0
        bfrc(ij, 5) = vard2d(i, j) * omask(i, j)
      end do
    end do

    !-- Synchronize all processes --
    call mpi_barrier(mpi_comm_world, ierr)

    !-- Free memory for variables --
    if (allocated(vard2d)) then
      deallocate(vard2d)
    end if

    write(nfstdo, *) "'--------------------------------------------------'"
    write(nfstdo, *) "'-- Finished reading Sed-Water boundary forcing  --'"
    write(nfstdo, *) "'--------------------------------------------------'"

  end subroutine digda_read_data

  subroutine digda_write_data()

    if ( myrank == iroot ) then

    write(nfstdo, *) "'--------------------------------------------------'"
    write(nfstdo, *) "'-- Started writing NetCDF file                  --'"
    write(nfstdo, *) "'--------------------------------------------------'"

    !-- Allocate the array with appropriate dimensions --
    allocate(vard3d(nxdim, nydim, sdzdim, 1))

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
  !   call write_netcdf_3d(filename, varname, varunit, vard3d, lon_values, lat_values, lz, 999)

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
  !   call write_netcdf_3d(filename, varname, varunit, vard3d, lon_values, lat_values, lz, 999)

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
  !   call write_netcdf_3d(filename, varname, varunit, vard3d, lon_values, lat_values, lz, 999)

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
  !   call write_netcdf_3d(filename, varname, varunit, vard3d, lon_values, lat_values, lz, 999)

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
  !   call write_netcdf_3d(filename, varname, varunit, vard3d, lon_values, lat_values, lz, 999)

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
  !   call write_netcdf_3d(filename, varname, varunit, vard3d, lon_values, lat_values, lz, 999)

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
  !   call write_netcdf_3d(filename, varname, varunit, vard3d, lon_values, lat_values, lz, 999)

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
  !   call write_netcdf_3d(filename, varname, varunit, vard3d, lon_values, lat_values, lz, 999)

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
  !   call write_netcdf_3d(filename, varname, varunit, vard3d, lon_values, lat_values, lz, 999)

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
  !   call write_netcdf_3d(filename, varname, varunit, vard3d, lon_values, lat_values, lz, 999)


    ! Free memory for variables
    if (allocated(vard3d)) then
      deallocate(vard3d)
    end if

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
