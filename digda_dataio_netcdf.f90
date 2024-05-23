program digda_datadigda_dataio_gtool
  use netcdf_utilities
  use module_sdmain
  implicit none

  character(len=100) :: filename, varname, varunit

  integer :: i, j, k, l

  double precision, dimension(:,:), allocatable :: bdep    !! ocean depth [m]
! double precision, dimension(nxdim, nydim) :: bdep    !! ocean depth [m]
  double precision, dimension(nxdim, nydim) :: omask   !! ocean depth [m]

  double precision, dimension(:,:), allocatable :: ficla   !! benthic clay flux [g/m3*cm/s]
  double precision, dimension(:,:), allocatable :: fipom   !! benthic pom flux [mMol/m3*cm/s]
  double precision, dimension(:,:), allocatable :: ficar   !! benthic caco3 flux [mMol/m3*cm/s]
  double precision, dimension(:,:), allocatable :: fiopa   !! benthic opal flux [mMol/m3*cm/s]
  double precision, dimension(:,:), allocatable :: fiage   !!
! double precision, dimension(nxdim, nydim) :: ficla   !! benthic clay flux [g/m3*cm/s]
! double precision, dimension(nxdim, nydim) :: fipom   !! benthic pom flux [mMol/m3*cm/s]
! double precision, dimension(nxdim, nydim) :: ficar   !! benthic caco3 flux [mMol/m3*cm/s]
! double precision, dimension(nxdim, nydim) :: fiopa   !! benthic opal flux [mMol/m3*cm/s]
! double precision, dimension(nxdim, nydim) :: fiage   !!

! double precision, dimension(:,:), allocatable :: outfocla  !! benthic clay flux [g/m3*cm/s]
! double precision, dimension(:,:), allocatable :: outfopom  !! benthic pom flux [mMol/m3*cm/s]
! double precision, dimension(:,:), allocatable :: outfocar  !! benthic caco3 flux [mMol/m3*cm/s]
! double precision, dimension(:,:), allocatable :: outfoopa  !! benthic opal flux [mMol/m3*cm/s]
! double precision, dimension(:,:), allocatable :: outfoage  !!
  double precision, dimension(nxdim, nydim) :: outfocla  !! benthic clay flux [g/m3*cm/s]
  double precision, dimension(nxdim, nydim) :: outfopom  !! benthic pom flux [mMol/m3*cm/s]
  double precision, dimension(nxdim, nydim) :: outfocar  !! benthic caco3 flux [mMol/m3*cm/s]
  double precision, dimension(nxdim, nydim) :: outfoopa  !! benthic opal flux [mMol/m3*cm/s]
  double precision, dimension(nxdim, nydim) :: outfoage  !!

  double precision, dimension(:,:), allocatable :: bwtem   !! temperature [mMol/m3]
  double precision, dimension(:,:), allocatable :: bwsal   !! salinity [mMol/m3]
  double precision, dimension(:,:), allocatable :: bwnit   !! phosphate [mMol/m3]
  double precision, dimension(:,:), allocatable :: bwdic   !! DIC [mMol/m3]
  double precision, dimension(:,:), allocatable :: bwalk   !! alkalinity [mMol/m3]
  double precision, dimension(:,:), allocatable :: bwdo2   !! oxygen [mMol/m3]
  double precision, dimension(:,:), allocatable :: bwsil   !! silicate [mMol/m3]o
  double precision, dimension(:,:), allocatable :: bw13dic !! 13C-DIC [mMol/m3]
  double precision, dimension(:,:), allocatable :: bw14dic !! 14C-DIC [mMol/m3]
! double precision, dimension(nxdim, nydim) :: bwtem   !! temperature [mMol/m3]
! double precision, dimension(nxdim, nydim) :: bwsal   !! salinity [mMol/m3]
! double precision, dimension(nxdim, nydim) :: bwnit   !! phosphate [mMol/m3]
! double precision, dimension(nxdim, nydim) :: bwdic   !! DIC [mMol/m3]
! double precision, dimension(nxdim, nydim) :: bwalk   !! alkalinity [mMol/m3]
! double precision, dimension(nxdim, nydim) :: bwdo2   !! oxygen [mMol/m3]
! double precision, dimension(nxdim, nydim) :: bwsil   !! silicate [mMol/m3]o
! double precision, dimension(nxdim, nydim) :: bw13dic !! 13C-DIC [mMol/m3]
! double precision, dimension(nxdim, nydim) :: bw14dic !! 14C-DIC [mMol/m3]

! double precision, dimension(:,:), allocatable :: bfrc !!
  double precision, dimension(nxydim, benfdm) :: bfrc !!
! double precision, dimension(:,:,:), allocatable :: sdtx !!
  double precision, dimension(nxydim, sdzdim, sedtdm) :: sdtx !!

  double precision, dimension(:), allocatable :: lon_values !!
  double precision, dimension(:), allocatable :: lat_values !!
  double precision, dimension(:), allocatable :: depth_values !!

  double precision, dimension(:,:,:,:), allocatable :: outscla !!
  double precision, dimension(:,:,:,:), allocatable :: outspom !!
  double precision, dimension(:,:,:,:), allocatable :: outscar !!
  double precision, dimension(:,:,:,:), allocatable :: outsopa !!
  double precision, dimension(:,:,:,:), allocatable :: outsage !!
  double precision, dimension(:,:,:,:), allocatable :: outcnit !!
  double precision, dimension(:,:,:,:), allocatable :: outcdic !!
  double precision, dimension(:,:,:,:), allocatable :: outcalk !!
  double precision, dimension(:,:,:,:), allocatable :: outcsil !!
  double precision, dimension(:,:,:,:), allocatable :: outcdo2 !!

  double precision, dimension(:,:,:,:), allocatable :: out13spom !!
  double precision, dimension(:,:,:,:), allocatable :: out13scar !!
  double precision, dimension(:,:,:,:), allocatable :: out13cdic !!
  double precision, dimension(:,:,:,:), allocatable :: out14spom !!
  double precision, dimension(:,:,:,:), allocatable :: out14scar !!
  double precision, dimension(:,:,:,:), allocatable :: out14cdic !!

  ! Allocate the array with appropriate dimensions
  allocate(bdep(nxdim, nydim))
  allocate(ficla(nxdim, nydim))
  allocate(fipom(nxdim, nydim))
  allocate(ficar(nxdim, nydim))
  allocate(fiopa(nxdim, nydim))
  allocate(fiage(nxdim, nydim))
  allocate(bwtem(nxdim, nydim))
  allocate(bwsal(nxdim, nydim))
  allocate(bwnit(nxdim, nydim))
  allocate(bwdic(nxdim, nydim))
  allocate(bwalk(nxdim, nydim))
  allocate(bwdo2(nxdim, nydim))
  allocate(bwsil(nxdim, nydim))
  allocate(bw13dic(nxdim, nydim))
  allocate(bw14dic(nxdim, nydim))

  allocate(outscla(nxdim, nydim, sdzdim, 1))
  allocate(outspom(nxdim, nydim, sdzdim, 1))
  allocate(outscar(nxdim, nydim, sdzdim, 1))
  allocate(outsopa(nxdim, nydim, sdzdim, 1))
  allocate(outsage(nxdim, nydim, sdzdim, 1))
  allocate(outcnit(nxdim, nydim, sdzdim, 1))
  allocate(outcdic(nxdim, nydim, sdzdim, 1))
  allocate(outcalk(nxdim, nydim, sdzdim, 1))
  allocate(outcsil(nxdim, nydim, sdzdim, 1))
  allocate(outcdo2(nxdim, nydim, sdzdim, 1))

  ! Call the read_netcdf subroutine
  filename = trim(bound_dir)//"bdep.nc"
  varname = "depth"
  call read_netcdf_2d(filename, varname, bdep, lon_values, lat_values)
  filename = trim(bound_dir)//"btem.nc"
  varname = "TEM"
  call read_netcdf_2d(filename, varname, bwtem, lon_values, lat_values)
  filename = trim(bound_dir)//"bsal.nc"
  varname = "SAL"
  call read_netcdf_2d(filename, varname, bwsal, lon_values, lat_values)
  filename = trim(bound_dir)//"bnit.nc"
  varname = "NIT"
  call read_netcdf_2d(filename, varname, bwnit, lon_values, lat_values)
  filename = trim(bound_dir)//"bdic.nc"
  varname = "DIC"
  call read_netcdf_2d(filename, varname, bwdic, lon_values, lat_values)
  filename = trim(bound_dir)//"balk.nc"
  varname = "ALK"
  call read_netcdf_2d(filename, varname, bwalk, lon_values, lat_values)
  filename = trim(bound_dir)//"bdo2.nc"
  varname = "DO2"
  call read_netcdf_2d(filename, varname, bwdo2, lon_values, lat_values)
  filename = trim(bound_dir)//"bsil.nc"
  varname = "SIL"
  call read_netcdf_2d(filename, varname, bwsil, lon_values, lat_values)
  filename = trim(bound_dir)//"bc13.nc"
  varname = "C13"
  call read_netcdf_2d(filename, varname, bw13dic, lon_values, lat_values)
  filename = trim(bound_dir)//"bc14.nc"
  varname = "C14"
  call read_netcdf_2d(filename, varname, bw14dic, lon_values, lat_values)
  filename = trim(bound_dir)//"ficla.nc"
  varname = "FICLA"
  call read_netcdf_2d(filename, varname, ficla, lon_values, lat_values)
  filename = trim(bound_dir)//"fipom.nc"
  varname = "FIPOM"
  call read_netcdf_2d(filename, varname, fipom, lon_values, lat_values)
  filename = trim(bound_dir)//"ficar.nc"
  varname = "FICAR"
  call read_netcdf_2d(filename, varname, ficar, lon_values, lat_values)
  filename = trim(bound_dir)//"fiopa.nc"
  varname = "FIOPA"
  call read_netcdf_2d(filename, varname, fiopa, lon_values, lat_values)

  fiage = 0.d0

  do j = 1, nydim
    do i = 1, nxdim
      if (bdep(i, j) == rmiss) then
        omask(i, j) = 0.D0
      else
        omask(i, j) = 1.D0
      end if
    end do
  end do

  do j = 1, nydim
    do i = 1, nxdim
      l = i + (j - 1) * nxdim
      if (ficla(i, j) == rmiss) then
        ficla(i, j) = 0.D0
      end if
      if (fipom(i, j) == rmiss) then
        fipom(i, j) = 0.D0
      end if
      if (ficar(i, j) == rmiss) then
        ficar(i, j) = 0.D0
      end if
      if (fiopa(i, j) == rmiss) then
        fiopa(i, j) = 0.D0
      end if
      if (fiage(i, j) == rmiss) then
        fiage(i, j) = 0.D0
      end if
      bfrc(l, 1) = ficla(i, j) * omask(i, j)
      bfrc(l, 2) = fipom(i, j) * omask(i, j)
      bfrc(l, 3) = ficar(i, j) * omask(i, j)
      bfrc(l, 4) = fiopa(i, j) * omask(i, j)
      bfrc(l, 5) = fiage(i, j) * omask(i, j)
      bfrc(l, 6) = bwtem(i, j)
      bfrc(l, 7) = bwsal(i, j)
      bfrc(l, 8) = bwnit(i, j)
      bfrc(l, 9) = bwdic(i, j)
      bfrc(l, 10) = bwalk(i, j)
      bfrc(l, 11) = bwsil(i, j)
      bfrc(l, 12) = bwdo2(i, j)
      bfrc(l, 13) = bw13dic(i, j)
      bfrc(l, 14) = bw14dic(i, j)
    end do
  end do

  call sd( &
 &     bdep, &
 &     bfrc, &
 &     sdtx )

  do k = 1, sdzdim
    do j = 1, nydim
      do i = 1, nxdim
        l = i + (j - 1) * nxdim
        outscla(i, j, k, 1) = sdtx(l, k, 1)
        outspom(i, j, k, 1) = sdtx(l, k, 2)
        outscar(i, j, k, 1) = sdtx(l, k, 3)
        outsopa(i, j, k, 1) = sdtx(l, k, 4)
        outsage(i, j, k, 1) = sdtx(l, k, 5)
        outcnit(i, j, k, 1) = sdtx(l, k, 6)
        outcdic(i, j, k, 1) = sdtx(l, k, 7)
        outcalk(i, j, k, 1) = sdtx(l, k, 8)
        outcsil(i, j, k, 1) = sdtx(l, k, 9)
        outcdo2(i, j, k, 1) = sdtx(l, k, 10)
!       out13spom(i, j, k, 1) = sdtx(l, k, 11)
!       out13scar(i, j, k, 1) = sdtx(l, k, 12)
!       out13cdic(i, j, k, 1) = sdtx(l, k, 13)
!       out14spom(i, j, k, 1) = sdtx(l, k, 14)
!       out14scar(i, j, k, 1) = sdtx(l, k, 15)
!       out14cdic(i, j, k, 1) = sdtx(l, k, 16)
      end do
    end do
  end do

  if ( myrank == iroot ) then

    ! Call the write_netcdf subroutine
    filename = "scla.nc"
    varname = "Clay"
    varunit = "g m-3"
    call write_netcdf_3d(filename, varname, varunit, outscla, lon_values, lat_values, lz, 999)

    filename = "spom.nc"
    varname = "OrgM"
    varunit = "mmol m-3"
    call write_netcdf_3d(filename, varname, varunit, outspom, lon_values, lat_values, lz, 999)

    filename = "scar.nc"
    varname = "CaCO3"
    varunit = "mmol m-3"
    call write_netcdf_3d(filename, varname, varunit, outscar, lon_values, lat_values, lz, 999)

    filename = "sopa.nc"
    varname = "SiO2"
    varunit = "mmol m-3"
    call write_netcdf_3d(filename, varname, varunit, outsopa, lon_values, lat_values, lz, 999)

    filename = "sage.nc"
    varname = "Age"
    varunit = "yr"
    call write_netcdf_3d(filename, varname, varunit, outsage, lon_values, lat_values, lz, 999)

    filename = "cnit.nc"
    varname = "Phosphate"
    varunit = "mmol m-3"
    call write_netcdf_3d(filename, varname, varunit, outcnit, lon_values, lat_values, lz, 999)

    filename = "cdic.nc"
    varname = "DIC"
    varunit = "mmol m-3"
    call write_netcdf_3d(filename, varname, varunit, outcdic, lon_values, lat_values, lz, 999)

    filename = "calk.nc"
    varname = "Alkalinity"
    varunit = "mmol m-3"
    call write_netcdf_3d(filename, varname, varunit, outcalk, lon_values, lat_values, lz, 999)

    filename = "cdo2.nc"
    varname = "Oxygen"
    varunit = "mmol m-3"
    call write_netcdf_3d(filename, varname, varunit, outcdo2, lon_values, lat_values, lz, 999)

    filename = "csil.nc"
    varname = "Silicate"
    varunit = "mmol m-3"
    call write_netcdf_3d(filename, varname, varunit, outcsil, lon_values, lat_values, lz, 999)

  end if

  ! Synchronize all processes
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  ! Free memory for variables
  if (allocated(bdep)) then
    deallocate(bdep)
  end if
  if (allocated(ficla)) then
    deallocate(ficla)
  end if
  if (allocated(fipom)) then
    deallocate(fipom)
  end if
  if (allocated(ficar)) then
    deallocate(ficar)
  end if
  if (allocated(fiopa)) then
    deallocate(fiopa)
  end if
  if (allocated(fiage)) then
    deallocate(fiage)
  end if
  if (allocated(bwtem)) then
    deallocate(bwtem)
  end if
  if (allocated(bwsal)) then
    deallocate(bwsal)
  end if
  if (allocated(bwnit)) then
    deallocate(bwnit)
  end if
  if (allocated(bwdic)) then
    deallocate(bwdic)
  end if
  if (allocated(bwalk)) then
    deallocate(bwalk)
  end if
  if (allocated(bwdo2)) then
    deallocate(bwdo2)
  end if
  if (allocated(bwsil)) then
    deallocate(bwsil)
  end if
  if (allocated(bw13dic)) then
    deallocate(bw13dic)
  end if
  if (allocated(bw14dic)) then
    deallocate(bw14dic)
  end if

  ! Free memory for variables
  if (allocated(outscla)) then
    deallocate(outscla)
  end if
  if (allocated(outspom)) then
    deallocate(outspom)
  end if
  if (allocated(outscar)) then
    deallocate(outscar)
  end if
  if (allocated(outsopa)) then
    deallocate(outsopa)
  end if
  if (allocated(outsage)) then
    deallocate(outsage)
  end if
  if (allocated(outcnit)) then
    deallocate(outcnit)
  end if
  if (allocated(outcdic)) then
    deallocate(outcdic)
  end if
  if (allocated(outcalk)) then
    deallocate(outcalk)
  end if
  if (allocated(outcsil)) then
    deallocate(outcsil)
  end if
  if (allocated(outcdo2)) then
    deallocate(outcdo2)
  end if

  ! Finalize MPI
  call MPI_FINALIZE(ierr)


end program digda_datadigda_dataio_gtool
