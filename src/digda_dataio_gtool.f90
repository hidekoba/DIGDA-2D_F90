program digda_dataio_gtool
  use gtool_utilities
  use module_digda_cprdc
  implicit none

  character(len=100) :: filename, varname, vartitle, varunit

  integer :: i, j, k, l

  double precision, dimension(:), allocatable :: bdep    !! ocean depth [m]
! double precision, dimension(nxydim) :: bdep    !! ocean depth [m]
  double precision, dimension(nxydim) :: omask   !! ocean depth [m]

  double precision, dimension(:), allocatable :: ficla   !! benthic clay flux [g/m3*cm/s]
  double precision, dimension(:), allocatable :: fipom   !! benthic pom flux [mMol/m3*cm/s]
  double precision, dimension(:), allocatable :: ficar   !! benthic caco3 flux [mMol/m3*cm/s]
  double precision, dimension(:), allocatable :: fiopa   !! benthic opal flux [mMol/m3*cm/s]
  double precision, dimension(:), allocatable :: fiage   !!
! double precision, dimension(nxydim) :: ficla   !! benthic clay flux [g/m3*cm/s]
! double precision, dimension(nxydim) :: fipom   !! benthic pom flux [mMol/m3*cm/s]
! double precision, dimension(nxydim) :: ficar   !! benthic caco3 flux [mMol/m3*cm/s]
! double precision, dimension(nxydim) :: fiopa   !! benthic opal flux [mMol/m3*cm/s]
! double precision, dimension(nxydim) :: fiage   !!

! double precision, dimension(:), allocatable :: outfocla  !! benthic clay flux [g/m3*cm/s]
! double precision, dimension(:), allocatable :: outfopom  !! benthic pom flux [mMol/m3*cm/s]
! double precision, dimension(:), allocatable :: outfocar  !! benthic caco3 flux [mMol/m3*cm/s]
! double precision, dimension(:), allocatable :: outfoopa  !! benthic opal flux [mMol/m3*cm/s]
! double precision, dimension(:), allocatable :: outfoage  !!
  double precision, dimension(nxydim) :: outfocla  !! benthic clay flux [g/m3*cm/s]
  double precision, dimension(nxydim) :: outfopom  !! benthic pom flux [mMol/m3*cm/s]
  double precision, dimension(nxydim) :: outfocar  !! benthic caco3 flux [mMol/m3*cm/s]
  double precision, dimension(nxydim) :: outfoopa  !! benthic opal flux [mMol/m3*cm/s]
  double precision, dimension(nxydim) :: outfoage  !!

  double precision, dimension(:), allocatable :: bwtem   !! temperature [mMol/m3]
  double precision, dimension(:), allocatable :: bwsal   !! salinity [mMol/m3]
  double precision, dimension(:), allocatable :: bwnit   !! phosphate [mMol/m3]
  double precision, dimension(:), allocatable :: bwdic   !! DIC [mMol/m3]
  double precision, dimension(:), allocatable :: bwalk   !! alkalinity [mMol/m3]
  double precision, dimension(:), allocatable :: bwdo2   !! oxygen [mMol/m3]
  double precision, dimension(:), allocatable :: bwsil   !! silicate [mMol/m3]o
  double precision, dimension(:), allocatable :: bw13dic !! 13C-DIC [mMol/m3]
  double precision, dimension(:), allocatable :: bw14dic !! 14C-DIC [mMol/m3]
! double precision, dimension(nxydim) :: bwtem   !! temperature [mMol/m3]
! double precision, dimension(nxydim) :: bwsal   !! salinity [mMol/m3]
! double precision, dimension(nxydim) :: bwnit   !! phosphate [mMol/m3]
! double precision, dimension(nxydim) :: bwdic   !! DIC [mMol/m3]
! double precision, dimension(nxydim) :: bwalk   !! alkalinity [mMol/m3]
! double precision, dimension(nxydim) :: bwdo2   !! oxygen [mMol/m3]
! double precision, dimension(nxydim) :: bwsil   !! silicate [mMol/m3]o
! double precision, dimension(nxydim) :: bw13dic !! 13C-DIC [mMol/m3]
! double precision, dimension(nxydim) :: bw14dic !! 14C-DIC [mMol/m3]

! double precision, dimension(:,:), allocatable :: bfrc !!
  double precision, dimension(nxydim, benfdm) :: bfrc !!
! double precision, dimension(:,:,:), allocatable :: sdtx !!
  double precision, dimension(nxydim, sdzdim, sedtdm) :: sdtx !!

  double precision, dimension(:), allocatable :: lon_values !!
  double precision, dimension(:), allocatable :: lat_values !!
  double precision, dimension(:), allocatable :: depth_values !!

  double precision, dimension(:,:,:), allocatable :: outscla !!
  double precision, dimension(:,:,:), allocatable :: outspom !!
  double precision, dimension(:,:,:), allocatable :: outscar !!
  double precision, dimension(:,:,:), allocatable :: outsopa !!
  double precision, dimension(:,:,:), allocatable :: outsage !!
  double precision, dimension(:,:,:), allocatable :: outcnit !!
  double precision, dimension(:,:,:), allocatable :: outcdic !!
  double precision, dimension(:,:,:), allocatable :: outcalk !!
  double precision, dimension(:,:,:), allocatable :: outcsil !!
  double precision, dimension(:,:,:), allocatable :: outcdo2 !!

  double precision, dimension(:,:,:), allocatable :: out13spom !!
  double precision, dimension(:,:,:), allocatable :: out13scar !!
  double precision, dimension(:,:,:), allocatable :: out13cdic !!
  double precision, dimension(:,:,:), allocatable :: out14spom !!
  double precision, dimension(:,:,:), allocatable :: out14scar !!
  double precision, dimension(:,:,:), allocatable :: out14cdic !!


  ! Allocate the array with appropriate dimensions
  allocate(bdep(nxydim))
  allocate(ficla(nxydim))
  allocate(fipom(nxydim))
  allocate(ficar(nxydim))
  allocate(fiopa(nxydim))
  allocate(fiage(nxydim))
  allocate(bwtem(nxydim))
  allocate(bwsal(nxydim))
  allocate(bwnit(nxydim))
  allocate(bwdic(nxydim))
  allocate(bwalk(nxydim))
  allocate(bwdo2(nxydim))
  allocate(bwsil(nxydim))
  allocate(bw13dic(nxydim))
  allocate(bw14dic(nxydim))

  allocate(outscla(nxdim, nydim, sdzdim))
  allocate(outspom(nxdim, nydim, sdzdim))
  allocate(outscar(nxdim, nydim, sdzdim))
  allocate(outsopa(nxdim, nydim, sdzdim))
  allocate(outsage(nxdim, nydim, sdzdim))
  allocate(outcnit(nxdim, nydim, sdzdim))
  allocate(outcdic(nxdim, nydim, sdzdim))
  allocate(outcalk(nxdim, nydim, sdzdim))
  allocate(outcsil(nxdim, nydim, sdzdim))
  allocate(outcdo2(nxdim, nydim, sdzdim))

  ! Call the read_gtool subroutine
  filename = "bound256192/bdep"
  call read_gtool(filename, nxydim, bdep)
  filename = "bound256192/btem"
  call read_gtool(filename, nxydim, bwtem)
  filename = "bound256192/bsal"
  call read_gtool(filename, nxydim, bwsal)
  filename = "bound256192/bnit"
  call read_gtool(filename, nxydim, bwnit)
  filename = "bound256192/bdic"
  call read_gtool(filename, nxydim, bwdic)
  filename = "bound256192/balk"
  call read_gtool(filename, nxydim, bwalk)
  filename = "bound256192/bdo2"
  call read_gtool(filename, nxydim, bwdo2)
  filename = "bound256192/bsil"
  call read_gtool(filename, nxydim, bwsil)
  filename = "bound256192/bc13"
  call read_gtool(filename, nxydim, bw13dic)
  filename = "bound256192/bc14"
  call read_gtool(filename, nxydim, bw14dic)
  filename = "bound256192/ficla"
  call read_gtool(filename, nxydim, ficla)
  filename = "bound256192/fipom"
  call read_gtool(filename, nxydim, fipom)
  filename = "bound256192/ficar"
  call read_gtool(filename, nxydim, ficar)
  filename = "bound256192/fiopa"
  call read_gtool(filename, nxydim, fiopa)

  fiage = 0.d0

  do i = 1, nxydim
    if (bdep(i) == rmiss) then
      omask(i) = 0.d0
    else
      omask(i) = 1.d0
    end if
  end do

  do i = 1, nxydim
    bfrc(i, 1) = ficla(i) * omask(i)
    bfrc(i, 2) = fipom(i) * omask(i)
    bfrc(i, 3) = ficar(i) * omask(i)
    bfrc(i, 4) = fiopa(i) * omask(i)
    bfrc(i, 5) = fiage(i) * omask(i)
    bfrc(i, 6) = bwtem(i)
    bfrc(i, 7) = bwsal(i)
    bfrc(i, 8) = bwnit(i)
    bfrc(i, 9) = bwdic(i)
    bfrc(i, 10) = bwalk(i)
    bfrc(i, 11) = bwsil(i)
    bfrc(i, 12) = bwdo2(i)
    bfrc(i, 13) = bw13dic(i)
    bfrc(i, 14) = bw14dic(i)
  end do

  call digda_cprdc( &
 &     bdep, &
 &     bfrc, &
 &     sdtx )

  do k = 1, sdzdim
    do j = 1, nydim
      do i = 1, nxdim
        l = i + (j - 1) * nxdim
        outscla(i, j, k) = sdtx(l, k, 1)
        outspom(i, j, k) = sdtx(l, k, 2)
        outscar(i, j, k) = sdtx(l, k, 3)
        outsopa(i, j, k) = sdtx(l, k, 4)
        outsage(i, j, k) = sdtx(l, k, 5)
        outcnit(i, j, k) = sdtx(l, k, 6)
        outcdic(i, j, k) = sdtx(l, k, 7)
        outcalk(i, j, k) = sdtx(l, k, 8)
        outcsil(i, j, k) = sdtx(l, k, 9)
        outcdo2(i, j, k) = sdtx(l, k, 10)
!       out13spom(i, j, k) = sdtx(l, k, 11)
!       out13scar(i, j, k) = sdtx(l, k, 12)
!       out13cdic(i, j, k) = sdtx(l, k, 13)
!       out14spom(i, j, k) = sdtx(l, k, 14)
!       out14scar(i, j, k) = sdtx(l, k, 15)
!       out14cdic(i, j, k) = sdtx(l, k, 16)
      end do
    end do
  end do

  ! Call the write_netcdf subroutine
  filename = "scla.nc"
  varitem = "scla"
  vartitle = "clay"
  varunit = "g m-3"
  call write_gtool(filename, varitem, vartitle, varunit, lon_values, lat_values, lz, outscla)

  filename = "spom.nc"
  varitem = "spom"
  vartitle = "organic carbon"
  varunit = "mmol m-3"
  call write_gtool(filename, varitem, vartitle, varunit, lon_values, lat_values, lz, outspom)

  filename = "scar.nc"
  varitem = "scar"
  vartitle = "CaCO3"
  varunit = "mmol m-3"
  call write_gtool(filename, varitem, vartitle, varunit, lon_values, lat_values, lz, outscar)

  filename = "sopa.nc"
  varitem = "sopa"
  vartitle = "SiO2"
  varunit = "mmol m-3"
  call write_gtool(filename, varitem, vartitle, varunit, lon_values, lat_values, lz, outsopa)

  filename = "sage.nc"
  varitem = "sage"
  vartitle = "agetracer"
  varunit = "yr"
  call write_gtool(filename, varitem, vartitle, varunit, lon_values, lat_values, lz, outsage)

  filename = "cnit.nc"
  varitem = "NIT"
  vartitle = "Phosphate"
  varunit = "mmol m-3"
  call write_gtool(filename, varitem, vartitle, varunit, lon_values, lat_values, lz, outcnit)

  filename = "cdic.nc"
  varitem = "DIC"
  vartitle = "DIC"
  varunit = "mmol m-3"
  call write_gtool(filename, varitem, vartitle, varunit, lon_values, lat_values, lz, outcdic)

  filename = "calk.nc"
  varitem = "ALK"
  vartitle = "Alkalinity"
  varunit = "mmol m-3"
  call write_gtool(filename, varitem, vartitle, varunit, lon_values, lat_values, lz, outcalk)

  filename = "cdo2.nc"
  varitem = "O2"
  vartitle = "Oxygen"
  varunit = "mmol m-3"
  call write_gtool(filename, varitem, vartitle, varunit, lon_values, lat_values, lz, outcdo2)

  filename = "csil.nc"
  varitem = "SIL"
  vartitle = "Silicate"
  varunit = "mmol m-3"
  call write_gtool(filename, varitem, vartitle, varunit, lon_values, lat_values, lz, outcsil)

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

end program digda_dataio_gtool
