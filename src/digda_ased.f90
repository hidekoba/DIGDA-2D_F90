program digda_ased

  use module_zcfile
  use module_ucloc
  use module_digda_dataio_netcdf
  use module_digda_cprdc
  implicit none

! character(len=100) :: filename, varname, varunit

! integer :: i, j, k, l

  double precision, dimension(nxydim) :: bdep !! ocean depth [m]

  double precision, dimension(nxydim, benfdm) :: bfrc
  double precision, dimension(nxydim, sdzdim, sedtdm) :: sdtx

  integer :: lenstd
  character :: crun*(ncf) = '(RUN NAME WAS NOT SET)'
  character :: cstdo*(ncf) = 'STDOUT'

  namelist /nmrun/ crun
  namelist /nmstdo/ cstdo

  integer :: ios

  !-- Initialize MPI --
  if (use_mpi) then
    call initialize_digda_mpi()
  else
    myrank = 0
  end if

  !-- Initial setup --
  call clcstr('SETUP')
  open( UNIT=nfparm, FILE='PARAMET', STATUS='OLD', &
      & ACCESS='SEQUENTIAL', FORM='FORMATTED'      )
  rewind(nfparm)
  read(nfparm, NML=nmstdo, IOSTAT=ios)
 !if (ios /= 0) then
 !  write(*, *) 'End of file or error reading NAMELIST file 01'
 !end if

  lenstd = index(cstdo, ' ')
  write(cstdo(lenstd:lenstd+2), '(A1,I2.2)') '.', myrank
  open( UNIT=nfstdo, FILE=cstdo, STATUS='UNKNOWN', &
      & ACCESS='SEQUENTIAL', FORM='FORMATTED' )
  write(nfstdo, *) "'--------------------------------------------------'"
  write(nfstdo, *) 'MESSAGE OUTPUT FOR RANK', myrank
  rewind(nfparm)
  read(nfparm, NML=nmrun, IOSTAT=ios)
  write(nfstdo, *) 'Run name :'//trim(crun)
 !if (ios /= 0) then
 !  write(*, *) 'End of file or error reading NAMELIST file 02'
 !else
 !  write(nfstdo, *) 'Run name :'//trim(crun)
 !end if
  write(nfstdo, *) "'--------------------------------------------------'"

  write(nfstdo, *) "'--------------------------------------------------'"
  write(nfstdo, *) "'-- DIGDA MAIN: Started                          --'"
  write(nfstdo, *) "'--------------------------------------------------'"

! open(nfparm, FILE="PARAMET", STATUS="UNKNOWN")
! open(nfstdo, FILE="STDOUT",  STATUS="UNKNOWN")
! if ( myrank == iroot ) then
  open(UNIT=nfsots, FILE="TSTEP", STATUS="UNKNOWN")
! end if

 !!-- Read the namelist --
 !call read_namelist()

  !-- Read RESTART file --
  call digda_read_restart(sdtx, lon_values, lat_values, depth_values, time_values)

  !-- Read the forcing --
  call digda_read_bottom_forcing(bdep, bfrc)

  call clcend('SETUP')

  !-- Start time-stepping --
  call clcstr('PREDIC')
  call digda_cprdc(bdep, bfrc, sdtx, lon_values, lat_values, depth_values)
  call clcend('PREDIC')

  !-- Write the restart file --
  call clcstr('FINOUT')
  call digda_write_restart(sdtx, lon_values, lat_values, depth_values, time_values)
  call digda_write_data()
  call clcend('FINOUT')
  call clcout()

  !-- Finalize MPI
  if (use_mpi) then
    call finalize_digda_mpi()
  end if

  write(nfstdo, *) "'--------------------------------------------------'"
  write(nfstdo, *) "'-- DIGDA MAIN: Finished                         --'"
  write(nfstdo, *) "'--------------------------------------------------'"

contains

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

end program digda_ased
