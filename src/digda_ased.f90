program digda_ased

  use module_zcfile
  use module_ucloc
  use module_digda_dataio_netcdf
  use module_digda_cprdc
  implicit none

  double precision, dimension(nxydim) :: bdep !! ocean depth [m]

  double precision, dimension(nxydim, benfdm) :: bfrc
  double precision, dimension(nxydim, sdzdim, sedtdm) :: sdtx

  double precision, dimension(:), allocatable :: lon_values
  double precision, dimension(:), allocatable :: lat_values
  double precision, dimension(:), allocatable :: depth_values
  integer(kind=4),  dimension(:), allocatable :: time_values

  integer :: lenstd
  character :: crun*(ncf) = '(RUN NAME WAS NOT SET)'
  character :: cstdo*(ncf) = 'STDOUT_DIGDA'

  namelist /nmrun/ crun
  namelist /nmstdo/ cstdo

  integer :: ios
  integer :: dt(8)

  !-- Initialize MPI --
  if (use_mpi) then
    call initialize_digda_mpi()
  else
    myrank = 0
  end if

  !-- Initial setup --
  call clcstr('SETUP')
  open( UNIT=nfparm, FILE='PARAMET_DIGDA', STATUS='OLD', &
      & ACCESS='SEQUENTIAL', FORM='FORMATTED'      )
  rewind(nfparm)
  read(nfparm, NML=nmstdo, IOSTAT=ios)

  lenstd = index(cstdo, ' ')
  write(cstdo(lenstd:lenstd+2), '(A1,I2.2)') '.', myrank
  open( UNIT=nfstdo, FILE=cstdo, STATUS='UNKNOWN', &
      & ACCESS='SEQUENTIAL', FORM='FORMATTED' )
  write(nfstdo, *) "'--------------------------------------------------'"
  write(nfstdo, *) 'MESSAGE OUTPUT FOR RANK', myrank
  rewind(nfparm)
  read(nfparm, NML=nmrun, IOSTAT=ios)
  write(nfstdo, *) 'Run name :'//trim(crun)
  write(nfstdo, *) "'--------------------------------------------------'"

  call date_and_time(values=dt)
  write(nfstdo, *) "'--------------------------------------------------'"
  write(nfstdo,'(A,I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2)') &
      &            "'-- DIGDA MAIN: Started at ", &
      &              dt(1), dt(2), dt(3), dt(5), dt(6), dt(7)
  write(nfstdo, *) "'-- DIGDA MAIN: Started                          --'"
  write(nfstdo, *) "'--------------------------------------------------'"

! if ( myrank == iroot ) then
  open(UNIT=nfsots, FILE="TSTEP_DIGDA", STATUS="UNKNOWN")
! end if

 !!-- Read the namelist --
 !call read_namelist()

  !-- Initialize the name of tracer, forcing
  call init_digda_grid()

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

  call date_and_time(values=dt)
  write(nfstdo, *) "'--------------------------------------------------'"
  write(nfstdo, *) "'--------------------------------------------------'"
  write(nfstdo,'(A,I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2)') &
      &            "'-- DIGDA MAIN: Finished at ", &
      &              dt(1), dt(2), dt(3), dt(5), dt(6), dt(7)
  write(nfstdo, *) "'-- DIGDA MAIN: Finished                         --'"
  write(nfstdo, *) "'--------------------------------------------------'"

end program digda_ased
