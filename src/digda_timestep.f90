!######################################################################
module module_digda_timestep
  use module_zcfile
  implicit none

! ! Time step settings
! ! Time_Step_1y
!!double precision, parameter :: sdtsyr = 365.0d0 / 365d0   !! [yr]
!!! Time_Step_1_2y
!!double precision, parameter :: sdtsyr = 182.5d0 / 365d0   !! [yr]
!!! Time_Step_50d
!!double precision, parameter :: sdtsyr = 50.0d0 / 365d0   !! [yr]
!!! Time_Step_30d
!!double precision, parameter :: sdtsyr = 1.0d0 / 12.0d0   !! [yr]
!!! Time_Step_10d
!!double precision, parameter :: sdtsyr = 10.0d0 / 365d0   !! [yr]
!!! Time_Step_5d
!!double precision, parameter :: sdtsyr = 5.0d0 / 365d0   !! [yr]
!!! Time_Step_2d
!!double precision, parameter :: sdtsyr = 2.0d0 / 365d0   !! [yr]
! ! Time_Step_1d
! double precision, parameter :: sdtsyr = 1.0d0 / 365d0   !! [yr]
!!! Time_Step_1_10d
!!double precision, parameter :: sdtsyr = 0.1d0 / 365d0   !! [yr]
!!! Time_Step_1_25d
!!double precision, parameter :: sdtsyr = 0.04d0 / 365d0   !! [yr]
!!! Time_Step_1_50d
!!double precision, parameter :: sdtsyr = 0.02d0 / 365d0   !! [yr]

! double precision :: ts = sdtsyr * 3.1536d+7  !! [s]

! ! Start and end year
! integer, parameter :: tstr = 1
! integer, parameter :: sdyrend = 500
! integer, parameter :: tend = sdyrend / sdtsyr !! Converts years to time steps
!!namelist /nmsdtm/ sdyrend

! ! Interval for outputting history files
! integer, parameter :: sdyrhst = 100
! integer, parameter :: thst = sdyrhst / sdtsyr !! Converts years to time steps

! ! Control whether history files and time steps are written
! logical :: wotstp = .true.  !! write out time step ?
! logical :: wohist = .true.  !! write out history file ?
! namelist /nmsdhi/ wotstp, wohist


  ! Time step settings
  ! Time_Step_1y
!!double precision :: sdtsyr = 365.0d0 / 365d0   !! [yr]
!!! Time_Step_1_2y
!!double precision :: sdtsyr = 182.5d0 / 365d0   !! [yr]
!!! Time_Step_50d
!!double precision :: sdtsyr = 50.0d0 / 365d0   !! [yr]
!!! Time_Step_30d
!!double precision :: sdtsyr = 1.0d0 / 12.0d0   !! [yr]
!!! Time_Step_10d
!!double precision :: sdtsyr = 10.0d0 / 365d0   !! [yr]
!!! Time_Step_5d
!!double precision :: sdtsyr = 5.0d0 / 365d0   !! [yr]
!!! Time_Step_2d
!!double precision :: sdtsyr = 2.0d0 / 365d0   !! [yr]
  ! Time_Step_1d
  double precision :: sdtsyr = 1.0d0 / 365d0   !! [yr]
!!! Time_Step_1_10d
!!double precision :: sdtsyr = 0.1d0 / 365d0   !! [yr]
!!! Time_Step_1_25d
!!double precision :: sdtsyr = 0.04d0 / 365d0   !! [yr]
!!! Time_Step_1_50d
!!double precision :: sdtsyr = 0.02d0 / 365d0   !! [yr]

  double precision :: tssec !! [s]

  ! Start and end year
  integer :: tstr = 1
  integer :: sdyrstr = 0
  integer :: sdyrend = 500
  integer :: tend
  logical :: wotstp = .true.  !! write out time step ?
  namelist /nmsdtm/ sdtsyr, sdyrstr, sdyrend, wotstp

  ! Interval for outputting history files
  integer :: sdyrhst = 100
  integer :: thst
  logical :: wohist = .true.  !! write out history file ?
  logical :: hstavg = .true.  !! output the average file as history file ?
  namelist /nmsdhi/ sdyrhst, wohist, hstavg

  contains

  subroutine initialize_timestep()

    integer :: ios

    ! Open and read the namelist file
    rewind(nfparm)
    read(nfparm, NML=nmsdtm, IOSTAT=ios)
    write(nfstdo, nmsdtm)
    rewind(nfparm)
    read(nfparm, NML=nmsdhi, IOSTAT=ios)
    write(nfstdo, nmsdhi)

    ! Calculate derived values
    tssec = sdtsyr * 3.1536d+7  !! [yr] --> [s]
    tend = (sdyrend - sdyrstr) / sdtsyr !! Converts years to time steps
    thst = sdyrhst / sdtsyr !! Converts years to time steps

  end subroutine initialize_timestep

end module module_digda_timestep
