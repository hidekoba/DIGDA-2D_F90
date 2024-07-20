!######################################################################
module module_digda_timestep
  implicit none

  ! Time step settings
  ! Time_Step_1y
!!double precision, parameter :: tsyr = 365.0d0 / 365d0   !! [yr]
!!! Time_Step_1_2y
!!double precision, parameter :: tsyr = 182.5d0 / 365d0   !! [yr]
!!! Time_Step_50d
!!double precision, parameter :: tsyr = 50.0d0 / 365d0   !! [yr]
!!! Time_Step_30d
!!double precision, parameter :: tsyr = 1.0d0 / 12.0d0   !! [yr]
!!! Time_Step_10d
!!double precision, parameter :: tsyr = 10.0d0 / 365d0   !! [yr]
!!! Time_Step_5d
!!double precision, parameter :: tsyr = 5.0d0 / 365d0   !! [yr]
!!! Time_Step_2d
!!double precision, parameter :: tsyr = 2.0d0 / 365d0   !! [yr]
  ! Time_Step_1d
  double precision, parameter :: tsyr = 1.0d0 / 365d0   !! [yr]
!!! Time_Step_1_10d
!!double precision, parameter :: tsyr = 0.1d0 / 365d0   !! [yr]
!!! Time_Step_1_25d
!!double precision, parameter :: tsyr = 0.04d0 / 365d0   !! [yr]
!!! Time_Step_1_50d
!!double precision, parameter :: tsyr = 0.02d0 / 365d0   !! [yr]

  double precision :: ts = tsyr * 3.1536d+7  !! [s]

  ! Start and end year
  integer, parameter :: tstr = 1
  integer, parameter :: tendyr = 2
  integer, parameter :: tend = tendyr / tsyr !! Converts years to time steps
 !namelist /nmtime/ tendyr

  ! Interval for outputting history files
  integer, parameter :: toutyr = 100
  integer, parameter :: tout = toutyr / tsyr !! Converts years to time steps

  ! Control whether history files and time steps are written
  logical :: wotstp = .true.  !! write out time step ?
  logical :: wohist = .true.  !! write out history file ?
  namelist /nmsdhi/ wotstp, wohist

end module module_digda_timestep
