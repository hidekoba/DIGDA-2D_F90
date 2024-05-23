!######################################################################
module module_digda_timestep
  implicit none

  ! Time step settings
  ! Time_Step_1y
!!double precision, parameter :: TSYR = 365.0D0 / 365D0   !! [yr]
!!! Time_Step_1_2y
!!double precision, parameter :: TSYR = 182.5D0 / 365D0   !! [yr]
!!! Time_Step_50d
!!double precision, parameter :: TSYR = 50.0D0 / 365D0   !! [yr]
!!! Time_Step_30d
!!double precision, parameter :: TSYR = 1.0D0 / 12.0D0   !! [yr]
!!! Time_Step_10d
!!double precision, parameter :: TSYR = 10.0D0 / 365D0   !! [yr]
!!! Time_Step_5d
!!double precision, parameter :: TSYR = 5.0D0 / 365D0   !! [yr]
!!! Time_Step_2d
!!double precision, parameter :: TSYR = 2.0D0 / 365D0   !! [yr]
  ! Time_Step_1d
  double precision, parameter :: TSYR = 1.0D0 / 365D0   !! [yr]
!!! Time_Step_1_10d
!!double precision, parameter :: TSYR = 0.1D0 / 365D0   !! [yr]
!!! Time_Step_1_25d
!!double precision, parameter :: TSYR = 0.04D0 / 365D0   !! [yr]
!!! Time_Step_1_50d
!!double precision, parameter :: TSYR = 0.02D0 / 365D0   !! [yr]

  double precision :: TS = TSYR * 3.1536D+7  !! [s]

  ! -- Start and end year
  integer, parameter :: TSTR = 1
  integer, parameter :: TENDYR = 3
  integer, parameter :: TEND = TENDYR / TSYR !! Converts years to time steps

  ! -- Interval for outputting history files
  integer, parameter :: TOUTYR = 1
  integer, parameter :: TOUT = TOUTYR / TSYR !! Converts years to time steps

end module module_digda_timestep
