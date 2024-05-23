!######################################################################
module module_digda_config
    implicit none

    ! Parallel digda_configuration
    logical, parameter :: PARALLEL = .TRUE.

    ! Grid digda_configuration
    logical, parameter :: ONE_DIM = .FALSE.
    logical, parameter :: Grid_128x64 = .FALSE.
    logical, parameter :: Grid_256x192 = .TRUE.
    logical, parameter :: Grid_360x180 = .FALSE.
    logical, parameter :: Grid_2160x1080 = .FALSE.

    ! Vertical grid digda_configuration
    logical, parameter :: ZGRID_CONSTANT = .FALSE.
    logical, parameter :: ZGRID_ARCHER = .FALSE.
    logical, parameter :: ZGRID_HEINZE = .TRUE.

    ! Time step digda_configuration
    logical, parameter :: Time_Step_1y = .FALSE.
    logical, parameter :: Time_Step_1_2y = .FALSE.
    logical, parameter :: Time_Step_50d = .FALSE.
    logical, parameter :: Time_Step_30d = .FALSE.
    logical, parameter :: Time_Step_10d = .FALSE.
    logical, parameter :: Time_Step_5d = .FALSE.
    logical, parameter :: Time_Step_2d = .FALSE.
    logical, parameter :: Time_Step_1d = .TRUE.
    logical, parameter :: Time_Step_1_10d = .FALSE.
    logical, parameter :: Time_Step_1_25d = .FALSE.
    logical, parameter :: Time_Step_1_50d = .FALSE.

end module module_digda_config
