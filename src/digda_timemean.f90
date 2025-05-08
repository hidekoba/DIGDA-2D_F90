!######################################################################

! --- information -----------------------------------------------------
!
!  Sediment model : Calculate time mean value
!
!  HISTORY
!     '15.01.09  H.Kobayashi :
!
! ---------------------------------------------------------------------

module module_digda_time_mean
    use module_digda_grid
    implicit none
    private
    public :: time_mean_type, initialize_time_mean, update_time_mean, get_time_mean, reset_time_mean

    type :: time_mean_type
        double precision, allocatable :: sum(:)       ! Accumulated sum
        integer, allocatable :: count(:)  ! Count of values
        integer :: nxy, nz                ! Grid dimensions
        integer :: ndim                   ! Number of dimensions (2 or 3)
    end type time_mean_type

contains

    subroutine initialize_time_mean(mean_data, nxy, nz, ndim)
        ! Initialize the time mean structure
        type(time_mean_type), intent(inout) :: mean_data
        integer, intent(in) :: nxy, nz, ndim

        mean_data%nxy = nxy
        mean_data%nz = nz
        mean_data%ndim = ndim

        if (ndim == 2) then
            allocate(mean_data%sum(nxy))
            allocate(mean_data%count(nxy))
        else
            allocate(mean_data%sum(nxy * nz))
            allocate(mean_data%count(nxy * nz))
        end if
        ! Reset the accumulated data
        call reset_time_mean(mean_data)
    end subroutine initialize_time_mean

    subroutine update_time_mean(mean_data, data)
        ! Update the time mean with new data
        type(time_mean_type), intent(inout) :: mean_data
        double precision, intent(in) :: data(:,:)

        mean_data%sum = mean_data%sum + reshape(data, shape(mean_data%sum))
        mean_data%count = mean_data%count + 1
    end subroutine update_time_mean

    subroutine get_time_mean(mean_data, avg_data)
        ! Compute and retrieve the time mean for the grid
        type(time_mean_type), intent(in) :: mean_data
        double precision, intent(out) :: avg_data(:,:)

        if (any(mean_data%count == 0)) then
            avg_data = 0.0  ! Avoid division by zero
        else
            if (mean_data%ndim == 2) then
                avg_data(:,1) = reshape(mean_data%sum / real(mean_data%count, kind=kind(1.0)), &
                                          [mean_data%nxy])
            else
                avg_data = reshape(mean_data%sum / real(mean_data%count, kind=kind(1.0)), &
                                   [mean_data%nxy, mean_data%nz])
            end if
        end if
    end subroutine get_time_mean

    subroutine reset_time_mean(mean_data)
        ! Reset the accumulated sum and count for the grid
        type(time_mean_type), intent(inout) :: mean_data

        mean_data%sum = 0.0
        mean_data%count = 0
    end subroutine reset_time_mean

end module module_digda_time_mean
