module netcdf_utilities
    use netcdf
    implicit none
contains

  subroutine read_netcdf_2dwot(filename, varname, var, xvalues, yvalues)
    implicit none

    character(len=*), intent(in) :: filename, varname
    double precision, dimension(:,:), allocatable, intent(out) :: var
    double precision, dimension(:), allocatable, intent(out) :: xvalues, yvalues

    integer(kind=4) :: ncid, varid, retval
    integer(kind=4) :: xid, xsize
    integer(kind=4) :: yid, ysize

    character(len=*), parameter :: xname = 'lon'
    character(len=*), parameter :: yname = 'lat'

    integer(kind=4), dimension(1) :: start_axis, count_axis
    integer(kind=4), dimension(2) :: start_var, count_var
    integer(kind=4) :: i

    ! Open the netCDF file
    retval = nf90_open(filename, nf90_nowrite, ncid)
    if (retval /= nf90_noerr) then
      print *, 'Failed to open file: ', trim(filename)
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Retrieve the dimension IDs and sizes for longitude, and latitude
    ! Get the longitude dimension ID and size
    retval = nf90_inq_dimid(ncid, xname, xid)
    if (retval == nf90_noerr) then
      retval = nf90_inquire_dimension(ncid, xid, len=xsize)
    else
      print *, 'Failed to get ID and length of longitude dimension.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif
    print *, 'Longitude size = ', xsize

    ! Get the latitude dimension ID and size
    retval = nf90_inq_dimid(ncid, yname, yid)
    if (retval == nf90_noerr) then
      retval = nf90_inquire_dimension(ncid, yid, len=ysize)
    else
      print *, 'Failed to get ID and length of latitude dimension.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif
    print *, 'Latitude size = ', ysize

    ! Deallocate arrays if they are already allocated
    if (allocated(xvalues)) deallocate(xvalues)
    if (allocated(yvalues)) deallocate(yvalues)
    if (allocated(var)) deallocate(var)

    ! Allocate memory for the dimension values
    allocate(xvalues(xsize))
    allocate(yvalues(ysize))

!   ! Retrieve the longitude values
!   start_axis = (/ 1 /)
!   count_axis = (/ xsize /)
!   retval = nf90_get_var(ncid, xid, xvalues, start=start_axis, count=count_axis)
!   if (retval /= nf90_noerr) then
!    !print *, 'Failed to read longitude values.'
!     print *, 'Failed to read longitude values.', varname
!     print *, 'Error: ', nf90_strerror(retval)
!     stop
!   endif
!  !print *, 'Longitude values:'
!  !do i = 1, xsize
!  !  print *, xvalues(i)
!  !end do

!   ! Retrieve the latitude values
!   start_axis = (/ 1 /)
!   count_axis = (/ ysize /)
!   retval = nf90_get_var(ncid, yid, yvalues, start=start_axis, count=count_axis)
!   if (retval /= nf90_noerr) then
!     print *, 'Failed to read latitude values.'
!     print *, 'Error: ', nf90_strerror(retval)
!     stop
!   endif
!  !print *, 'Latitude values:'
!  !do i = 1, ysize
!  !  print *, yvalues(i)
!  !end do

    ! Allocate memory for the 2D data variable based on retrieved dimensions
    allocate(var(xsize, ysize))
    if (.not. allocated(var)) then
      print *, 'Allocation of ', varname, ' failed.'
      stop
    endif

    ! Retrieve the variable ID and read the data into the allocated array
    retval = nf90_inq_varid(ncid, varname, varid)
    if (retval /= nf90_noerr) then
      print *, 'Failed to retrieve variable ID: ', varname
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif
    start_var = (/ 1, 1 /)
    count_var = (/ xsize, ysize /)
    retval = nf90_get_var(ncid, varid, var, start=start_var, count=count_var)
    if (retval /= nf90_noerr) then
      print *, 'Failed to read data for variable: ', varname
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Close the netCDF file
    retval = nf90_close(ncid)
    if (retval /= nf90_noerr) then
      print *, 'Error occurred while closing the file: ', trim(filename)
      stop
    endif

   !! Deallocate arrays if they are already allocated
   !if (allocated(xvalues)) deallocate(xvalues)
   !if (allocated(yvalues)) deallocate(yvalues)
   !if (allocated(var)) deallocate(var)

  end subroutine read_netcdf_2dwot


  subroutine read_netcdf_2dwt(filename, varname, var, xvalues, yvalues, tvalues)
    implicit none

    character(len=*), intent(in) :: filename, varname
    double precision, dimension(:,:,:), allocatable, intent(out) :: var
    double precision, dimension(:), allocatable, intent(out) :: xvalues, yvalues
    integer(kind=4), dimension(:), allocatable, intent(out) :: tvalues

    integer(kind=4) :: ncid, varid, retval
    integer(kind=4) :: xid, xsize
    integer(kind=4) :: yid, ysize
    integer(kind=4) :: tid, tsize

    character(len=*), parameter :: xname = 'lon'
    character(len=*), parameter :: yname = 'lat'
    character(len=*), parameter :: tname = 'time'

    integer(kind=4), dimension(1) :: start_axis, count_axis
    integer(kind=4), dimension(3) :: start_var, count_var

    ! Open the netCDF file
    retval = nf90_open(filename, nf90_nowrite, ncid)
    if (retval /= nf90_noerr) then
      print *, 'Failed to open file: ', trim(filename)
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Retrieve the dimension IDs and sizes for longitude, latitude, and depth
    ! Get the longitude dimension ID and size
    retval = nf90_inq_dimid(ncid, xname, xid)
    if (retval == nf90_noerr) then
      retval = nf90_inquire_dimension(ncid, xid, len=xsize)
    else
      print *, 'Failed to get ID and length of longitude dimension.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif
    print *, 'Longitude size = ', xsize

    ! Get the latitude dimension ID and size
    retval = nf90_inq_dimid(ncid, yname, yid)
    if (retval == nf90_noerr) then
      retval = nf90_inquire_dimension(ncid, yid, len=ysize)
    else
      print *, 'Failed to get ID and length of latitude dimension.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif
    print *, 'Latitude size = ', ysize

    ! Get the time dimension ID and size
    retval = nf90_inq_dimid(ncid, tname, tid)
    if (retval == nf90_noerr) then
      retval = nf90_inquire_dimension(ncid, tid, len=tsize)
    else
      print *, 'Failed to get ID and length of time dimension.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif
    print *, 'time size = ', tsize

    ! Deallocate arrays if they are already allocated
    if (allocated(xvalues)) deallocate(xvalues)
    if (allocated(yvalues)) deallocate(yvalues)
    if (allocated(tvalues)) deallocate(tvalues)
    if (allocated(var)) deallocate(var)

    ! Allocate memory for the dimension values
    allocate(xvalues(xsize))
    allocate(yvalues(ysize))
    allocate(tvalues(tsize))

    ! Retrieve the longitude values
    start_axis = (/ 1 /)
    count_axis = (/ xsize /)
    retval = nf90_get_var(ncid, xid, xvalues, start=start_axis, count=count_axis)
    if (retval /= nf90_noerr) then
     !print *, 'Failed to read longitude values.'
      print *, 'Failed to read longitude values.', varname
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Retrieve the latitude values
    start_axis = (/ 1 /)
    count_axis = (/ ysize /)
    retval = nf90_get_var(ncid, yid, yvalues, start=start_axis, count=count_axis)
    if (retval /= nf90_noerr) then
      print *, 'Failed to read latitude values.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Retrieve the time values
    start_axis = (/ 1 /)
    count_axis = (/ tsize /)
    retval = nf90_get_var(ncid, tid, tvalues, start=start_axis, count=count_axis)
    if (retval /= nf90_noerr) then
      print *, 'Failed to read time values.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Allocate memory for the 4D data variable based on retrieved dimensions
    allocate(var(xsize, ysize, tsize))
    if (.not. allocated(var)) then
      print *, 'Allocation of ', varname, ' failed.'
      stop
    endif

    ! Retrieve the variable ID and read the data into the allocated array
    retval = nf90_inq_varid(ncid, varname, varid)
    if (retval /= nf90_noerr) then
      print *, 'Failed to retrieve variable ID: ', varname
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif
    start_var = (/ 1, 1, 1 /)
    count_var = (/ xsize, ysize, tsize /)
    retval = nf90_get_var(ncid, varid, var, start=start_var, count=count_var)
    if (retval /= nf90_noerr) then
      print *, 'Failed to read data for variable: ', varname
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Close the netCDF file
    retval = nf90_close(ncid)
    if (retval /= nf90_noerr) then
      print *, 'Error occurred while closing the file: ', trim(filename)
      stop
    endif

   !! Deallocate arrays if they are already allocated
   !if (allocated(xvalues)) deallocate(xvalues)
   !if (allocated(yvalues)) deallocate(yvalues)
   !if (allocated(zvalues)) deallocate(zvalues)
   !if (allocated(tvalues)) deallocate(tvalues)
   !if (allocated(var)) deallocate(var)

  end subroutine read_netcdf_2dwt


  subroutine read_netcdf_3dwot(filename, varname, var, xvalues, yvalues, zvalues)
    implicit none

    character(len=*), intent(in) :: filename, varname
    double precision, dimension(:,:,:), allocatable, intent(out) :: var
    double precision, dimension(:), allocatable, intent(out) :: xvalues, yvalues, zvalues

    integer(kind=4) :: ncid, varid, retval
    integer(kind=4) :: xid, xsize
    integer(kind=4) :: yid, ysize
    integer(kind=4) :: zid, zsize

    character(len=*), parameter :: xname = 'lon'
    character(len=*), parameter :: yname = 'lat'
    character(len=*), parameter :: zname = 'depth'

    integer(kind=4), dimension(1) :: start_axis, count_axis
    integer(kind=4), dimension(3) :: start_var, count_var

    ! Open the netCDF file
    retval = nf90_open(filename, nf90_nowrite, ncid)
    if (retval /= nf90_noerr) then
      print *, 'Failed to open file: ', trim(filename)
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Retrieve the dimension IDs and sizes for longitude, latitude, and depth
    ! Get the longitude dimension ID and size
    retval = nf90_inq_dimid(ncid, xname, xid)
    if (retval == nf90_noerr) then
      retval = nf90_inquire_dimension(ncid, xid, len=xsize)
    else
      print *, 'Failed to get ID and length of longitude dimension.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif
   !print *, 'Longitude size = ', xsize

    ! Get the latitude dimension ID and size
    retval = nf90_inq_dimid(ncid, yname, yid)
    if (retval == nf90_noerr) then
      retval = nf90_inquire_dimension(ncid, yid, len=ysize)
    else
      print *, 'Failed to get ID and length of latitude dimension.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif
   !print *, 'Latitude size = ', ysize

    ! Get the depth dimension ID and size
    retval = nf90_inq_dimid(ncid, zname, zid)
    if (retval == nf90_noerr) then
      retval = nf90_inquire_dimension(ncid, zid, len=zsize)
    else
      print *, 'Failed to get ID and length of depth dimension.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif
   !print *, 'Depth size = ', zsize

    ! Deallocate arrays if they are already allocated
    if (allocated(xvalues)) deallocate(xvalues)
    if (allocated(yvalues)) deallocate(yvalues)
    if (allocated(zvalues)) deallocate(zvalues)
    if (allocated(var)) deallocate(var)

    ! Allocate memory for the dimension values
    allocate(xvalues(xsize))
    allocate(yvalues(ysize))
    allocate(zvalues(zsize))

    ! Retrieve the longitude values
    start_axis = (/ 1 /)
    count_axis = (/ xsize /)
    retval = nf90_get_var(ncid, xid, xvalues, start=start_axis, count=count_axis)
    if (retval /= nf90_noerr) then
      print *, 'Failed to read longitude values.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Retrieve the latitude values
    start_axis = (/ 1 /)
    count_axis = (/ ysize /)
    retval = nf90_get_var(ncid, yid, yvalues, start=start_axis, count=count_axis)
    if (retval /= nf90_noerr) then
      print *, 'Failed to read latitude values.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Retrieve the depth values
    start_axis = (/ 1 /)
    count_axis = (/ zsize /)
    retval = nf90_get_var(ncid, zid, zvalues, start=start_axis, count=count_axis)
    if (retval /= nf90_noerr) then
      print *, 'Failed to read depth values.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Allocate memory for the 3D data variable based on retrieved dimensions
    allocate(var(xsize, ysize, zsize))
    if (.not. allocated(var)) then
      print *, 'Allocation of ', varname, ' failed.'
      stop
    endif

    ! Retrieve the variable ID and read the data into the allocated array
    retval = nf90_inq_varid(ncid, varname, varid)
    if (retval /= nf90_noerr) then
      print *, 'Failed to retrieve variable ID: ', varname
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif
    start_var = (/ 1, 1, 1 /)
    count_var = (/ xsize, ysize, zsize /)
    retval = nf90_get_var(ncid, varid, var, start=start_var, count=count_var)
    if (retval /= nf90_noerr) then
      print *, 'Failed to read data for variable: ', varname
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Close the netCDF file
    retval = nf90_close(ncid)
    if (retval /= nf90_noerr) then
      print *, 'Error occurred while closing the file: ', trim(filename)
      stop
    endif

   !! Deallocate arrays if they are already allocated
   !if (allocated(xvalues)) deallocate(xvalues)
   !if (allocated(yvalues)) deallocate(yvalues)
   !if (allocated(zvalues)) deallocate(zvalues)
   !if (allocated(var)) deallocate(var)

  end subroutine read_netcdf_3dwot


  subroutine read_netcdf_3dwt(filename, varname, var, xvalues, yvalues, zvalues, tvalues)
    implicit none

    character(len=*), intent(in) :: filename, varname
    double precision, dimension(:,:,:,:), allocatable, intent(out) :: var
    double precision, dimension(:), allocatable, intent(out) :: xvalues, yvalues, zvalues
    integer(kind=4), dimension(:), allocatable, intent(out) :: tvalues

    integer(kind=4) :: ncid, varid, retval
    integer(kind=4) :: xid, xsize
    integer(kind=4) :: yid, ysize
    integer(kind=4) :: zid, zsize
    integer(kind=4) :: tid, tsize

    character(len=*), parameter :: xname = 'lon'
    character(len=*), parameter :: yname = 'lat'
    character(len=*), parameter :: zname = 'depth'
    character(len=*), parameter :: tname = 'time'

    integer(kind=4), dimension(1) :: start_axis, count_axis
    integer(kind=4), dimension(4) :: start_var, count_var

    ! Open the netCDF file
    retval = nf90_open(filename, nf90_nowrite, ncid)
    if (retval /= nf90_noerr) then
      print *, 'Failed to open file: ', trim(filename)
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Retrieve the dimension IDs and sizes for longitude, latitude, and depth
    ! Get the longitude dimension ID and size
    retval = nf90_inq_dimid(ncid, xname, xid)
    if (retval == nf90_noerr) then
      retval = nf90_inquire_dimension(ncid, xid, len=xsize)
    else
      print *, 'Failed to get ID and length of longitude dimension.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif
   !print *, 'Longitude size = ', xsize

    ! Get the latitude dimension ID and size
    retval = nf90_inq_dimid(ncid, yname, yid)
    if (retval == nf90_noerr) then
      retval = nf90_inquire_dimension(ncid, yid, len=ysize)
    else
      print *, 'Failed to get ID and length of latitude dimension.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif
   !print *, 'Latitude size = ', ysize

    ! Get the depth dimension ID and size
    retval = nf90_inq_dimid(ncid, zname, zid)
    if (retval == nf90_noerr) then
      retval = nf90_inquire_dimension(ncid, zid, len=zsize)
    else
      print *, 'Failed to get ID and length of depth dimension.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif
   !print *, 'Depth size = ', zsize

    ! Get the time dimension ID and size
    retval = nf90_inq_dimid(ncid, tname, tid)
    if (retval == nf90_noerr) then
      retval = nf90_inquire_dimension(ncid, tid, len=tsize)
    else
      print *, 'Failed to get ID and length of time dimension.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif
   !print *, 'time size = ', tsize

    ! Deallocate arrays if they are already allocated
    if (allocated(xvalues)) deallocate(xvalues)
    if (allocated(yvalues)) deallocate(yvalues)
    if (allocated(zvalues)) deallocate(zvalues)
    if (allocated(tvalues)) deallocate(tvalues)
    if (allocated(var)) deallocate(var)

    ! Allocate memory for the dimension values
    allocate(xvalues(xsize))
    allocate(yvalues(ysize))
    allocate(zvalues(zsize))
    allocate(tvalues(tsize))

    ! Retrieve the longitude values
    start_axis = (/ 1 /)
    count_axis = (/ xsize /)
    retval = nf90_get_var(ncid, xid, xvalues, start=start_axis, count=count_axis)
    if (retval /= nf90_noerr) then
      print *, 'Failed to read longitude values.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Retrieve the latitude values
    start_axis = (/ 1 /)
    count_axis = (/ ysize /)
    retval = nf90_get_var(ncid, yid, yvalues, start=start_axis, count=count_axis)
    if (retval /= nf90_noerr) then
      print *, 'Failed to read latitude values.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Retrieve the depth values
    start_axis = (/ 1 /)
    count_axis = (/ zsize /)
    retval = nf90_get_var(ncid, zid, zvalues, start=start_axis, count=count_axis)
    if (retval /= nf90_noerr) then
      print *, 'Failed to read depth values.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Retrieve the time values
    start_axis = (/ 1 /)
    count_axis = (/ tsize /)
    retval = nf90_get_var(ncid, tid, tvalues, start=start_axis, count=count_axis)
    if (retval /= nf90_noerr) then
      print *, 'Failed to read time values.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Allocate memory for the 4D data variable based on retrieved dimensions
    allocate(var(xsize, ysize, zsize, tsize))
    if (.not. allocated(var)) then
      print *, 'Allocation of ', varname, ' failed.'
      stop
    endif

    ! Retrieve the variable ID and read the data into the allocated array
    retval = nf90_inq_varid(ncid, varname, varid)
    if (retval /= nf90_noerr) then
      print *, 'Failed to retrieve variable ID: ', varname
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif
    start_var = (/ 1, 1, 1, 1 /)
    count_var = (/ xsize, ysize, zsize, tsize /)
    retval = nf90_get_var(ncid, varid, var, start=start_var, count=count_var)
    if (retval /= nf90_noerr) then
      print *, 'Failed to read data for variable: ', varname
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Close the netCDF file
    retval = nf90_close(ncid)
    if (retval /= nf90_noerr) then
      print *, 'Error occurred while closing the file: ', trim(filename)
      stop
    endif

   !! Deallocate arrays if they are already allocated
   !if (allocated(xvalues)) deallocate(xvalues)
   !if (allocated(yvalues)) deallocate(yvalues)
   !if (allocated(zvalues)) deallocate(zvalues)
   !if (allocated(tvalues)) deallocate(tvalues)
   !if (allocated(var)) deallocate(var)

  end subroutine read_netcdf_3dwt

  subroutine read_netcdf_auto(filename, varname, var, xvalues, yvalues, zvalues, tvalues)
    use netcdf
    implicit none

    character(len=*), intent(in) :: filename, varname
    double precision, allocatable, intent(out) :: var(:,:,:,:)
    double precision, allocatable, intent(out) :: xvalues(:), yvalues(:), zvalues(:)
    integer(kind=4),  allocatable, intent(out) :: tvalues(:)

    integer(kind=4) :: ncid, varid, retval
   !integer(kind=4) :: dimids(nf90_max_var_dims)
    integer(kind=4), allocatable :: dimids(:)
    integer(kind=4) :: ndims
    integer(kind=4) :: xid, yid, zid, tid
    integer(kind=4) :: xsize, ysize, zsize, tsize
    integer(kind=4), dimension(1) :: start_axis, count_axis
    integer(kind=4), dimension(4) :: start_var, count_var
    integer(kind=4) :: i
    logical :: has_lon, has_lat, has_depth, has_time

    character(len=*), parameter :: xname = "lon"
    character(len=*), parameter :: yname = "lat"
    character(len=*), parameter :: zname = "depth"
    character(len=*), parameter :: tname = "time"

    integer(kind=4) :: nvars
    character(len=192) :: varname_file

    ! Initialize variables
    has_lon = .false.
    has_lat = .false.
    has_depth = .false.
    has_time = .false.
    xsize = 1
    ysize = 1
    zsize = 1
    tsize = 1
    ndims = 0
    nvars = 0

    ! Deallocate arrays if they are already allocated
    if (allocated(dimids)) deallocate(dimids)
    if (allocated(xvalues)) deallocate(xvalues)
    if (allocated(yvalues)) deallocate(yvalues)
    if (allocated(zvalues)) deallocate(zvalues)
    if (allocated(tvalues)) deallocate(tvalues)
    if (allocated(var)) deallocate(var)

    ! Open the NetCDF file
    retval = nf90_open(filename, nf90_nowrite, ncid)
    if (retval /= nf90_noerr) then
      print *, 'Failed to open file: ', trim(filename)
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Retrieve the variable ID and read the data into the allocated array
    retval = nf90_inq_varid(ncid, varname, varid)
    if (retval /= nf90_noerr) then
      print *, 'Failed to retrieve variable ID: ', varname
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Inquire the dimensions of the variable
    allocate(dimids(nf90_max_var_dims))  ! Allocate enough space for dimids
    retval = nf90_inquire_variable(ncid, varid, ndims=ndims, dimids=dimids)
    if (retval /= nf90_noerr) then
      print *, 'Failed to inquire variable: ', varname
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Check each dimension and determine its size
    do i = 1, ndims
      select case (i)
        case (1)  ! Longitude (x)
          retval = nf90_inq_dimid(ncid, xname, xid)
          if (retval == nf90_noerr) then
            retval = nf90_inquire_dimension(ncid, xid, len=xsize)
            has_lon = .true.
          endif
        case (2)  ! Latitude (y)
          retval = nf90_inq_dimid(ncid, yname, yid)
          if (retval == nf90_noerr) then
            retval = nf90_inquire_dimension(ncid, yid, len=ysize)
            has_lat = .true.
          endif
        case (3)  ! Depth or Time
          retval = nf90_inq_dimid(ncid, zname, zid)
          if (retval == nf90_noerr) then
            retval = nf90_inquire_dimension(ncid, zid, len=zsize)
            has_depth = .true.
          else
            retval = nf90_inq_dimid(ncid, tname, tid)
            if (retval == nf90_noerr) then
              retval = nf90_inquire_dimension(ncid, tid, len=tsize)
              has_time = .true.
            endif
          endif
        case (4)  ! Time (if it's the fourth dimension)
          retval = nf90_inq_dimid(ncid, tname, tid)
          if (retval == nf90_noerr) then
            retval = nf90_inquire_dimension(ncid, tid, len=tsize)
            has_time = .true.
          endif
      end select
    end do

    ! Allocate arrays based on the dimensions
    if (has_lon) allocate(xvalues(xsize))
    if (has_lat) allocate(yvalues(ysize))
    if (has_depth) allocate(zvalues(zsize))
    if (has_time) allocate(tvalues(tsize))
    allocate(var(max(1, xsize), max(1, ysize), max(1, zsize), max(1, tsize)))

    ! Get dimension IDs to read the dimension values
    retval = nf90_inquire(ncid, ndims, nvars)
    if (retval /= nf90_noerr) then
      print *, "Error inquiring file: ", nf90_strerror(retval)
      stop
    endif

    do i = 1, nvars
      retval = nf90_inquire_variable(ncid, i, varname_file)
      if (retval /= nf90_noerr) then
        print *, "Error reading variable name at index ", i, nf90_strerror(retval)
        cycle
      endif
      if (trim(varname_file) == trim(xname)) then
        xid = i
      endif
      if (trim(varname_file) == trim(yname)) then
        yid = i
      endif
      if (trim(varname_file) == trim(zname)) then
        zid = i
      endif
      if (trim(varname_file) == trim(tname)) then
        tid = i
      endif
    end do

    ! Read longitude values if available
    if (has_lon) then
      start_axis = (/ 1 /)
      count_axis = (/ xsize /)
      retval = nf90_get_var(ncid, xid, xvalues, start=start_axis, count=count_axis)
      if (retval /= nf90_noerr) then
        print *, 'Failed to read longitude values.'
        print *, nf90_strerror(retval)
        stop
      endif
    endif

    ! Read latitude values if available
    if (has_lat) then
      start_axis = (/ 1 /)
      count_axis = (/ ysize /)
      retval = nf90_get_var(ncid, yid, yvalues, start=start_axis, count=count_axis)
      if (retval /= nf90_noerr) then
        print *, 'Failed to read latitude values.'
        print *, nf90_strerror(retval)
        stop
      endif
    endif

    ! Read depth values if available
    if (has_depth) then
      start_axis = (/ 1 /)
      count_axis = (/ zsize /)
      retval = nf90_get_var(ncid, zid, zvalues, start=start_axis, count=count_axis)
      if (retval /= nf90_noerr) then
        print *, 'Failed to read depth values.'
        print *, nf90_strerror(retval)
        stop
      endif
    endif

    ! Read time values if available
    if (has_time) then
      start_axis = (/ 1 /)
      count_axis = (/ tsize /)
      retval = nf90_get_var(ncid, tid, tvalues, start=start_axis, count=count_axis)
      if (retval /= nf90_noerr) then
        print *, 'Failed to read time values.'
        print *, nf90_strerror(retval)
        stop
      endif
    endif

    ! Read variable data
    retval = nf90_inq_varid(ncid, varname, varid)
    if (retval /= nf90_noerr) then
      print *, 'Failed to retrieve variable ID: ', varname
      print *, nf90_strerror(retval)
      stop
    endif
    start_var = (/ 1, 1, 1, 1 /)
    count_var = (/ xsize, ysize, zsize, tsize /)
    retval = nf90_get_var(ncid, varid, var, start=start_var, count=count_var)
    if (retval /= nf90_noerr) then
      print *, 'Failed to read variable data.'
      print *, nf90_strerror(retval)
      stop
    endif

    ! Close the NetCDF file
    retval = nf90_close(ncid)
    if (retval /= nf90_noerr) then
      print *, 'Failed to close file: ', trim(filename)
      print *, nf90_strerror(retval)
      stop
    endif

  end subroutine read_netcdf_auto

  subroutine write_netcdf_2dwot(filename, varname, varunit, var, xvalues, yvalues)
    implicit none

    character(len=*), intent(in) :: filename, varname, varunit
    double precision, dimension(:,:), intent(in) :: var
    double precision, dimension(:), intent(in) :: xvalues, yvalues

    integer(kind=4) :: ncid, varid, retval
    integer(kind=4) :: xid, xsize
    integer(kind=4) :: yid, ysize
    integer(kind=4) :: xvarid, yvarid

    character(len=*), parameter :: xname = 'lon'
    character(len=*), parameter :: yname = 'lat'

    double precision, parameter :: fill_value = 999.0

    ! Determine the size of each dimension from the input arrays
    xsize = size(xvalues)
    ysize = size(yvalues)

    ! Create the netCDF file
    retval = nf90_create(filename, nf90_clobber, ncid)
    if (retval /= nf90_noerr) then
      print *, 'Failed to create file: ', trim(filename)
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Define the dimensions
    retval = nf90_def_dim(ncid, xname, xsize, xid)
    if (retval /= nf90_noerr) then
      print *, 'Failed to define dimension ', xname
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    retval = nf90_def_dim(ncid, yname, ysize, yid)
    if (retval /= nf90_noerr) then
      print *, 'Failed to define dimension ', yname
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Define the dimension variables
    retval = nf90_def_var(ncid, xname, nf90_double, [xid], xvarid)
    if (retval /= nf90_noerr) then
      print *, 'Failed to define variable: ', xname
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    retval = nf90_def_var(ncid, yname, nf90_double, [yid], yvarid)
    if (retval /= nf90_noerr) then
      print *, 'Failed to define variable: ', yname
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Add attributes to the dimension variables
    retval = nf90_put_att(ncid, xvarid, 'standard_name', 'longitude')
    retval = nf90_put_att(ncid, xvarid, 'long_name', 'longitude')
    retval = nf90_put_att(ncid, xvarid, 'units', 'degrees_east')
    retval = nf90_put_att(ncid, xvarid, 'axis', 'X')

    retval = nf90_put_att(ncid, yvarid, 'standard_name', 'latitude')
    retval = nf90_put_att(ncid, yvarid, 'long_name', 'latitude')
    retval = nf90_put_att(ncid, yvarid, 'units', 'degrees_north')
    retval = nf90_put_att(ncid, yvarid, 'axis', 'Y')

    ! Define the variable
    retval = nf90_def_var(ncid, varname, nf90_double, [xid, yid], varid)
    if (retval /= nf90_noerr) then
      print *, 'Failed to define variable: ', varname
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Add attributes to the variable
    retval = nf90_put_att(ncid, varid, 'units', varunit)
    retval = nf90_put_att(ncid, varid, '_FillValue', fill_value)
    retval = nf90_put_att(ncid, varid, 'missing_value', fill_value)
    retval = nf90_put_att(ncid, varid, 'cell_methods', 'time: mean')
    retval = nf90_put_att(ncid, varid, 'original_name', varname)

    ! End define mode
    retval = nf90_enddef(ncid)
    if (retval /= nf90_noerr) then
      print *, 'Failed to end define mode.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Write the data for the dimensions
    retval = nf90_put_var(ncid, xvarid, xvalues)
    if (retval /= nf90_noerr) then
      print *, 'Failed to write data: ', xname
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    retval = nf90_put_var(ncid, yvarid, yvalues)
    if (retval /= nf90_noerr) then
      print *, 'Failed to write data: ', yname
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Write the data for the variable
    retval = nf90_put_var(ncid, varid, var)
    if (retval /= nf90_noerr) then
      print *, 'Failed to write data for variable: ', varname
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Close the netCDF file
    retval = nf90_close(ncid)
    if (retval /= nf90_noerr) then
      print *, 'Error occurred while closing the file: ', trim(filename)
      stop
    endif

  end subroutine write_netcdf_2dwot


  subroutine write_netcdf_2dwt(filename, varname, varunit, var, xvalues, yvalues, tvalues)
    implicit none
    character(len=*), intent(in) :: filename, varname, varunit
    double precision, dimension(:,:,:), intent(in) :: var
    double precision, dimension(:), intent(in) :: xvalues, yvalues
    integer(kind=4),  dimension(:), intent(in) :: tvalues
    double precision, dimension(:), allocatable :: tvalues_real
    integer(kind=4) :: ncid, varid, retval
    integer(kind=4) :: xid, xsize
    integer(kind=4) :: yid, ysize
    integer(kind=4) :: tid, tsize
    integer(kind=4) :: xvarid, yvarid, tvarid
    integer(kind=4) :: time_index
    logical :: file_exists
    character(len=*), parameter :: xname = 'lon'
    character(len=*), parameter :: yname = 'lat'
    character(len=*), parameter :: tname = 'time'
    double precision, parameter :: fill_value = 999.0
    integer(kind=4), dimension(1) :: start_time, count_time
    integer(kind=4), dimension(3) :: start_var, count_var

    ! Deallocate arrays if they are already allocated
    if (allocated(tvalues_real)) deallocate(tvalues_real)

    ! Determine the size of each dimension from the input arrays
    xsize = size(xvalues)
    ysize = size(yvalues)
    tsize = size(tvalues)
    allocate(tvalues_real(tsize))
    tvalues_real = real(tvalues, kind=8)

    ! Check if file already exists
    inquire(file=filename, exist=file_exists)

    if (file_exists) then
      ! Open existing file in write mode
      retval = nf90_open(filename, nf90_write, ncid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to open existing file: ', trim(filename)
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      ! Get the time dimension
      retval = nf90_inq_dimid(ncid, tname, tid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to get time dimension.'
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      ! Get the current size of the time dimension
      retval = nf90_inquire_dimension(ncid, tid, len=tsize)
      if (retval /= nf90_noerr) then
        print *, 'Failed to get time dimension size.'
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      ! Get the time variable ID
      retval = nf90_inq_varid(ncid, tname, tvarid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to get time variable ID.'
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      ! Get the variable ID
      retval = nf90_inq_varid(ncid, varname, varid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to get variable ID.'
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      ! Set the new time index for the new data
      time_index = tsize + 1

    else
      ! Create the netCDF file
      retval = nf90_create(filename, nf90_clobber, ncid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to create file: ', trim(filename)
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      ! Define the dimensions
      retval = nf90_def_dim(ncid, xname, xsize, xid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to define dimension ', xname
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      retval = nf90_def_dim(ncid, yname, ysize, yid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to define dimension ', yname
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      retval = nf90_def_dim(ncid, tname, nf90_unlimited, tid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to define dimension ', tname
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      ! Get the current size of the time dimension
      retval = nf90_inquire_dimension(ncid, tid, len=tsize)
      if (retval /= nf90_noerr) then
        print *, 'Failed to get time dimension size.'
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      ! Define the dimension variables
      retval = nf90_def_var(ncid, xname, nf90_double, [xid], xvarid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to define variable: ', xname
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      retval = nf90_def_var(ncid, yname, nf90_double, [yid], yvarid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to define variable: ', yname
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      retval = nf90_def_var(ncid, tname, nf90_double, [tid], tvarid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to define variable: ', tname
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      ! Add attributes to the dimension variables
      retval = nf90_put_att(ncid, xvarid, 'standard_name', 'longitude')
      retval = nf90_put_att(ncid, xvarid, 'long_name', 'longitude')
      retval = nf90_put_att(ncid, xvarid, 'units', 'degrees_east')
      retval = nf90_put_att(ncid, xvarid, 'axis', 'X')

      retval = nf90_put_att(ncid, yvarid, 'standard_name', 'latitude')
      retval = nf90_put_att(ncid, yvarid, 'long_name', 'latitude')
      retval = nf90_put_att(ncid, yvarid, 'units', 'degrees_north')
      retval = nf90_put_att(ncid, yvarid, 'axis', 'Y')

      retval = nf90_put_att(ncid, tvarid, 'standard_name', 'time')
      retval = nf90_put_att(ncid, tvarid, 'long_name', 'time')
      retval = nf90_put_att(ncid, tvarid, 'units', 'year')

      ! Define the variable
      retval = nf90_def_var(ncid, varname, nf90_double, [xid, yid, tid], varid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to define variable: ', varname
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      ! Add attributes to the variable
      retval = nf90_put_att(ncid, varid, 'units', varunit)
      retval = nf90_put_att(ncid, varid, '_FillValue', fill_value)
      retval = nf90_put_att(ncid, varid, 'missing_value', fill_value)
      retval = nf90_put_att(ncid, varid, 'cell_methods', 'time: mean')
      retval = nf90_put_att(ncid, varid, 'original_name', varname)

      ! End define mode
      retval = nf90_enddef(ncid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to end define mode.'
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      ! Write the data for the dimensions
      retval = nf90_put_var(ncid, xvarid, xvalues)
      if (retval /= nf90_noerr) then
        print *, 'Failed to write data: ', xname
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      retval = nf90_put_var(ncid, yvarid, yvalues)
      if (retval /= nf90_noerr) then
        print *, 'Failed to write data: ', yname
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      ! Initialize time index for the new file
      time_index = 1

      ! Write the initial time value
      retval = nf90_put_var(ncid, tvarid, tvalues_real)
      if (retval /= nf90_noerr) then
        print *, 'Failed to write time value.'
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

    endif

    ! Write the new time value for the existing file
    if (file_exists) then
      ! Extend the time dimension and write the new time value
      start_time = (/ time_index /)
      count_time = (/ 1 /)
      retval = nf90_put_var(ncid, tvarid, tvalues_real, start=start_time)
      if (retval /= nf90_noerr) then
        print *, 'Failed to write new time value.'
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif
    endif

    ! Write the data for the variable
    start_var = (/ 1, 1, time_index /)
    count_var = (/ xsize, ysize, 1 /)
    retval = nf90_put_var(ncid, varid, var, start=start_var, count=count_var)
    if (retval /= nf90_noerr) then
      print *, 'Failed to write data for variable: ', varname
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Close the netCDF file
    retval = nf90_close(ncid)
    if (retval /= nf90_noerr) then
      print *, 'Error occurred while closing the file: ', trim(filename)
      stop
    endif

  end subroutine write_netcdf_2dwt


  subroutine write_netcdf_3dwot(filename, varname, varunit, var, xvalues, yvalues, zvalues)
    implicit none

    character(len=*), intent(in) :: filename, varname, varunit
    double precision, dimension(:,:,:), intent(in) :: var
    double precision, dimension(:), intent(in) :: xvalues, yvalues, zvalues

    integer(kind=4) :: ncid, varid, retval
    integer(kind=4) :: xid, xsize
    integer(kind=4) :: yid, ysize
    integer(kind=4) :: zid, zsize
    integer(kind=4) :: xvarid, yvarid, zvarid

    character(len=*), parameter :: xname = 'lon'
    character(len=*), parameter :: yname = 'lat'
    character(len=*), parameter :: zname = 'depth'

    double precision, parameter :: fill_value = 999.0

    ! Determine the size of each dimension from the input arrays
    xsize = size(xvalues)
    ysize = size(yvalues)
    zsize = size(zvalues)

    ! Create the netCDF file
    retval = nf90_create(filename, nf90_clobber, ncid)
    if (retval /= nf90_noerr) then
      print *, 'Failed to create file: ', trim(filename)
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Define the dimensions
    retval = nf90_def_dim(ncid, xname, xsize, xid)
    if (retval /= nf90_noerr) then
      print *, 'Failed to define dimension ', xname
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    retval = nf90_def_dim(ncid, yname, ysize, yid)
    if (retval /= nf90_noerr) then
      print *, 'Failed to define dimension ', yname
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    retval = nf90_def_dim(ncid, zname, zsize, zid)
    if (retval /= nf90_noerr) then
      print *, 'Failed to define dimension ', zname
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Define the dimension variables
    retval = nf90_def_var(ncid, xname, nf90_double, [xid], xvarid)
    if (retval /= nf90_noerr) then
      print *, 'Failed to define variable: ', xname
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    retval = nf90_def_var(ncid, yname, nf90_double, [yid], yvarid)
    if (retval /= nf90_noerr) then
      print *, 'Failed to define variable: ', yname
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    retval = nf90_def_var(ncid, zname, nf90_double, [zid], zvarid)
    if (retval /= nf90_noerr) then
      print *, 'Failed to define variable: ', zname
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Add attributes to the dimension variables
    retval = nf90_put_att(ncid, xvarid, 'standard_name', 'longitude')
    retval = nf90_put_att(ncid, xvarid, 'long_name', 'longitude')
    retval = nf90_put_att(ncid, xvarid, 'units', 'degrees_east')
    retval = nf90_put_att(ncid, xvarid, 'axis', 'X')

    retval = nf90_put_att(ncid, yvarid, 'standard_name', 'latitude')
    retval = nf90_put_att(ncid, yvarid, 'long_name', 'latitude')
    retval = nf90_put_att(ncid, yvarid, 'units', 'degrees_north')
    retval = nf90_put_att(ncid, yvarid, 'axis', 'Y')

    retval = nf90_put_att(ncid, zvarid, 'standard_name', 'depth')
    retval = nf90_put_att(ncid, zvarid, 'long_name', 'depth')
    retval = nf90_put_att(ncid, zvarid, 'units', 'm')
    retval = nf90_put_att(ncid, zvarid, 'axis', 'Z')

    ! Define the variable
    retval = nf90_def_var(ncid, varname, nf90_double, [xid, yid, zid], varid)
    if (retval /= nf90_noerr) then
      print *, 'Failed to define variable: ', varname
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Add attributes to the variable
    retval = nf90_put_att(ncid, varid, 'units', varunit)
    retval = nf90_put_att(ncid, varid, '_FillValue', fill_value)
    retval = nf90_put_att(ncid, varid, 'missing_value', fill_value)
    retval = nf90_put_att(ncid, varid, 'cell_methods', 'time: mean')
    retval = nf90_put_att(ncid, varid, 'original_name', varname)

    ! End define mode
    retval = nf90_enddef(ncid)
    if (retval /= nf90_noerr) then
      print *, 'Failed to end define mode.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Write the data for the dimensions
    retval = nf90_put_var(ncid, xvarid, xvalues)
    if (retval /= nf90_noerr) then
      print *, 'Failed to write data: ', xname
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    retval = nf90_put_var(ncid, yvarid, yvalues)
    if (retval /= nf90_noerr) then
      print *, 'Failed to write data: ', yname
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    retval = nf90_put_var(ncid, zvarid, zvalues)
    if (retval /= nf90_noerr) then
      print *, 'Failed to write data: ', zname
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Write the data for the variable
    retval = nf90_put_var(ncid, varid, var)
    if (retval /= nf90_noerr) then
      print *, 'Failed to write data for variable: ', varname
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Close the netCDF file
    retval = nf90_close(ncid)
    if (retval /= nf90_noerr) then
      print *, 'Error occurred while closing the file: ', trim(filename)
      stop
    endif

  end subroutine write_netcdf_3dwot


  subroutine write_netcdf_3dwt(filename, varname, varunit, var, xvalues, yvalues, zvalues, tvalues)
    implicit none
    character(len=*), intent(in) :: filename, varname, varunit
    double precision, dimension(:,:,:,:), intent(in) :: var
    double precision, dimension(:), intent(in) :: xvalues, yvalues, zvalues
    integer(kind=4),  dimension(:), intent(in) :: tvalues
    double precision, dimension(:), allocatable :: tvalues_real
    integer(kind=4) :: ncid, varid, retval
    integer(kind=4) :: xid, xsize
    integer(kind=4) :: yid, ysize
    integer(kind=4) :: zid, zsize
    integer(kind=4) :: tid, tsize
    integer(kind=4) :: xvarid, yvarid, zvarid, tvarid
    integer(kind=4) :: time_index
    logical :: file_exists
    character(len=*), parameter :: xname = 'lon'
    character(len=*), parameter :: yname = 'lat'
    character(len=*), parameter :: zname = 'depth'
    character(len=*), parameter :: tname = 'time'
    double precision, parameter :: fill_value = 999.0
    integer(kind=4), dimension(1) :: start_time, count_time
    integer(kind=4), dimension(4) :: start_var, count_var

    ! Deallocate arrays if they are already allocated
    if (allocated(tvalues_real)) deallocate(tvalues_real)

    ! Determine the size of each dimension from the input arrays
    xsize = size(xvalues)
    ysize = size(yvalues)
    zsize = size(zvalues)
    tsize = size(tvalues)
    allocate(tvalues_real(tsize))
    tvalues_real = real(tvalues, kind=8)

    ! Check if file already exists
    inquire(file=filename, exist=file_exists)

    if (file_exists) then
      ! Open existing file in write mode
      retval = nf90_open(filename, nf90_write, ncid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to open existing file: ', trim(filename)
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      ! Get the time dimension
      retval = nf90_inq_dimid(ncid, tname, tid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to get time dimension.'
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      ! Get the current size of the time dimension
      retval = nf90_inquire_dimension(ncid, tid, len=tsize)
      if (retval /= nf90_noerr) then
        print *, 'Failed to get time dimension size.'
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      ! Get the time variable ID
      retval = nf90_inq_varid(ncid, tname, tvarid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to get time variable ID.'
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      ! Get the variable ID
      retval = nf90_inq_varid(ncid, varname, varid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to get variable ID.'
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      ! Set the new time index for the new data
      time_index = tsize + 1

    else
      ! Create the netCDF file
      retval = nf90_create(filename, nf90_clobber, ncid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to create file: ', trim(filename)
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      ! Define the dimensions
      retval = nf90_def_dim(ncid, xname, xsize, xid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to define dimension ', xname
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      retval = nf90_def_dim(ncid, yname, ysize, yid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to define dimension ', yname
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      retval = nf90_def_dim(ncid, zname, zsize, zid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to define dimension ', zname
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      retval = nf90_def_dim(ncid, tname, nf90_unlimited, tid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to define dimension ', tname
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      ! Get the current size of the time dimension
      retval = nf90_inquire_dimension(ncid, tid, len=tsize)
      if (retval /= nf90_noerr) then
        print *, 'Failed to get time dimension size.'
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      ! Define the dimension variables
      retval = nf90_def_var(ncid, xname, nf90_double, [xid], xvarid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to define variable: ', xname
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      retval = nf90_def_var(ncid, yname, nf90_double, [yid], yvarid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to define variable: ', yname
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      retval = nf90_def_var(ncid, zname, nf90_double, [zid], zvarid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to define variable: ', zname
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      retval = nf90_def_var(ncid, tname, nf90_double, [tid], tvarid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to define variable: ', tname
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      ! Add attributes to the dimension variables
      retval = nf90_put_att(ncid, xvarid, 'standard_name', 'longitude')
      retval = nf90_put_att(ncid, xvarid, 'long_name', 'longitude')
      retval = nf90_put_att(ncid, xvarid, 'units', 'degrees_east')
      retval = nf90_put_att(ncid, xvarid, 'axis', 'X')

      retval = nf90_put_att(ncid, yvarid, 'standard_name', 'latitude')
      retval = nf90_put_att(ncid, yvarid, 'long_name', 'latitude')
      retval = nf90_put_att(ncid, yvarid, 'units', 'degrees_north')
      retval = nf90_put_att(ncid, yvarid, 'axis', 'Y')

      retval = nf90_put_att(ncid, zvarid, 'standard_name', 'depth')
      retval = nf90_put_att(ncid, zvarid, 'long_name', 'depth')
      retval = nf90_put_att(ncid, zvarid, 'units', 'm')
      retval = nf90_put_att(ncid, zvarid, 'axis', 'Z')

      retval = nf90_put_att(ncid, tvarid, 'standard_name', 'time')
      retval = nf90_put_att(ncid, tvarid, 'long_name', 'time')
      retval = nf90_put_att(ncid, tvarid, 'units', 'year')

      ! Define the variable
      retval = nf90_def_var(ncid, varname, nf90_double, [xid, yid, zid, tid], varid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to define variable: ', varname
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      ! Add attributes to the variable
      retval = nf90_put_att(ncid, varid, 'units', varunit)
      retval = nf90_put_att(ncid, varid, '_FillValue', fill_value)
      retval = nf90_put_att(ncid, varid, 'missing_value', fill_value)
      retval = nf90_put_att(ncid, varid, 'cell_methods', 'time: mean')
      retval = nf90_put_att(ncid, varid, 'original_name', varname)

      ! End define mode
      retval = nf90_enddef(ncid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to end define mode.'
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      ! Write the data for the dimensions
      retval = nf90_put_var(ncid, xvarid, xvalues)
      if (retval /= nf90_noerr) then
        print *, 'Failed to write data: ', xname
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      retval = nf90_put_var(ncid, yvarid, yvalues)
      if (retval /= nf90_noerr) then
        print *, 'Failed to write data: ', yname
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      retval = nf90_put_var(ncid, zvarid, zvalues)
      if (retval /= nf90_noerr) then
        print *, 'Failed to write data: ', zname
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      ! Initialize time index for the new file
      time_index = 1

      ! Write the initial time value
      retval = nf90_put_var(ncid, tvarid, tvalues_real)
      if (retval /= nf90_noerr) then
        print *, 'Failed to write time value.'
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

    endif

    ! Write the new time value for the existing file
    if (file_exists) then
      ! Extend the time dimension and write the new time value
      start_time = (/ time_index /)
      count_time = (/ 1 /)
      retval = nf90_put_var(ncid, tvarid, tvalues_real, start=start_time)
      if (retval /= nf90_noerr) then
        print *, 'Failed to write new time value.'
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif
    endif

    ! Write the data for the variable
    start_var = (/ 1, 1, 1, time_index /)
    count_var = (/ xsize, ysize, zsize, 1 /)
    retval = nf90_put_var(ncid, varid, var, start=start_var, count=count_var)
    if (retval /= nf90_noerr) then
      print *, 'Failed to write data for variable: ', varname
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Close the netCDF file
    retval = nf90_close(ncid)
    if (retval /= nf90_noerr) then
      print *, 'Error occurred while closing the file: ', trim(filename)
      stop
    endif

  end subroutine write_netcdf_3dwt


!!subroutine write_restart(filename, varname, varunit, var, xvalues, yvalues, zvalues, tvalues)
!!  use netcdf
!!  implicit none

!!  character(len=*), intent(in) :: filename, varname, varunit
!!  double precision, dimension(:,:,:,:), intent(in) :: var
!!  double precision, dimension(:), intent(in) :: xvalues, yvalues, zvalues
!!  integer(kind=4),  dimension(:), intent(in) :: tvalues

!!  double precision, dimension(:), allocatable :: tvalues_real
!!  integer(kind=4) :: ncid, varid, retval
!!  integer(kind=4) :: xid, xsize
!!  integer(kind=4) :: yid, ysize
!!  integer(kind=4) :: zid, zsize
!!  integer(kind=4) :: tid, tsize
!!  integer(kind=4) :: xvarid, yvarid, zvarid, tvarid
!!  logical :: file_exists
!!  character(len=*), parameter :: xname = 'lon'
!!  character(len=*), parameter :: yname = 'lat'
!!  character(len=*), parameter :: zname = 'depth'
!!  character(len=*), parameter :: tname = 'time'
!!  double precision, parameter :: fill_value = 999.0
!!  integer(kind=4), dimension(4) :: start_var, count_var

!!  ! Deallocate arrays if they are already allocated
!!  if (allocated(tvalues_real)) deallocate(tvalues_real)

!!  ! Determine the size of each dimension from the input arrays
!!  xsize = size(xvalues)
!!  ysize = size(yvalues)
!!  zsize = size(zvalues)
!!  tsize = size(tvalues)
!!  allocate(tvalues_real(tsize))
!!  tvalues_real = real(tvalues, kind=8)

!!  ! Check if file already exists
!!  inquire(file=filename, exist=file_exists)

!!  if (file_exists) then
!!    ! Open existing file in write mode
!!    retval = nf90_open(filename, nf90_write, ncid)
!!    if (retval /= nf90_noerr) then
!!      print *, 'Failed to open existing file: ', trim(filename)
!!      print *, 'Error: ', nf90_strerror(retval)
!!      stop
!!    endif

!!    ! Get the time dimension
!!    retval = nf90_inq_dimid(ncid, tname, tid)
!!    if (retval /= nf90_noerr) then
!!      print *, 'Failed to get time dimension.'
!!      print *, 'Error: ', nf90_strerror(retval)
!!      stop
!!    endif

!!    ! Get the current size of the time dimension
!!    retval = nf90_inquire_dimension(ncid, tid, len=tsize)
!!    if (retval /= nf90_noerr) then
!!      print *, 'Failed to get time dimension size.'
!!      print *, 'Error: ', nf90_strerror(retval)
!!      stop
!!    endif

!!    ! Check if the variable already exists
!!    retval = nf90_inq_varid(ncid, varname, varid)
!!    if (retval /= nf90_noerr) then
!!      ! If the variable does not exist, define it
!!      retval = nf90_redef(ncid)
!!      if (retval /= nf90_noerr) then
!!        print *, 'Failed to enter define mode.'
!!        print *, 'Error: ', nf90_strerror(retval)
!!        stop
!!      endif

!!      ! Define the data variable
!!      retval = nf90_def_var(ncid, varname, nf90_double, [xid, yid, zid, tid], varid)
!!      if (retval /= nf90_noerr) then
!!        print *, 'Failed to define variable: ', varname
!!        print *, 'Error: ', nf90_strerror(retval)
!!        stop
!!      endif

!!      ! Add units attribute to the data variable
!!      retval = nf90_put_att(ncid, varid, 'units', varunit)
!!      if (retval /= nf90_noerr) then
!!        print *, 'Failed to add units attribute to variable ', varname
!!        print *, 'Error: ', nf90_strerror(retval)
!!        stop
!!      endif

!!      ! End define mode
!!      retval = nf90_enddef(ncid)
!!      if (retval /= nf90_noerr) then
!!        print *, 'Failed to end define mode.'
!!        print *, 'Error: ', nf90_strerror(retval)
!!        stop
!!      endif
!!    endif

!!  else
!!    ! Create the netCDF file
!!    retval = nf90_create(filename, nf90_clobber, ncid)
!!    if (retval /= nf90_noerr) then
!!      print *, 'Failed to create file: ', trim(filename)
!!      print *, 'Error: ', nf90_strerror(retval)
!!      stop
!!    endif

!!    ! Define the dimensions
!!    retval = nf90_def_dim(ncid, xname, xsize, xid)
!!    if (retval /= nf90_noerr) then
!!      print *, 'Failed to define dimension ', xname
!!      print *, 'Error: ', nf90_strerror(retval)
!!      stop
!!    endif

!!    retval = nf90_def_dim(ncid, yname, ysize, yid)
!!    if (retval /= nf90_noerr) then
!!      print *, 'Failed to define dimension ', yname
!!      print *, 'Error: ', nf90_strerror(retval)
!!      stop
!!    endif

!!    retval = nf90_def_dim(ncid, zname, zsize, zid)
!!    if (retval /= nf90_noerr) then
!!      print *, 'Failed to define dimension ', zname
!!      print *, 'Error: ', nf90_strerror(retval)
!!      stop
!!    endif

!!    retval = nf90_def_dim(ncid, tname, nf90_unlimited, tid)
!!    if (retval /= nf90_noerr) then
!!      print *, 'Failed to define dimension ', tname
!!      print *, 'Error: ', nf90_strerror(retval)
!!      stop
!!    endif

!!    ! Define the dimension variables
!!    retval = nf90_def_var(ncid, xname, nf90_double, [xid], xvarid)
!!    if (retval /= nf90_noerr) then
!!      print *, 'Failed to define variable: ', xname
!!      print *, 'Error: ', nf90_strerror(retval)
!!      stop
!!    endif

!!    retval = nf90_def_var(ncid, yname, nf90_double, [yid], yvarid)
!!    if (retval /= nf90_noerr) then
!!      print *, 'Failed to define variable: ', yname
!!      print *, 'Error: ', nf90_strerror(retval)
!!      stop
!!    endif

!!    retval = nf90_def_var(ncid, zname, nf90_double, [zid], zvarid)
!!    if (retval /= nf90_noerr) then
!!      print *, 'Failed to define variable: ', zname
!!      print *, 'Error: ', nf90_strerror(retval)
!!      stop
!!    endif

!!    retval = nf90_def_var(ncid, tname, nf90_double, [tid], tvarid)
!!    if (retval /= nf90_noerr) then
!!      print *, 'Failed to define variable: ', tname
!!      print *, 'Error: ', nf90_strerror(retval)
!!      stop
!!    endif

!!    ! Define the data variable
!!    retval = nf90_def_var(ncid, varname, nf90_double, [xid, yid, zid, tid], varid)
!!    if (retval /= nf90_noerr) then
!!      print *, 'Failed to define variable: ', varname
!!      print *, 'Error: ', nf90_strerror(retval)
!!      stop
!!    endif

!!    ! Add units attribute to the data variable
!!    retval = nf90_put_att(ncid, varid, 'units', varunit)
!!    if (retval /= nf90_noerr) then
!!      print *, 'Failed to add units attribute to variable ', varname
!!      print *, 'Error: ', nf90_strerror(retval)
!!      stop
!!    endif

!!    ! End define mode
!!    retval = nf90_enddef(ncid)
!!    if (retval /= nf90_noerr) then
!!      print *, 'Failed to end define mode.'
!!      print *, 'Error: ', nf90_strerror(retval)
!!      stop
!!    endif

!!    ! Write the dimension variables
!!    retval = nf90_put_var(ncid, xvarid, xvalues)
!!    if (retval /= nf90_noerr) then
!!      print *, 'Failed to write longitude values.'
!!      print *, 'Error: ', nf90_strerror(retval)
!!      stop
!!    endif

!!    retval = nf90_put_var(ncid, yvarid, yvalues)
!!    if (retval /= nf90_noerr) then
!!      print *, 'Failed to write latitude values.'
!!      print *, 'Error: ', nf90_strerror(retval)
!!      stop
!!    endif

!!    retval = nf90_put_var(ncid, zvarid, zvalues)
!!    if (retval /= nf90_noerr) then
!!      print *, 'Failed to write depth values.'
!!      print *, 'Error: ', nf90_strerror(retval)
!!      stop
!!    endif

!!    retval = nf90_put_var(ncid, tvarid, tvalues_real)
!!    if (retval /= nf90_noerr) then
!!      print *, 'Failed to write time values.'
!!      print *, 'Error: ', nf90_strerror(retval)
!!      stop
!!    endif
!!  endif

!!  ! Write data to the variable
!!  start_var = [1, 1, 1, 1]
!!  count_var = [xsize, ysize, zsize, tsize]

!!  retval = nf90_put_var(ncid, varid, var, start=start_var, count=count_var)
!!  if (retval /= nf90_noerr) then
!!    print *, 'Failed to write data to variable.'
!!    print *, 'Error: ', nf90_strerror(retval)
!!    stop
!!  endif

!!  ! Close the NetCDF file
!!  retval = nf90_close(ncid)
!!  if (retval /= nf90_noerr) then
!!    print *, 'Error occurred while closing the file.'
!!    print *, 'Error: ', nf90_strerror(retval)
!!    stop
!!  endif

!!  print *, "Successfully wrote to file: ", trim(filename)

!!end subroutine write_restart

  subroutine combine_netcdf_files(infile1, infile2, outfile)
    implicit none
    character(len=*), intent(in) :: infile1, infile2
    character(len=*), intent(in) :: outfile
    integer(kind=4) :: ncid1, ncid2, ncid_out
    integer(kind=4) :: varid1, varid2, varid_out
    integer(kind=4) :: dimid1, dimid2, dimid_out
    integer(kind=4) :: xid, xsize
    integer(kind=4) :: yid, ysize
    integer(kind=4) :: zid, zsize
    integer(kind=4) :: tid, tsize
    integer(kind=4) :: xvarid, yvarid, zvarid, tvarid
    integer, dimension(4) :: dimids

    character(len=32) :: varname1, varname2
    double precision, dimension(:,:,:,:), allocatable :: data1, data2

    character(len=*), parameter :: xname = 'lon'
    character(len=*), parameter :: yname = 'lat'
    character(len=*), parameter :: zname = 'depth'
    character(len=*), parameter :: tname = 'time'

    integer(kind=4) :: retval

    ! Deallocate arrays if they are already allocated
    if (allocated(data1)) deallocate(data1)
    if (allocated(data2)) deallocate(data2)

    ! Open the first NetCDF file
    retval = nf90_open(infile1, nf90_nowrite, ncid1)
    if (retval /= nf90_noerr) then
      print *, 'Failed to open existing file: ', trim(infile1)
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Open the second NetCDF file
    retval = nf90_open(infile2, nf90_nowrite, ncid2)
    if (retval /= nf90_noerr) then
      print *, 'Failed to open existing file: ', trim(infile2)
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Create the new NetCDF file
    retval = nf90_create(outfile, nf90_clobber, ncid_out)
    if (retval /= nf90_noerr) then
      print *, 'Failed to create file: ', trim(outfile)
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Get the time dimension
    retval = nf90_inq_dimid(ncid1, tname, tid)
    if (retval /= nf90_noerr) then
      print *, 'Failed to get time dimension.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Get the current size of the time dimension
    retval = nf90_inquire_dimension(ncid1, tid, len=tsize)
    if (retval /= nf90_noerr) then
      print *, 'Failed to get time dimension size.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Get the time variable ID
    retval = nf90_inq_varid(ncid1, tname, tvarid)
    if (retval /= nf90_noerr) then
      print *, 'Failed to get time variable ID.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Get the variable ID
    retval = nf90_inq_varid(ncid1, varname1, varid1)
    if (retval /= nf90_noerr) then
      print *, 'Failed to get variable ID.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Get the variable ID
    retval = nf90_inq_varid(ncid2, varname2, varid2)
    if (retval /= nf90_noerr) then
      print *, 'Failed to get variable ID.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Define dimensions and variables for the new file
    ! Get lon dimension
    retval = nf90_inq_dimid(ncid1, xname, xid)
    if (retval /= nf90_noerr) then
      print *, 'Failed to get ID of longitude dimension.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif
    retval = nf90_inquire_dimension(ncid1, xid, len=xsize)
    if (retval /= nf90_noerr) then
      print *, 'Failed to get length of longitude dimension.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Get lat dimension
    retval = nf90_inq_dimid(ncid1, yname, yid)
    if (retval /= nf90_noerr) then
      print *, 'Failed to get ID of latitude dimension.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif
    retval = nf90_inquire_dimension(ncid1, yid, len=ysize)
    if (retval /= nf90_noerr) then
      print *, 'Failed to get length of latitude dimension.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Get depth dimension
    retval = nf90_inq_dimid(ncid1, zname, zid)
    if (retval /= nf90_noerr) then
      print *, 'Failed to get ID of depth dimension.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif
    retval = nf90_inquire_dimension(ncid1, zid, len=zsize)
    if (retval /= nf90_noerr) then
      print *, 'Failed to get length of depth dimension.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Get time dimension
    retval = nf90_inq_dimid(ncid1, tname, tid)
    if (retval /= nf90_noerr) then
      print *, 'Failed to get ID of time dimension.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif
    retval = nf90_inquire_dimension(ncid1, tid, len=tsize)
    if (retval /= nf90_noerr) then
      print *, 'Failed to get length of time dimension.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Define dimensions in the new file
    retval = nf90_def_dim(ncid_out, xname, xsize, xid)
    if (retval /= nf90_noerr) then
      print *, 'Failed to define longitude dimension.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    retval = nf90_def_dim(ncid_out, yname, ysize, yid)
    if (retval /= nf90_noerr) then
      print *, 'Failed to define latitude dimension.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    retval = nf90_def_dim(ncid_out, zname, zsize, zid)
    if (retval /= nf90_noerr) then
      print *, 'Failed to define depth dimension.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    retval = nf90_def_dim(ncid_out, tname, tsize, tid)
    if (retval /= nf90_noerr) then
      print *, 'Failed to define time dimension.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Define variables in the new file
    dimids = (/tid, zid, yid, xid/)

    ! Get variable IDs from the input files
    retval = nf90_inq_varid(ncid1, varname1, varid1)
    if (retval /= nf90_noerr) then
      print *, 'Failed to retrieve variable ID: ', varname1
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    retval = nf90_inq_varid(ncid2, varname2, varid2)
    if (retval /= nf90_noerr) then
      print *, 'Failed to retrieve variable ID: ', varname2
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Define the variables in the new file
    retval = nf90_def_var(ncid_out, varname1, nf90_double, dimids, varid_out)
    if (retval /= nf90_noerr) then
      print *, 'Failed to define variable: ', varname1
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    retval = nf90_def_var(ncid_out, varname2, nf90_double, dimids, varid_out)
    if (retval /= nf90_noerr) then
      print *, 'Failed to define variable: ', varname2
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! End define mode
    retval = nf90_enddef(ncid_out)
    if (retval /= nf90_noerr) then
      print *, 'Failed to end define mode.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Allocate and read data from the first file
    allocate(data1(tsize, zsize, ysize, xsize))
    retval = nf90_get_var(ncid1, varid1, data1)
    if (retval /= nf90_noerr) then
      print *, 'Failed to read variable: ', varname1
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Allocate and read data from the second file
    allocate(data2(tsize, zsize, ysize, xsize))
    retval = nf90_get_var(ncid2, varid2, data2)
    if (retval /= nf90_noerr) then
      print *, 'Failed to read variable: ', varname2
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Write data to the new file
    retval = nf90_put_var(ncid_out, varid1, data1)
    if (retval /= nf90_noerr) then
      print *, 'Failed to write variable: ', varname1
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    retval = nf90_put_var(ncid_out, varid2, data2)
    if (retval /= nf90_noerr) then
      print *, 'Failed to write variable: ', varname2
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Close all files
    retval = nf90_close(ncid1)
    if (retval /= nf90_noerr) then
      print *, 'Error occurred while closing the file: ', trim(infile1)
      stop
    endif

    retval = nf90_close(ncid2)
    if (retval /= nf90_noerr) then
      print *, 'Error occurred while closing the file: ', trim(infile2)
      stop
    endif

    retval = nf90_close(ncid_out)
    if (retval /= nf90_noerr) then
      print *, 'Error occurred while closing the file: ', trim(outfile)
      stop
    endif

    ! Deallocate arrays
    if (allocated(data1)) deallocate(data1)
    if (allocated(data2)) deallocate(data2)

  end subroutine combine_netcdf_files

!!  subroutine combine_multiple_netcdf_files(infiles, outfile)
!!    implicit none
!!    character(len=*), dimension(:), intent(in) :: infiles
!!    character(len=*), intent(in) :: outfile
!!    integer(kind=4) :: ncid_in, ncid_out
!!    integer(kind=4) :: varid_in, varid_out
!!    integer(kind=4) :: dimid_lon, dimid_lat, dimid_depth, dimid_time
!!    integer(kind=4) :: lon_size, lat_size, depth_size, time_size
!!    integer(kind=4) :: retval
!!    integer, dimension(4) :: dimids
!!    character(len=32) :: varname
!!    double precision, dimension(:,:,:,:), allocatable :: data
!!    character(len=*), parameter :: lon_name = 'lon', lat_name = 'lat', depth_name = 'depth', time_name = 'time'
!!    integer :: i, nfiles
!!
!!    nfiles = size(infiles)
!!
!!    ! Open the first input NetCDF file to get dimension sizes and variable names
!!    retval = nf90_open(trim(infiles(1)), nf90_nowrite, ncid_in)
!!    if (retval /= nf90_noerr) then
!!        print *, 'Failed to open existing file: ', trim(infiles(1))
!!        print *, 'Error: ', nf90_strerror(retval)
!!        stop
!!    endif
!!
!!    ! Get dimension IDs and sizes from the first input file
!!    call get_dimension_info(ncid_in, lon_name, dimid_lon, lon_size)
!!    call get_dimension_info(ncid_in, lat_name, dimid_lat, lat_size)
!!    call get_dimension_info(ncid_in, depth_name, dimid_depth, depth_size)
!!    call get_dimension_info(ncid_in, time_name, dimid_time, time_size)
!!
!!    ! Create the output NetCDF file
!!    retval = nf90_create(outfile, nf90_clobber, ncid_out)
!!    if (retval /= nf90_noerr) then
!!        print *, 'Failed to create file: ', trim(outfile)
!!        print *, 'Error: ', nf90_strerror(retval)
!!        stop
!!    endif
!!
!!    ! Define dimensions in the output file
!!    retval = nf90_def_dim(ncid_out, lon_name, lon_size, dimid_lon)
!!    if (retval /= nf90_noerr) call handle_error(retval)
!!    retval = nf90_def_dim(ncid_out, lat_name, lat_size, dimid_lat)
!!    if (retval /= nf90_noerr) call handle_error(retval)
!!    retval = nf90_def_dim(ncid_out, depth_name, depth_size, dimid_depth)
!!    if (retval /= nf90_noerr) call handle_error(retval)
!!    retval = nf90_def_dim(ncid_out, time_name, time_size, dimid_time)
!!    if (retval /= nf90_noerr) call handle_error(retval)
!!
!!    dimids = (/dimid_time, dimid_depth, dimid_lat, dimid_lon/)
!!
!!    ! Define variables in the output file
!!    do i = 1, nfiles
!!        retval = nf90_open(trim(infiles(i)), nf90_nowrite, ncid_in)
!!        if (retval /= nf90_noerr) then
!!            print *, 'Failed to open existing file: ', trim(infiles(i))
!!            print *, 'Error: ', nf90_strerror(retval)
!!            stop
!!        endif
!!
!!        call get_variable_name(ncid_in, varname)
!!
!!        retval = nf90_def_var(ncid_out, varname, nf90_double, dimids, varid_out)
!!        if (retval /= nf90_noerr) call handle_error(retval)
!!
!!        retval = nf90_enddef(ncid_out)
!!        if (retval /= nf90_noerr) call handle_error(retval)
!!
!!        allocate(data(time_size, depth_size, lat_size, lon_size))
!!        retval = nf90_inq_varid(ncid_in, varname, varid_in)
!!        if (retval /= nf90_noerr) call handle_error(retval)
!!        retval = nf90_get_var(ncid_in, varid_in, data)
!!        if (retval /= nf90_noerr) call handle_error(retval)
!!
!!        retval = nf90_put_var(ncid_out, varid_out, data)
!!        if (retval /= nf90_noerr) call handle_error(retval)
!!
!!        retval = nf90_close(ncid_in)
!!        if (retval /= nf90_noerr) call handle_error(retval)
!!
!!        if (allocated(data)) deallocate(data)
!!    end do
!!
!!    ! Close the output file
!!    retval = nf90_close(ncid_out)
!!    if (retval /= nf90_noerr) call handle_error(retval)
!!
!!contains
!!
!!    subroutine get_dimension_info(ncid, dimname, dimid, dimsize)
!!        implicit none
!!        integer, intent(in) :: ncid
!!        character(len=*), intent(in) :: dimname
!!        integer, intent(out) :: dimid
!!        integer, intent(out) :: dimsize
!!        integer :: retval
!!
!!        retval = nf90_inq_dimid(ncid, dimname, dimid)
!!        if (retval /= nf90_noerr) call handle_error(retval)
!!        retval = nf90_inquire_dimension(ncid, dimid, len=dimsize)
!!        if (retval /= nf90_noerr) call handle_error(retval)
!!    end subroutine get_dimension_info
!!
!!    subroutine get_variable_name(ncid, varname)
!!        implicit none
!!        integer, intent(in) :: ncid
!!        character(len=32), intent(out) :: varname
!!        integer :: nvars, retval
!!
!!        retval = nf90_inq_nvars(ncid, nvars)
!!        if (retval /= nf90_noerr) call handle_error(retval)
!!
!!        ! Assuming there is only one variable in each file for simplicity
!!        retval = nf90_inq_varname(ncid, 1, varname)
!!        if (retval /= nf90_noerr) call handle_error(retval)
!!    end subroutine get_variable_name
!!
!!    subroutine handle_error(retval)
!!        integer, intent(in) :: retval
!!        print *, 'Error:', trim(nf90_strerror(retval))
!!        stop 2
!!    end subroutine handle_error
!!
!!  end subroutine combine_multiple_netcdf_files

end module netcdf_utilities
