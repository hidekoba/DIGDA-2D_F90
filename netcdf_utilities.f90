module netcdf_utilities
    use netcdf
    implicit none
contains

  subroutine read_netcdf_2d(filename, varname, var, lon_values, lat_values)
    implicit none

    character(len=*), intent(in) :: filename, varname
    double precision, dimension(:,:), allocatable, intent(out) :: var

    integer(kind=4) :: ncid, varid, retval
    integer(kind=4) :: lonid, lonsize
    integer(kind=4) :: latid, latsize

    character(len=*), parameter :: lonname = 'lon'
    character(len=*), parameter :: latname = 'lat'

    double precision, allocatable, dimension(:) :: lon_values, lat_values
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
    retval = nf90_inq_dimid(ncid, lonname, lonid)
    if (retval == nf90_noerr) then
      retval = nf90_inquire_dimension(ncid, lonid, len=lonsize)
    else
      print *, 'Failed to get ID and length of longitude dimension.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif
   !print *, 'Longitude size = ', lonsize

    ! Get the latitude dimension ID and size
    retval = nf90_inq_dimid(ncid, latname, latid)
    if (retval == nf90_noerr) then
      retval = nf90_inquire_dimension(ncid, latid, len=latsize)
    else
      print *, 'Failed to get ID and length of latitude dimension.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif
   !print *, 'Latitude size = ', latsize

    ! Deallocate arrays if they are already allocated
    if (allocated(lon_values)) deallocate(lon_values)
    if (allocated(lat_values)) deallocate(lat_values)
    if (allocated(var)) deallocate(var)

    ! Allocate memory for the dimension values
    allocate(lon_values(lonsize))
    allocate(lat_values(latsize))

    ! Retrieve the longitude values
    retval = nf90_get_var(ncid, lonid, lon_values)
    if (retval /= nf90_noerr) then
      print *, 'Failed to read longitude values.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif
   !print *, 'Longitude values:'
   !do i = 1, lonsize
   !  print *, lon_values(i)
   !end do

    ! Retrieve the latitude values
    retval = nf90_get_var(ncid, latid, lat_values)
    if (retval /= nf90_noerr) then
      print *, 'Failed to read latitude values.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif
   !print *, 'Latitude values:'
   !do i = 1, latsize
   !  print *, lat_values(i)
   !end do

    ! Allocate memory for the 2D data variable based on retrieved dimensions
    allocate(var(lonsize, latsize))
    if (.not. allocated(var)) then
      print *, 'Allocation of ', varname, ' failed.'
      stop
    endif

    ! Retrieve the variable ID and read the data into the allocated array
    retval = nf90_inq_varid(ncid, varname, varid)
    if (retval /= nf90_noerr) then
      print *, 'Failed to retrieve variable ID for ', varname, '.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif
    retval = nf90_get_var(ncid, varid, var)
    if (retval /= nf90_noerr) then
      print *, 'Failed to read data for variable ', varname, '.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Close the netCDF file
    retval = nf90_close(ncid)
    if (retval /= nf90_noerr) then
      print *, 'Error occurred while closing the file.'
      stop
    endif

  end subroutine read_netcdf_2d


  subroutine read_netcdf_3d(filename, varname, var, lon_values, lat_values, depth_values)
    implicit none

    character(len=*), intent(in) :: filename, varname
    double precision, dimension(:,:,:), allocatable, intent(out) :: var
    double precision, dimension(:), allocatable, intent(out) :: lon_values, lat_values, depth_values

    integer(kind=4) :: ncid, varid, retval
    integer(kind=4) :: lonid, lonsize
    integer(kind=4) :: latid, latsize
    integer(kind=4) :: depthid, depthsize

    character(len=*), parameter :: lonname = 'lon'
    character(len=*), parameter :: latname = 'lat'
    character(len=*), parameter :: depthname = 'depth'

    ! Open the netCDF file
    retval = nf90_open(filename, nf90_nowrite, ncid)
    if (retval /= nf90_noerr) then
      print *, 'Failed to open file: ', trim(filename)
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Retrieve the dimension IDs and sizes for longitude, latitude, and depth
    ! Get the longitude dimension ID and size
    retval = nf90_inq_dimid(ncid, lonname, lonid)
    if (retval == nf90_noerr) then
      retval = nf90_inquire_dimension(ncid, lonid, len=lonsize)
    else
      print *, 'Failed to get ID and length of longitude dimension.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif
   !print *, 'Longitude size = ', lonsize

    ! Get the latitude dimension ID and size
    retval = nf90_inq_dimid(ncid, latname, latid)
    if (retval == nf90_noerr) then
      retval = nf90_inquire_dimension(ncid, latid, len=latsize)
    else
      print *, 'Failed to get ID and length of latitude dimension.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif
   !print *, 'Latitude size = ', latsize

    ! Get the depth dimension ID and size
    retval = nf90_inq_dimid(ncid, depthname, depthid)
    if (retval == nf90_noerr) then
      retval = nf90_inquire_dimension(ncid, depthid, len=depthsize)
    else
      print *, 'Failed to get ID and length of depth dimension.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif
   !print *, 'Depth size = ', depthsize

    ! Deallocate arrays if they are already allocated
    if (allocated(lon_values)) deallocate(lon_values)
    if (allocated(lat_values)) deallocate(lat_values)
    if (allocated(depth_values)) deallocate(depth_values)
    if (allocated(var)) deallocate(var)

    ! Allocate memory for the dimension values
    allocate(lon_values(lonsize))
    allocate(lat_values(latsize))
    allocate(depth_values(depthsize))

    ! Retrieve the longitude values
    retval = nf90_get_var(ncid, lonid, lon_values)
    if (retval /= nf90_noerr) then
      print *, 'Failed to read longitude values.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Retrieve the latitude values
    retval = nf90_get_var(ncid, latid, lat_values)
    if (retval /= nf90_noerr) then
      print *, 'Failed to read latitude values.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Retrieve the depth values
    retval = nf90_get_var(ncid, depthid, depth_values)
    if (retval /= nf90_noerr) then
      print *, 'Failed to read depth values.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Allocate memory for the 3D data variable based on retrieved dimensions
    allocate(var(lonsize, latsize, depthsize))
    if (.not. allocated(var)) then
      print *, 'Allocation of ', varname, ' failed.'
      stop
    endif

    ! Retrieve the variable ID and read the data into the allocated array
    retval = nf90_inq_varid(ncid, varname, varid)
    if (retval /= nf90_noerr) then
      print *, 'Failed to retrieve variable ID for ', varname, '.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif
    retval = nf90_get_var(ncid, varid, var)
    if (retval /= nf90_noerr) then
      print *, 'Failed to read data for variable ', varname, '.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Close the netCDF file
    retval = nf90_close(ncid)
    if (retval /= nf90_noerr) then
      print *, 'Error occurred while closing the file.'
      stop
    endif

  end subroutine read_netcdf_3d

!!subroutine write_netcdf_2d(filename, varname, varunit, var, lon_values, lat_values)
!!  implicit none

!!  character(len=*), intent(in) :: filename, varname, varunit
!!  double precision, dimension(:,:), intent(in) :: var
!!  double precision, dimension(:), intent(in) :: lon_values, lat_values

!!  integer(kind=4) :: ncid, varid, retval
!!  integer(kind=4) :: lonid, lonsize
!!  integer(kind=4) :: latid, latsize
!!  integer(kind=4) :: lonvarid, latvarid

!!  character(len=*), parameter :: lonname = 'lon'
!!  character(len=*), parameter :: latname = 'lat'

!!  double precision, parameter :: fill_value = 999.0

!!  ! Determine the size of each dimension from the input arrays
!!  lonsize = size(lon_values)
!!  latsize = size(lat_values)

!!  ! Create the netCDF file
!!  retval = nf90_create(filename, nf90_clobber, ncid)
!!  if (retval /= nf90_noerr) then
!!    print *, 'Failed to create file.'
!!    print *, 'Error: ', nf90_strerror(retval)
!!    stop
!!  endif

!!  ! Define the dimensions
!!  retval = nf90_def_dim(ncid, lonname, lonsize, lonid)
!!  if (retval /= nf90_noerr) then
!!    print *, 'Failed to define dimension ', lonname
!!    print *, 'Error: ', nf90_strerror(retval)
!!    stop
!!  endif

!!  retval = nf90_def_dim(ncid, latname, latsize, latid)
!!  if (retval /= nf90_noerr) then
!!    print *, 'Failed to define dimension ', latname
!!    print *, 'Error: ', nf90_strerror(retval)
!!    stop
!!  endif

!!  ! Define the dimension variables
!!  retval = nf90_def_var(ncid, lonname, nf90_double, [lonid], lonvarid)
!!  if (retval /= nf90_noerr) then
!!    print *, 'Failed to define variable ', lonname
!!    print *, 'Error: ', nf90_strerror(retval)
!!    stop
!!  endif

!!  retval = nf90_def_var(ncid, latname, nf90_double, [latid], latvarid)
!!  if (retval /= nf90_noerr) then
!!    print *, 'Failed to define variable ', latname
!!    print *, 'Error: ', nf90_strerror(retval)
!!    stop
!!  endif

!!  ! Add attributes to the dimension variables
!!  retval = nf90_put_att(ncid, lonvarid, 'standard_name', 'longitude')
!!  retval = nf90_put_att(ncid, lonvarid, 'long_name', 'longitude')
!!  retval = nf90_put_att(ncid, lonvarid, 'units', 'degrees_east')
!!  retval = nf90_put_att(ncid, lonvarid, 'axis', 'X')

!!  retval = nf90_put_att(ncid, latvarid, 'standard_name', 'latitude')
!!  retval = nf90_put_att(ncid, latvarid, 'long_name', 'latitude')
!!  retval = nf90_put_att(ncid, latvarid, 'units', 'degrees_north')
!!  retval = nf90_put_att(ncid, latvarid, 'axis', 'Y')

!!  ! Define the variable
!!  retval = nf90_def_var(ncid, varname, nf90_double, [lonid, latid], varid)
!!  if (retval /= nf90_noerr) then
!!    print *, 'Failed to define variable ', varname
!!    print *, 'Error: ', nf90_strerror(retval)
!!    stop
!!  endif

!!  ! Add attributes to the variable
!!  retval = nf90_put_att(ncid, varid, 'units', varunit)
!!  retval = nf90_put_att(ncid, varid, '_FillValue', fill_value)
!!  retval = nf90_put_att(ncid, varid, 'missing_value', fill_value)
!!  retval = nf90_put_att(ncid, varid, 'cell_methods', 'time: mean')
!!  retval = nf90_put_att(ncid, varid, 'original_name', varname)

!!  ! End define mode
!!  retval = nf90_enddef(ncid)
!!  if (retval /= nf90_noerr) then
!!    print *, 'Failed to end define mode.'
!!    print *, 'Error: ', nf90_strerror(retval)
!!    stop
!!  endif

!!  ! Write the data for the dimensions
!!  retval = nf90_put_var(ncid, lonvarid, lon_values)
!!  if (retval /= nf90_noerr) then
!!    print *, 'Failed to write data for ', lonname
!!    print *, 'Error: ', nf90_strerror(retval)
!!    stop
!!  endif

!!  retval = nf90_put_var(ncid, latvarid, lat_values)
!!  if (retval /= nf90_noerr) then
!!    print *, 'Failed to write data for ', latname
!!    print *, 'Error: ', nf90_strerror(retval)
!!    stop
!!  endif

!!  ! Write the data for the variable
!!  retval = nf90_put_var(ncid, varid, var)
!!  if (retval /= nf90_noerr) then
!!    print *, 'Failed to write data for variable ', varname, '.'
!!    print *, 'Error: ', nf90_strerror(retval)
!!    stop
!!  endif

!!  ! Close the netCDF file
!!  retval = nf90_close(ncid)
!!  if (retval /= nf90_noerr) then
!!    print *, 'Error occurred while closing the file.'
!!    stop
!!  endif

!!end subroutine write_netcdf_2d

!!subroutine write_netcdf_3d(filename, varname, varunit, var, lon_values, lat_values, depth_values)
!!  implicit none

!!  character(len=*), intent(in) :: filename, varname, varunit
!!  double precision, dimension(:,:,:), intent(in) :: var
!!  double precision, dimension(:), intent(in) :: lon_values, lat_values, depth_values

!!  integer(kind=4) :: ncid, varid, retval
!!  integer(kind=4) :: lonid, lonsize
!!  integer(kind=4) :: latid, latsize
!!  integer(kind=4) :: depid, depthsize
!!  integer(kind=4) :: lonvarid, latvarid, depvarid

!!  character(len=*), parameter :: lonname = 'lon'
!!  character(len=*), parameter :: latname = 'lat'
!!  character(len=*), parameter :: depname = 'depth'

!!  double precision, parameter :: fill_value = 999.0

!!  ! Determine the size of each dimension from the input arrays
!!  lonsize = size(lon_values)
!!  latsize = size(lat_values)
!!  depthsize = size(depth_values)

!!  ! Create the netCDF file
!!  retval = nf90_create(filename, nf90_clobber, ncid)
!!  if (retval /= nf90_noerr) then
!!    print *, 'Failed to create file.'
!!    print *, 'Error: ', nf90_strerror(retval)
!!    stop
!!  endif

!!  ! Define the dimensions
!!  retval = nf90_def_dim(ncid, lonname, lonsize, lonid)
!!  if (retval /= nf90_noerr) then
!!    print *, 'Failed to define dimension ', lonname
!!    print *, 'Error: ', nf90_strerror(retval)
!!    stop
!!  endif

!!  retval = nf90_def_dim(ncid, latname, latsize, latid)
!!  if (retval /= nf90_noerr) then
!!    print *, 'Failed to define dimension ', latname
!!    print *, 'Error: ', nf90_strerror(retval)
!!    stop
!!  endif

!!  retval = nf90_def_dim(ncid, depname, depthsize, depid)
!!  if (retval /= nf90_noerr) then
!!    print *, 'Failed to define dimension ', depname
!!    print *, 'Error: ', nf90_strerror(retval)
!!    stop
!!  endif

!!  ! Define the dimension variables
!!  retval = nf90_def_var(ncid, lonname, nf90_double, [lonid], lonvarid)
!!  if (retval /= nf90_noerr) then
!!    print *, 'Failed to define variable ', lonname
!!    print *, 'Error: ', nf90_strerror(retval)
!!    stop
!!  endif

!!  retval = nf90_def_var(ncid, latname, nf90_double, [latid], latvarid)
!!  if (retval /= nf90_noerr) then
!!    print *, 'Failed to define variable ', latname
!!    print *, 'Error: ', nf90_strerror(retval)
!!    stop
!!  endif

!!  retval = nf90_def_var(ncid, depname, nf90_double, [depid], depvarid)
!!  if (retval /= nf90_noerr) then
!!    print *, 'Failed to define variable ', depname
!!    print *, 'Error: ', nf90_strerror(retval)
!!    stop
!!  endif

!!  ! Add attributes to the dimension variables
!!  retval = nf90_put_att(ncid, lonvarid, 'standard_name', 'longitude')
!!  retval = nf90_put_att(ncid, lonvarid, 'long_name', 'longitude')
!!  retval = nf90_put_att(ncid, lonvarid, 'units', 'degrees_east')
!!  retval = nf90_put_att(ncid, lonvarid, 'axis', 'X')

!!  retval = nf90_put_att(ncid, latvarid, 'standard_name', 'latitude')
!!  retval = nf90_put_att(ncid, latvarid, 'long_name', 'latitude')
!!  retval = nf90_put_att(ncid, latvarid, 'units', 'degrees_north')
!!  retval = nf90_put_att(ncid, latvarid, 'axis', 'Y')

!!  retval = nf90_put_att(ncid, depvarid, 'standard_name', 'depth')
!!  retval = nf90_put_att(ncid, depvarid, 'long_name', 'depth')
!!  retval = nf90_put_att(ncid, depvarid, 'units', 'm')
!!  retval = nf90_put_att(ncid, depvarid, 'axis', 'Z')

!!  ! Define the variable
!!  retval = nf90_def_var(ncid, varname, nf90_double, [lonid, latid, depid], varid)
!!  if (retval /= nf90_noerr) then
!!    print *, 'Failed to define variable ', varname
!!    print *, 'Error: ', nf90_strerror(retval)
!!    stop
!!  endif

!!  ! Add attributes to the variable
!!  retval = nf90_put_att(ncid, varid, 'units', varunit)
!!  retval = nf90_put_att(ncid, varid, '_FillValue', fill_value)
!!  retval = nf90_put_att(ncid, varid, 'missing_value', fill_value)
!!  retval = nf90_put_att(ncid, varid, 'cell_methods', 'time: mean')
!!  retval = nf90_put_att(ncid, varid, 'original_name', varname)

!!  ! End define mode
!!  retval = nf90_enddef(ncid)
!!  if (retval /= nf90_noerr) then
!!    print *, 'Failed to end define mode.'
!!    print *, 'Error: ', nf90_strerror(retval)
!!    stop
!!  endif

!!  ! Write the data for the dimensions
!!  retval = nf90_put_var(ncid, lonvarid, lon_values)
!!  if (retval /= nf90_noerr) then
!!    print *, 'Failed to write data for ', lonname
!!    print *, 'Error: ', nf90_strerror(retval)
!!    stop
!!  endif

!!  retval = nf90_put_var(ncid, latvarid, lat_values)
!!  if (retval /= nf90_noerr) then
!!    print *, 'Failed to write data for ', latname
!!    print *, 'Error: ', nf90_strerror(retval)
!!    stop
!!  endif

!!  retval = nf90_put_var(ncid, depvarid, depth_values)
!!  if (retval /= nf90_noerr) then
!!    print *, 'Failed to write data for ', depname
!!    print *, 'Error: ', nf90_strerror(retval)
!!    stop
!!  endif

!!  ! Write the data for the variable
!!  retval = nf90_put_var(ncid, varid, var)
!!  if (retval /= nf90_noerr) then
!!    print *, 'Failed to write data for variable ', varname, '.'
!!    print *, 'Error: ', nf90_strerror(retval)
!!    stop
!!  endif

!!  ! Close the netCDF file
!!  retval = nf90_close(ncid)
!!  if (retval /= nf90_noerr) then
!!    print *, 'Error occurred while closing the file.'
!!    stop
!!  endif

!!end subroutine write_netcdf_3d

  subroutine write_netcdf_2d(filename, varname, varunit, var, lon_values, lat_values, time_value)
    implicit none
    character(len=*), intent(in) :: filename, varname, varunit
    double precision, dimension(:,:,:), intent(in) :: var
    double precision, dimension(:), intent(in) :: lon_values, lat_values
    integer(kind=4), intent(in) :: time_value
    double precision :: time_value_real
    integer(kind=4) :: ncid, varid, retval
    integer(kind=4) :: lonid, lonsize
    integer(kind=4) :: latid, latsize
    integer(kind=4) :: timeid, timesize
    integer(kind=4) :: lonvarid, latvarid, timevarid
    integer(kind=4) :: time_index
    logical :: file_exists
    character(len=*), parameter :: lonname = 'lon'
    character(len=*), parameter :: latname = 'lat'
    character(len=*), parameter :: timename = 'time'
    double precision, parameter :: fill_value = 999.0
    integer(kind=4), dimension(1) :: start_time, count_time
    integer(kind=4), dimension(3) :: start_var, count_var

    ! Determine the size of each dimension from the input arrays
    lonsize = size(lon_values)
    latsize = size(lat_values)
    time_value_real = real(time_value, kind=8)

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
      retval = nf90_inq_dimid(ncid, timename, timeid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to get time dimension.'
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      ! Get the current size of the time dimension
      retval = nf90_inquire_dimension(ncid, timeid, len=timesize)
      if (retval /= nf90_noerr) then
        print *, 'Failed to get time dimension size.'
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      ! Get the time variable ID
      retval = nf90_inq_varid(ncid, timename, timevarid)
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
      time_index = timesize + 1

    else
      ! Create the netCDF file
      retval = nf90_create(filename, nf90_clobber, ncid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to create file.'
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      ! Define the dimensions
      retval = nf90_def_dim(ncid, lonname, lonsize, lonid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to define dimension ', lonname
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      retval = nf90_def_dim(ncid, latname, latsize, latid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to define dimension ', latname
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      retval = nf90_def_dim(ncid, timename, nf90_unlimited, timeid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to define dimension ', timename
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      ! Get the current size of the time dimension
      retval = nf90_inquire_dimension(ncid, timeid, len=timesize)
      if (retval /= nf90_noerr) then
        print *, 'Failed to get time dimension size.'
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      ! Define the dimension variables
      retval = nf90_def_var(ncid, lonname, nf90_double, [lonid], lonvarid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to define variable ', lonname
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      retval = nf90_def_var(ncid, latname, nf90_double, [latid], latvarid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to define variable ', latname
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      retval = nf90_def_var(ncid, timename, nf90_double, [timeid], timevarid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to define variable ', timename
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      ! Add attributes to the dimension variables
      retval = nf90_put_att(ncid, lonvarid, 'standard_name', 'longitude')
      retval = nf90_put_att(ncid, lonvarid, 'long_name', 'longitude')
      retval = nf90_put_att(ncid, lonvarid, 'units', 'degrees_east')
      retval = nf90_put_att(ncid, lonvarid, 'axis', 'X')

      retval = nf90_put_att(ncid, latvarid, 'standard_name', 'latitude')
      retval = nf90_put_att(ncid, latvarid, 'long_name', 'latitude')
      retval = nf90_put_att(ncid, latvarid, 'units', 'degrees_north')
      retval = nf90_put_att(ncid, latvarid, 'axis', 'Y')

      retval = nf90_put_att(ncid, timevarid, 'standard_name', 'time')
      retval = nf90_put_att(ncid, timevarid, 'long_name', 'time')
      retval = nf90_put_att(ncid, timevarid, 'units', 'year')

      ! Define the variable
      retval = nf90_def_var(ncid, varname, nf90_double, [lonid, latid, timeid], varid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to define variable ', varname
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
      retval = nf90_put_var(ncid, lonvarid, lon_values)
      if (retval /= nf90_noerr) then
        print *, 'Failed to write data for ', lonname
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      retval = nf90_put_var(ncid, latvarid, lat_values)
      if (retval /= nf90_noerr) then
        print *, 'Failed to write data for ', latname
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      ! Initialize time index for the new file
      time_index = 1

      ! Write the initial time value
      retval = nf90_put_var(ncid, timevarid, time_value_real)
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
      retval = nf90_put_var(ncid, timevarid, time_value_real, start=start_time)
      if (retval /= nf90_noerr) then
        print *, 'Failed to write new time value.'
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif
    endif

    ! Write the data for the variable
    start_var = (/ 1, 1, time_index /)
    count_var = (/ lonsize, latsize, 1 /)
    retval = nf90_put_var(ncid, varid, var, start=start_var, count=count_var)
    if (retval /= nf90_noerr) then
      print *, 'Failed to write data for variable ', varname, '.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Close the netCDF file
    retval = nf90_close(ncid)
    if (retval /= nf90_noerr) then
      print *, 'Error occurred while closing the file.'
      stop
    endif

  end subroutine write_netcdf_2d

  subroutine write_netcdf_3d(filename, varname, varunit, var, lon_values, lat_values, depth_values, time_value)
    implicit none
    character(len=*), intent(in) :: filename, varname, varunit
    double precision, dimension(:,:,:,:), intent(in) :: var
    double precision, dimension(:), intent(in) :: lon_values, lat_values, depth_values
    integer(kind=4), intent(in) :: time_value
    double precision :: time_value_real
    integer(kind=4) :: ncid, varid, retval
    integer(kind=4) :: lonid, lonsize
    integer(kind=4) :: latid, latsize
    integer(kind=4) :: depid, depthsize
    integer(kind=4) :: timeid, timesize
    integer(kind=4) :: lonvarid, latvarid, depvarid, timevarid
    integer(kind=4) :: time_index
    logical :: file_exists
    character(len=*), parameter :: lonname = 'lon'
    character(len=*), parameter :: latname = 'lat'
    character(len=*), parameter :: depname = 'depth'
    character(len=*), parameter :: timename = 'time'
    double precision, parameter :: fill_value = 999.0
    integer(kind=4), dimension(1) :: start_time, count_time
    integer(kind=4), dimension(4) :: start_var, count_var

    ! Determine the size of each dimension from the input arrays
    lonsize = size(lon_values)
    latsize = size(lat_values)
    depthsize = size(depth_values)
    time_value_real = real(time_value, kind=8)

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
      retval = nf90_inq_dimid(ncid, timename, timeid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to get time dimension.'
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      ! Get the current size of the time dimension
      retval = nf90_inquire_dimension(ncid, timeid, len=timesize)
      if (retval /= nf90_noerr) then
        print *, 'Failed to get time dimension size.'
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      ! Get the time variable ID
      retval = nf90_inq_varid(ncid, timename, timevarid)
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
      time_index = timesize + 1

    else
      ! Create the netCDF file
      retval = nf90_create(filename, nf90_clobber, ncid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to create file.'
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      ! Define the dimensions
      retval = nf90_def_dim(ncid, lonname, lonsize, lonid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to define dimension ', lonname
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      retval = nf90_def_dim(ncid, latname, latsize, latid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to define dimension ', latname
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      retval = nf90_def_dim(ncid, depname, depthsize, depid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to define dimension ', depname
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      retval = nf90_def_dim(ncid, timename, nf90_unlimited, timeid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to define dimension ', timename
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      ! Get the current size of the time dimension
      retval = nf90_inquire_dimension(ncid, timeid, len=timesize)
      if (retval /= nf90_noerr) then
        print *, 'Failed to get time dimension size.'
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      ! Define the dimension variables
      retval = nf90_def_var(ncid, lonname, nf90_double, [lonid], lonvarid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to define variable ', lonname
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      retval = nf90_def_var(ncid, latname, nf90_double, [latid], latvarid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to define variable ', latname
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      retval = nf90_def_var(ncid, depname, nf90_double, [depid], depvarid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to define variable ', depname
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      retval = nf90_def_var(ncid, timename, nf90_double, [timeid], timevarid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to define variable ', timename
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      ! Add attributes to the dimension variables
      retval = nf90_put_att(ncid, lonvarid, 'standard_name', 'longitude')
      retval = nf90_put_att(ncid, lonvarid, 'long_name', 'longitude')
      retval = nf90_put_att(ncid, lonvarid, 'units', 'degrees_east')
      retval = nf90_put_att(ncid, lonvarid, 'axis', 'X')

      retval = nf90_put_att(ncid, latvarid, 'standard_name', 'latitude')
      retval = nf90_put_att(ncid, latvarid, 'long_name', 'latitude')
      retval = nf90_put_att(ncid, latvarid, 'units', 'degrees_north')
      retval = nf90_put_att(ncid, latvarid, 'axis', 'Y')

      retval = nf90_put_att(ncid, depvarid, 'standard_name', 'depth')
      retval = nf90_put_att(ncid, depvarid, 'long_name', 'depth')
      retval = nf90_put_att(ncid, depvarid, 'units', 'm')
      retval = nf90_put_att(ncid, depvarid, 'axis', 'Z')

      retval = nf90_put_att(ncid, timevarid, 'standard_name', 'time')
      retval = nf90_put_att(ncid, timevarid, 'long_name', 'time')
      retval = nf90_put_att(ncid, timevarid, 'units', 'year')

      ! Define the variable
      retval = nf90_def_var(ncid, varname, nf90_double, [lonid, latid, depid, timeid], varid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to define variable ', varname
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
      retval = nf90_put_var(ncid, lonvarid, lon_values)
      if (retval /= nf90_noerr) then
        print *, 'Failed to write data for ', lonname
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      retval = nf90_put_var(ncid, latvarid, lat_values)
      if (retval /= nf90_noerr) then
        print *, 'Failed to write data for ', latname
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      retval = nf90_put_var(ncid, depvarid, depth_values)
      if (retval /= nf90_noerr) then
        print *, 'Failed to write data for ', depname
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      ! Initialize time index for the new file
      time_index = 1

      ! Write the initial time value
      retval = nf90_put_var(ncid, timevarid, time_value_real)
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
      retval = nf90_put_var(ncid, timevarid, time_value_real, start=start_time)
      if (retval /= nf90_noerr) then
        print *, 'Failed to write new time value.'
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif
    endif

    ! Write the data for the variable
    start_var = (/ 1, 1, 1, time_index /)
    count_var = (/ lonsize, latsize, depthsize, 1 /)
    retval = nf90_put_var(ncid, varid, var, start=start_var, count=count_var)
    if (retval /= nf90_noerr) then
      print *, 'Failed to write data for variable ', varname, '.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Close the netCDF file
    retval = nf90_close(ncid)
    if (retval /= nf90_noerr) then
      print *, 'Error occurred while closing the file.'
      stop
    endif

  end subroutine write_netcdf_3d

end module netcdf_utilities
