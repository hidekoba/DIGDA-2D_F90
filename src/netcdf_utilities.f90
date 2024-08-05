module netcdf_utilities
    use netcdf
    implicit none
contains

  subroutine read_netcdf_2d(filename, varname, var, xvalues, yvalues)
    implicit none

    character(len=*), intent(in) :: filename, varname
    double precision, dimension(:,:), allocatable, intent(out) :: var

    integer(kind=4) :: ncid, varid, retval
    integer(kind=4) :: xid, xsize
    integer(kind=4) :: yid, ysize

    character(len=*), parameter :: xname = 'lon'
    character(len=*), parameter :: yname = 'lat'

    double precision, allocatable, dimension(:) :: xvalues, yvalues
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

    ! Deallocate arrays if they are already allocated
    if (allocated(xvalues)) deallocate(xvalues)
    if (allocated(yvalues)) deallocate(yvalues)
    if (allocated(var)) deallocate(var)

    ! Allocate memory for the dimension values
    allocate(xvalues(xsize))
    allocate(yvalues(ysize))

    ! Retrieve the longitude values
    retval = nf90_get_var(ncid, xid, xvalues)
    if (retval /= nf90_noerr) then
      print *, 'Failed to read longitude values.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif
   !print *, 'Longitude values:'
   !do i = 1, xsize
   !  print *, xvalues(i)
   !end do

    ! Retrieve the latitude values
    retval = nf90_get_var(ncid, yid, yvalues)
    if (retval /= nf90_noerr) then
      print *, 'Failed to read latitude values.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif
   !print *, 'Latitude values:'
   !do i = 1, ysize
   !  print *, yvalues(i)
   !end do

    ! Allocate memory for the 2D data variable based on retrieved dimensions
    allocate(var(xsize, ysize))
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

   !! Deallocate arrays if they are already allocated
   !if (allocated(xvalues)) deallocate(xvalues)
   !if (allocated(yvalues)) deallocate(yvalues)
   !if (allocated(var)) deallocate(var)

  end subroutine read_netcdf_2d


  subroutine read_netcdf_3d(filename, varname, var, xvalues, yvalues, zvalues)
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
    retval = nf90_get_var(ncid, xid, xvalues)
    if (retval /= nf90_noerr) then
      print *, 'Failed to read longitude values.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Retrieve the latitude values
    retval = nf90_get_var(ncid, yid, yvalues)
    if (retval /= nf90_noerr) then
      print *, 'Failed to read latitude values.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Retrieve the depth values
    retval = nf90_get_var(ncid, zid, zvalues)
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

   !! Deallocate arrays if they are already allocated
   !if (allocated(xvalues)) deallocate(xvalues)
   !if (allocated(yvalues)) deallocate(yvalues)
   !if (allocated(zvalues)) deallocate(zvalues)
   !if (allocated(var)) deallocate(var)

  end subroutine read_netcdf_3d

  subroutine read_netcdf_4d(filename, varname, var, xvalues, yvalues, zvalues, tvalues)
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
    retval = nf90_get_var(ncid, xid, xvalues)
    if (retval /= nf90_noerr) then
      print *, 'Failed to read longitude values.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Retrieve the latitude values
    retval = nf90_get_var(ncid, yid, yvalues)
    if (retval /= nf90_noerr) then
      print *, 'Failed to read latitude values.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Retrieve the depth values
    retval = nf90_get_var(ncid, zid, zvalues)
    if (retval /= nf90_noerr) then
      print *, 'Failed to read depth values.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Retrieve the time values
    retval = nf90_get_var(ncid, tid, tvalues)
    if (retval /= nf90_noerr) then
      print *, 'Failed to read depth values.'
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

   !! Deallocate arrays if they are already allocated
   !if (allocated(xvalues)) deallocate(xvalues)
   !if (allocated(yvalues)) deallocate(yvalues)
   !if (allocated(zvalues)) deallocate(zvalues)
   !if (allocated(tvalues)) deallocate(tvalues)
   !if (allocated(var)) deallocate(var)

  end subroutine read_netcdf_4d

!!subroutine write_netcdf_2d(filename, varname, varunit, var, xvalues, yvalues)
!!  implicit none

!!  character(len=*), intent(in) :: filename, varname, varunit
!!  double precision, dimension(:,:), intent(in) :: var
!!  double precision, dimension(:), intent(in) :: xvalues, yvalues

!!  integer(kind=4) :: ncid, varid, retval
!!  integer(kind=4) :: xid, xsize
!!  integer(kind=4) :: yid, ysize
!!  integer(kind=4) :: xvarid, yvarid

!!  character(len=*), parameter :: xname = 'lon'
!!  character(len=*), parameter :: yname = 'lat'

!!  double precision, parameter :: fill_value = 999.0

!!  ! Determine the size of each dimension from the input arrays
!!  xsize = size(xvalues)
!!  ysize = size(yvalues)

!!  ! Create the netCDF file
!!  retval = nf90_create(filename, nf90_clobber, ncid)
!!  if (retval /= nf90_noerr) then
!!    print *, 'Failed to create file.'
!!    print *, 'Error: ', nf90_strerror(retval)
!!    stop
!!  endif

!!  ! Define the dimensions
!!  retval = nf90_def_dim(ncid, xname, xsize, xid)
!!  if (retval /= nf90_noerr) then
!!    print *, 'Failed to define dimension ', xname
!!    print *, 'Error: ', nf90_strerror(retval)
!!    stop
!!  endif

!!  retval = nf90_def_dim(ncid, yname, ysize, yid)
!!  if (retval /= nf90_noerr) then
!!    print *, 'Failed to define dimension ', yname
!!    print *, 'Error: ', nf90_strerror(retval)
!!    stop
!!  endif

!!  ! Define the dimension variables
!!  retval = nf90_def_var(ncid, xname, nf90_double, [xid], xvarid)
!!  if (retval /= nf90_noerr) then
!!    print *, 'Failed to define variable ', xname
!!    print *, 'Error: ', nf90_strerror(retval)
!!    stop
!!  endif

!!  retval = nf90_def_var(ncid, yname, nf90_double, [yid], yvarid)
!!  if (retval /= nf90_noerr) then
!!    print *, 'Failed to define variable ', yname
!!    print *, 'Error: ', nf90_strerror(retval)
!!    stop
!!  endif

!!  ! Add attributes to the dimension variables
!!  retval = nf90_put_att(ncid, xvarid, 'standard_name', 'longitude')
!!  retval = nf90_put_att(ncid, xvarid, 'long_name', 'longitude')
!!  retval = nf90_put_att(ncid, xvarid, 'units', 'degrees_east')
!!  retval = nf90_put_att(ncid, xvarid, 'axis', 'X')

!!  retval = nf90_put_att(ncid, yvarid, 'standard_name', 'latitude')
!!  retval = nf90_put_att(ncid, yvarid, 'long_name', 'latitude')
!!  retval = nf90_put_att(ncid, yvarid, 'units', 'degrees_north')
!!  retval = nf90_put_att(ncid, yvarid, 'axis', 'Y')

!!  ! Define the variable
!!  retval = nf90_def_var(ncid, varname, nf90_double, [xid, yid], varid)
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
!!  retval = nf90_put_var(ncid, xvarid, xvalues)
!!  if (retval /= nf90_noerr) then
!!    print *, 'Failed to write data for ', xname
!!    print *, 'Error: ', nf90_strerror(retval)
!!    stop
!!  endif

!!  retval = nf90_put_var(ncid, yvarid, yvalues)
!!  if (retval /= nf90_noerr) then
!!    print *, 'Failed to write data for ', yname
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

!!subroutine write_netcdf_3d(filename, varname, varunit, var, xvalues, yvalues, zvalues)
!!  implicit none

!!  character(len=*), intent(in) :: filename, varname, varunit
!!  double precision, dimension(:,:,:), intent(in) :: var
!!  double precision, dimension(:), intent(in) :: xvalues, yvalues, zvalues

!!  integer(kind=4) :: ncid, varid, retval
!!  integer(kind=4) :: xid, xsize
!!  integer(kind=4) :: yid, ysize
!!  integer(kind=4) :: zid, zsize
!!  integer(kind=4) :: xvarid, yvarid, zvarid

!!  character(len=*), parameter :: xname = 'lon'
!!  character(len=*), parameter :: yname = 'lat'
!!  character(len=*), parameter :: depname = 'depth'

!!  double precision, parameter :: fill_value = 999.0

!!  ! Determine the size of each dimension from the input arrays
!!  xsize = size(xvalues)
!!  ysize = size(yvalues)
!!  zsize = size(zvalues)

!!  ! Create the netCDF file
!!  retval = nf90_create(filename, nf90_clobber, ncid)
!!  if (retval /= nf90_noerr) then
!!    print *, 'Failed to create file.'
!!    print *, 'Error: ', nf90_strerror(retval)
!!    stop
!!  endif

!!  ! Define the dimensions
!!  retval = nf90_def_dim(ncid, xname, xsize, xid)
!!  if (retval /= nf90_noerr) then
!!    print *, 'Failed to define dimension ', xname
!!    print *, 'Error: ', nf90_strerror(retval)
!!    stop
!!  endif

!!  retval = nf90_def_dim(ncid, yname, ysize, yid)
!!  if (retval /= nf90_noerr) then
!!    print *, 'Failed to define dimension ', yname
!!    print *, 'Error: ', nf90_strerror(retval)
!!    stop
!!  endif

!!  retval = nf90_def_dim(ncid, depname, zsize, zid)
!!  if (retval /= nf90_noerr) then
!!    print *, 'Failed to define dimension ', depname
!!    print *, 'Error: ', nf90_strerror(retval)
!!    stop
!!  endif

!!  ! Define the dimension variables
!!  retval = nf90_def_var(ncid, xname, nf90_double, [xid], xvarid)
!!  if (retval /= nf90_noerr) then
!!    print *, 'Failed to define variable ', xname
!!    print *, 'Error: ', nf90_strerror(retval)
!!    stop
!!  endif

!!  retval = nf90_def_var(ncid, yname, nf90_double, [yid], yvarid)
!!  if (retval /= nf90_noerr) then
!!    print *, 'Failed to define variable ', yname
!!    print *, 'Error: ', nf90_strerror(retval)
!!    stop
!!  endif

!!  retval = nf90_def_var(ncid, depname, nf90_double, [zid], zvarid)
!!  if (retval /= nf90_noerr) then
!!    print *, 'Failed to define variable ', depname
!!    print *, 'Error: ', nf90_strerror(retval)
!!    stop
!!  endif

!!  ! Add attributes to the dimension variables
!!  retval = nf90_put_att(ncid, xvarid, 'standard_name', 'longitude')
!!  retval = nf90_put_att(ncid, xvarid, 'long_name', 'longitude')
!!  retval = nf90_put_att(ncid, xvarid, 'units', 'degrees_east')
!!  retval = nf90_put_att(ncid, xvarid, 'axis', 'X')

!!  retval = nf90_put_att(ncid, yvarid, 'standard_name', 'latitude')
!!  retval = nf90_put_att(ncid, yvarid, 'long_name', 'latitude')
!!  retval = nf90_put_att(ncid, yvarid, 'units', 'degrees_north')
!!  retval = nf90_put_att(ncid, yvarid, 'axis', 'Y')

!!  retval = nf90_put_att(ncid, zvarid, 'standard_name', 'depth')
!!  retval = nf90_put_att(ncid, zvarid, 'long_name', 'depth')
!!  retval = nf90_put_att(ncid, zvarid, 'units', 'm')
!!  retval = nf90_put_att(ncid, zvarid, 'axis', 'Z')

!!  ! Define the variable
!!  retval = nf90_def_var(ncid, varname, nf90_double, [xid, yid, zid], varid)
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
!!  retval = nf90_put_var(ncid, xvarid, xvalues)
!!  if (retval /= nf90_noerr) then
!!    print *, 'Failed to write data for ', xname
!!    print *, 'Error: ', nf90_strerror(retval)
!!    stop
!!  endif

!!  retval = nf90_put_var(ncid, yvarid, yvalues)
!!  if (retval /= nf90_noerr) then
!!    print *, 'Failed to write data for ', yname
!!    print *, 'Error: ', nf90_strerror(retval)
!!    stop
!!  endif

!!  retval = nf90_put_var(ncid, zvarid, zvalues)
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

  subroutine write_netcdf_2d(filename, varname, varunit, var, xvalues, yvalues, tvalues)
    implicit none
    character(len=*), intent(in) :: filename, varname, varunit
    double precision, dimension(:,:,:), intent(in) :: var
    double precision, dimension(:), intent(in) :: xvalues, yvalues
    integer(kind=4), dimension(:), intent(in) :: tvalues
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
        print *, 'Failed to create file.'
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
        print *, 'Failed to define variable ', xname
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      retval = nf90_def_var(ncid, yname, nf90_double, [yid], yvarid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to define variable ', yname
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      retval = nf90_def_var(ncid, tname, nf90_double, [tid], tvarid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to define variable ', tname
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
      retval = nf90_put_var(ncid, xvarid, xvalues)
      if (retval /= nf90_noerr) then
        print *, 'Failed to write data for ', xname
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      retval = nf90_put_var(ncid, yvarid, yvalues)
      if (retval /= nf90_noerr) then
        print *, 'Failed to write data for ', yname
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

  subroutine write_netcdf_3d(filename, varname, varunit, var, xvalues, yvalues, zvalues, tvalues)
    implicit none
    character(len=*), intent(in) :: filename, varname, varunit
    double precision, dimension(:,:,:,:), intent(in) :: var
    double precision, dimension(:), intent(in) :: xvalues, yvalues, zvalues
    integer(kind=4), dimension(:), intent(in) :: tvalues
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
    character(len=*), parameter :: depname = 'depth'
    character(len=*), parameter :: tname = 'time'
    double precision, parameter :: fill_value = 999.0
    integer(kind=4), dimension(1) :: start_time, count_time
    integer(kind=4), dimension(4) :: start_var, count_var

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
        print *, 'Failed to create file.'
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

      retval = nf90_def_dim(ncid, depname, zsize, zid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to define dimension ', depname
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
        print *, 'Failed to define variable ', xname
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      retval = nf90_def_var(ncid, yname, nf90_double, [yid], yvarid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to define variable ', yname
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      retval = nf90_def_var(ncid, depname, nf90_double, [zid], zvarid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to define variable ', depname
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      retval = nf90_def_var(ncid, tname, nf90_double, [tid], tvarid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to define variable ', tname
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
      retval = nf90_put_var(ncid, xvarid, xvalues)
      if (retval /= nf90_noerr) then
        print *, 'Failed to write data for ', xname
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      retval = nf90_put_var(ncid, yvarid, yvalues)
      if (retval /= nf90_noerr) then
        print *, 'Failed to write data for ', yname
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      retval = nf90_put_var(ncid, zvarid, zvalues)
      if (retval /= nf90_noerr) then
        print *, 'Failed to write data for ', depname
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

  subroutine write_restart(filename, varname, varunit, var, xvalues, yvalues, zvalues, tvalues)
    use netcdf
    implicit none

    character(len=*), intent(in) :: filename, varname, varunit
    double precision, dimension(:,:,:,:), intent(in) :: var
    double precision, dimension(:), intent(in) :: xvalues, yvalues, zvalues
    integer(kind=4), dimension(:), intent(in) :: tvalues

    double precision, dimension(:), allocatable :: tvalues_real
    integer(kind=4) :: ncid, varid, retval
    integer(kind=4) :: xid, xsize
    integer(kind=4) :: yid, ysize
    integer(kind=4) :: zid, zsize
    integer(kind=4) :: tid, tsize
    integer(kind=4) :: xvarid, yvarid, zvarid, tvarid
    logical :: file_exists
    character(len=*), parameter :: xname = 'lon'
    character(len=*), parameter :: yname = 'lat'
    character(len=*), parameter :: zname = 'depth'
    character(len=*), parameter :: tname = 'time'
    double precision, parameter :: fill_value = 999.0
    integer(kind=4), dimension(4) :: start_var, count_var

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

      ! Check if the variable already exists
      retval = nf90_inq_varid(ncid, varname, varid)
      if (retval /= nf90_noerr) then
        ! If the variable does not exist, define it
        retval = nf90_redef(ncid)
        if (retval /= nf90_noerr) then
          print *, 'Failed to enter define mode.'
          print *, 'Error: ', nf90_strerror(retval)
          stop
        endif

        ! Define the data variable
        retval = nf90_def_var(ncid, varname, nf90_double, [xid, yid, zid, tid], varid)
        if (retval /= nf90_noerr) then
          print *, 'Failed to define variable ', varname
          print *, 'Error: ', nf90_strerror(retval)
          stop
        endif

        ! Add units attribute to the data variable
        retval = nf90_put_att(ncid, varid, 'units', varunit)
        if (retval /= nf90_noerr) then
          print *, 'Failed to add units attribute to variable ', varname
          print *, 'Error: ', nf90_strerror(retval)
          stop
        endif

        ! End define mode
        retval = nf90_enddef(ncid)
        if (retval /= nf90_noerr) then
          print *, 'Failed to end define mode.'
          print *, 'Error: ', nf90_strerror(retval)
          stop
        endif
      endif

    else
      ! Create the netCDF file
      retval = nf90_create(filename, nf90_clobber, ncid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to create file.'
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

      ! Define the dimension variables
      retval = nf90_def_var(ncid, xname, nf90_double, [xid], xvarid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to define variable ', xname
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      retval = nf90_def_var(ncid, yname, nf90_double, [yid], yvarid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to define variable ', yname
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      retval = nf90_def_var(ncid, zname, nf90_double, [zid], zvarid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to define variable ', zname
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      retval = nf90_def_var(ncid, tname, nf90_double, [tid], tvarid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to define variable ', tname
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      ! Define the data variable
      retval = nf90_def_var(ncid, varname, nf90_double, [xid, yid, zid, tid], varid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to define variable ', varname
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      ! Add units attribute to the data variable
      retval = nf90_put_att(ncid, varid, 'units', varunit)
      if (retval /= nf90_noerr) then
        print *, 'Failed to add units attribute to variable ', varname
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      ! End define mode
      retval = nf90_enddef(ncid)
      if (retval /= nf90_noerr) then
        print *, 'Failed to end define mode.'
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      ! Write the dimension variables
      retval = nf90_put_var(ncid, xvarid, xvalues)
      if (retval /= nf90_noerr) then
        print *, 'Failed to write longitude values.'
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      retval = nf90_put_var(ncid, yvarid, yvalues)
      if (retval /= nf90_noerr) then
        print *, 'Failed to write latitude values.'
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      retval = nf90_put_var(ncid, zvarid, zvalues)
      if (retval /= nf90_noerr) then
        print *, 'Failed to write depth values.'
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif

      retval = nf90_put_var(ncid, tvarid, tvalues_real)
      if (retval /= nf90_noerr) then
        print *, 'Failed to write time values.'
        print *, 'Error: ', nf90_strerror(retval)
        stop
      endif
    endif

    ! Write data to the variable
    start_var = [1, 1, 1, 1]
    count_var = [xsize, ysize, zsize, tsize]

    retval = nf90_put_var(ncid, varid, var, start=start_var, count=count_var)
    if (retval /= nf90_noerr) then
      print *, 'Failed to write data to variable.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    ! Close the NetCDF file
    retval = nf90_close(ncid)
    if (retval /= nf90_noerr) then
      print *, 'Error occurred while closing the file.'
      print *, 'Error: ', nf90_strerror(retval)
      stop
    endif

    print *, "Successfully wrote to file: ", trim(filename)

  end subroutine write_restart

end module netcdf_utilities
