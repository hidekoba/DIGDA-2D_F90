module netcdf_utilities
    use netcdf
    implicit none
contains

    subroutine read_netcdf_2d(filename, varname, var)
        implicit none

        character(len=*), intent(in) :: filename, varname
        real(kind=8), dimension(:,:), allocatable, intent(out) :: var

        integer(kind=4) :: ncid, varid, retval
        integer(kind=4) :: lonid, lonsize
        integer(kind=4) :: latid, latsize

        character(len=*), parameter :: lonname = 'lon'
        character(len=*), parameter :: latname = 'lat'

        integer(kind=4) :: i, j, k

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
        else if (retval /= nf90_noerr) then
            print *, 'Failed to get ID and length of longitude dimension.'
            print *, 'Error: ', nf90_strerror(retval)
            lonsize = 1
        endif
        print *, 'Longitude size = ', lonsize

        ! Get the latitude dimension ID and size
        retval = nf90_inq_dimid(ncid, latname, latid)
        if (retval == nf90_noerr) then
            retval = nf90_inquire_dimension(ncid, latid, len=latsize)
        else if (retval /= nf90_noerr) then
            print *, 'Failed to get ID and length of latitude dimension.'
            print *, 'Error: ', nf90_strerror(retval)
            latsize = 1
        endif
        print *, 'Latitude size = ', latsize

        ! Allocate memory for the 3D data variable based on retrieved dimensions
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

    subroutine read_netcdf_3d(filename, varname, var)
        implicit none

        character(len=*), intent(in) :: filename, varname
        real(kind=8), dimension(:,:,:), allocatable, intent(out) :: var

        integer(kind=4) :: ncid, varid, retval
        integer(kind=4) :: lonid, lonsize
        integer(kind=4) :: latid, latsize
        integer(kind=4) :: depid, depsize

        character(len=*), parameter :: lonname = 'lon'
        character(len=*), parameter :: latname = 'lat'
        character(len=*), parameter :: depname = 'dep'

        integer(kind=4) :: i, j, k

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
        else if (retval /= nf90_noerr) then
            print *, 'Failed to get ID and length of longitude dimension.'
            print *, 'Error: ', nf90_strerror(retval)
            lonsize = 1
        endif
        print *, 'Longitude size = ', lonsize

        ! Get the latitude dimension ID and size
        retval = nf90_inq_dimid(ncid, latname, latid)
        if (retval == nf90_noerr) then
            retval = nf90_inquire_dimension(ncid, latid, len=latsize)
        else if (retval /= nf90_noerr) then
            print *, 'Failed to get ID and length of latitude dimension.'
            print *, 'Error: ', nf90_strerror(retval)
            latsize = 1
        endif
        print *, 'Latitude size = ', latsize

        ! Get the depth dimension ID and size
        retval = nf90_inq_dimid(ncid, depname, depid)
        if (retval == nf90_noerr) then
            retval = nf90_inquire_dimension(ncid, depid, len=depsize)
        else if (retval /= nf90_noerr) then
            print *, 'Failed to get ID and length of depth dimension.'
            print *, 'Error: ', nf90_strerror(retval)
            depsize = 1
        endif
        print *, 'Depth size = ', depsize

        ! Allocate memory for the 3D data variable based on retrieved dimensions
        allocate(var(lonsize, latsize, depsize))
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

    subroutine write_netcdf_2d(filename, varname, var)
        implicit none

        character(len=*), intent(in) :: filename, varname
        real(kind=8), dimension(:,:), intent(in) :: var

        integer(kind=4) :: ncid, varid, retval
        integer(kind=4) :: lonid, lonsize
        integer(kind=4) :: latid, latsize

        character(len=*), parameter :: lonname = 'lon'
        character(len=*), parameter :: latname = 'lat'

        ! Determine the size of each dimension from the input variable
        lonsize = size(var, 1)
        latsize = size(var, 2)

        ! Create the netCDF file
        retval = nf90_create(filename, nf90_clobber, ncid)
        if (retval /= nf90_noerr) then
            print *, 'Failed to create file.'
            print *, 'Error: ', nf90_strerror(retval)
            stop
        endif

        ! Define the dimensions
        retval = nf90_def_dim(ncid, lonname, lonsize, lonid)
        retval = nf90_def_dim(ncid, latname, latsize, latid)

        ! Define the variable
        retval = nf90_def_var(ncid, varname, nf90_double, [lonid, latid], varid)

        ! End define mode
        retval = nf90_enddef(ncid)

        ! Write the data
        retval = nf90_put_var(ncid, varid, var)
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

    subroutine write_netcdf_3d(filename, varname, var)
        implicit none

        character(len=*), intent(in) :: filename, varname
        real(kind=8), dimension(:,:,:), intent(in) :: var

        integer(kind=4) :: ncid, varid, retval
        integer(kind=4) :: lonid, lonsize
        integer(kind=4) :: latid, latsize
        integer(kind=4) :: depid, depsize

        character(len=*), parameter :: lonname = 'lon'
        character(len=*), parameter :: latname = 'lat'
        character(len=*), parameter :: depname = 'depth'

        ! Determine the size of each dimension from the input variable
        lonsize = size(var, 1)
        latsize = size(var, 2)
        depsize = size(var, 3)

        ! Create the netCDF file
        retval = nf90_create(filename, nf90_clobber, ncid)
        if (retval /= nf90_noerr) then
            print *, 'Failed to create file.'
            print *, 'Error: ', nf90_strerror(retval)
            stop
        endif

        ! Define the dimensions
        retval = nf90_def_dim(ncid, lonname, lonsize, lonid)
        retval = nf90_def_dim(ncid, latname, latsize, latid)
        retval = nf90_def_dim(ncid, depname, depsize, depid)

        ! Define the variable
        retval = nf90_def_var(ncid, varname, nf90_double, [lonid, latid, depid], varid)

        ! End define mode
        retval = nf90_enddef(ncid)

        ! Write the data
        retval = nf90_put_var(ncid, varid, var)
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
