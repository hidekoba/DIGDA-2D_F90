module gtool_utilities
    implicit none

contains

    subroutine read_gtool(filename, nsize, var)
        implicit none

        character(len=*), intent(in) :: filename
       !character(len=*), intent(in) :: varname
        integer(kind=4),intent(in) :: nsize
        real(8), allocatable, intent(out) :: var(:)
        character(len=16) :: chead(64)
        integer(kind=4) :: unit, iostat

        ! Variable declarations
        open(newunit=unit, file=filename, form='unformatted', status='old', iostat=iostat)
        write(*,*)'iostat' , iostat
        if (iostat /= 0) then
          print *, "Error opening file: ", filename
          stop
        endif

        ! Read header information
        read(unit, iostat=iostat) chead
        if (iostat /= 0) then
          print *, "Error reading header from file: ", filename
          stop
        endif

        ! Allocate the variable array
        allocate(var(nsize))

        ! Read the data into the variable array
        read(unit, iostat=iostat) var
        if (iostat /= 0) then
          print *, "Error reading data from file: ", filename
          stop
        endif

        ! close the file
        close(unit)
    end subroutine read_gtool

    subroutine write_gtool(filename, citem, ctitle, cunit, nx, ny, nz, var)
       implicit none

       character(len=*), intent(in) :: filename
       integer, intent(in) :: nx, ny, nz
       real(kind=8), intent(in) :: var(nx, ny, nz)
       character(len=16) :: chead(64)
       character(len=16) :: citem
       character(len=16) :: ctitle
       character(len=16) :: cunit
       integer :: unit, i, j, k
       real(kind=8) :: missing_value
       character(len=16) :: tmpStr

       ! Initialize header with blank spaces
       chead = ' '

       ! Example of setting up some header information
       missing_value = -999.0

       ! Setting up header values (as an example)
       chead(1) = '9010'
       chead(3) = citem
       chead(12) = '0'
       chead(13) = '1'
       chead(14) = ctitle
       chead(16) = cunit
       chead(25) = '0'
       chead(27) = ' 09990101 000000 '
       chead(28) = '0'
       chead(29) = 'TX'
       chead(30) = '1'
       write(tmpStr, '(I10)') nx
       chead(31) = trim(adjustl(tmpStr))
       chead(32) = 'TY'
       chead(33) = '1'
       write(tmpStr, '(I10)') ny
       chead(34) = trim(adjustl(tmpStr))
       chead(35) = 'TZ'
       chead(36) = '1'
       write(tmpStr, '(I10)') nz
       chead(37) = trim(adjustl(tmpStr))
       chead(38) = 'UR8'
       chead(50) = '09990101000000'
       write(tmpStr, '(E15.7)') missing_value
       chead(39:43) = trim(adjustl(tmpStr))
       chead(44) = '1'
       chead(46) = '0'
       chead(47) = '0.'
       chead(64) = '0'

       ! Open file for writing
       open(newunit=unit, file=filename, form='unformatted', status='replace', action='write')

       ! Write header and data
       write(unit) chead
      !write(unit) (((real(var(i,j,k)),i=1,nx),j=1,ny),k=1,nz)
       write(unit) var

       ! Close the file
       close(unit)

    end subroutine write_gtool

end module gtool_utilities
