module gtool_utilities
    implicit none

contains

  subroutine read_gtool(filename, nsize, ntsize, chead, var)
    implicit none

    character(len=*), intent(in) :: filename
    integer,intent(in) :: nsize, ntsize
    double precision, allocatable, intent(out) :: var(:,:)
    character(len=16), intent(out) :: chead(64, ntsize)
    integer :: unit, iostat
    integer :: t

    ! Variable declarations
    open(newunit=unit, file=filename, form='unformatted', status='old', iostat=iostat)
    write(*,*)'iostat' , iostat
    if (iostat /= 0) then
      print *, "Error opening file: ", filename
      stop
    endif

    ! Allocate the variable array
    allocate(var(nsize, ntsize))

    do t = 1, ntsize

      ! Read header information
      read(unit, iostat=iostat) chead(:,t)
      if (iostat /= 0) then
        print *, "Error reading header from file: ", filename
        stop
      endif

      ! Read the data into the variable array
      read(unit, iostat=iostat) var(:,t)
      if (iostat /= 0) then
        print *, "Error reading data from file: ", filename
        stop
      endif

    end do

    ! close the file
    close(unit)
  end subroutine read_gtool

  subroutine write_gtool(filename, chead_orig, citem, ctitle, cunit, nx, ny, nz, nt, var)
    implicit none

    character(len=*), intent(in) :: filename
    integer, intent(in) :: nx, ny, nz, nt
    double precision, intent(in) :: var(nx*ny*nz, nt)
    character(len=16), intent(in) :: chead_orig(64, nt)
    character(len=16) :: chead(64)
    character(len=16), intent(in) :: citem
    character(len=16), intent(in) :: ctitle
    character(len=16), intent(in) :: cunit
    integer :: unit, i, j, k, t
    double precision :: missing_value
    character(len=16) :: tmpstr
    logical :: file_first

    file_first = .true.

    ! Open file for writing
    if (file_first) then
      open(newunit=unit, file=filename, form='unformatted', status='replace', action='write')
      file_first = .false.
    else
      open(newunit=unit, file=filename, form='unformatted', status='old', action='write', position='append')
    end if

    do t = 1, nt

    ! Initialize header with blank spaces
    chead = ' '

    ! Initialize header
   !chead = chead_orig(:, t)

    ! Example of setting up some header information
    missing_value = -999.0

    ! Setting up header values (as an example)
    chead(1) = '            9010'
    chead(3) = citem
    chead(14) = ctitle
    chead(16) = cunit
   !chead(27) = ' 09990101 000000 '
    chead(27) = chead_orig(27, t)
    chead(29) = '              TX'
    chead(30) = '               1'
    write(tmpstr, '(I10)') nx
    chead(31) = trim(adjustl(tmpstr))
    chead(32) = '              TY'
    chead(33) = '               1'
    write(tmpstr, '(I10)') ny
    chead(34) = trim(adjustl(tmpstr))
    chead(35) = '              TZ'
    chead(36) = '               1'
    write(tmpstr, '(I10)') nz
    chead(37) = trim(adjustl(tmpstr))
    chead(38) = 'UR8'
   !chead(50) = '09990101000000'
    chead(50) = chead_orig(50, t)
    write(tmpstr, '(E15.7)') missing_value
    chead(39:43) = trim(adjustl(tmpstr))
    chead(44) = '1'
    chead(46) = '0'
    chead(47) = '0.'
    chead(64) = '0'


    ! Write header and data
    write(unit) chead
   !write(unit) (((real(var(i,j,k)),i=1,nx),j=1,ny),k=1,nz)
    write(unit) var(:, t)

    end do

    ! Close the file
    close(unit)

  end subroutine write_gtool

  subroutine write_gtool_beta( filename, chead_orig, citem, ctitle,  cunit, nx, ny, nz, nt, var)
    implicit none

    character(len=*), intent(in) :: filename
    integer, intent(in) :: nx, ny, nz, nt
    double precision, intent(in) :: var(nx*ny, nz, nt)
    double precision :: vartmp(nx*ny*nz, nt)
    character(len=16), intent(in) :: chead_orig(64, nt)
    character(len=16) :: chead(64)
    character(len=16), intent(in) :: citem
    character(len=16), intent(in) :: ctitle
    character(len=16), intent(in) :: cunit
    integer :: unit, i, j, k, t, ij, ijk
    logical :: file_first

    file_first = .true.

    ! Open file for writing
    if (file_first) then
      open(newunit=unit, file=filename, form='unformatted', status='replace', action='write')
      file_first = .false.
    else
      open(newunit=unit, file=filename, form='unformatted', status='old', action='write', position='append')
    end if

    ! Initialize header with blank spaces
    chead = ' '

    do t = 1, nt

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          ijk = (k-1) * nx * ny + (j-1) * nx + i
          ij = (j-1) * nx + i
          vartmp(ijk, t) = var(ij, k, t)
        end do
      end do
    end do

    chead(1) = '            9010'  !!  1 : IDFM  : I10   :   GTOOL3 version                 : 9010
    chead(2) = ''                  !!  2 : DSET  : A16   :   Dataset name                   :
    chead(3) = citem               !!  3 : ITEM  : A16   :   Identification (variable) name :
    chead(4) = ''                  !!  4 : EDIT1 : A16   :   Edit abbreviation 01           :
    chead(5) = ''                  !!  5 : EDIT2 : A16   :   Edit abbreviation 02           :
    chead(6) = ''                  !!  6 : EDIT3 : A16   :   Edit abbreviation 03           :
    chead(7) = ''                  !!  7 : EDIT4 : A16   :   Edit abbreviation 04           :
    chead(8) = ''                  !!  8 : EDIT5 : A16   :   Edit abbreviation 05           :
    chead(9) = ''                  !!  9 : EDIT6 : A16   :   Edit abbreviation 06           :
    chead(10) = ''                 !! 10 : EDIT7 : A16   :   Edit abbreviation 07           :
    chead(11) = ''                 !! 11 : EDIT8 : A16   :   Edit abbreviation 08           :
    chead(12) = ''                 !! 12 : FNUM  : I10   :   File number                    : 0
    chead(13) = ''                 !! 13 : DNUM  : I10   :   Data number                    : 0
    chead(14) = ctitle             !! 14 : TITL1 : A16   :   Title 1                        :
    chead(15) = ''                 !! 15 : TITL2 : A16   :   Title 2                        :
    chead(16) = cunit              !! 16 : UNIT  : A16   :   Unit                           :
    chead(17) = ''                 !! 17 : ETTL1 : A16   :   Edit title 01                  :
    chead(18) = ''                 !! 18 : ETTL2 : A16   :   Edit title 02                  :
    chead(19) = ''                 !! 19 : ETTL3 : A16   :   Edit title 03                  :
    chead(20) = ''                 !! 20 : ETTL4 : A16   :   Edit title 04                  :
    chead(21) = ''                 !! 21 : ETTL5 : A16   :   Edit title 05                  :
    chead(22) = ''                 !! 22 : ETTL6 : A16   :   Edit title 06                  :
    chead(23) = ''                 !! 23 : ETTL7 : A16   :   Edit title 07                  :
    chead(24) = ''                 !! 24 : ETTL8 : A16   :   Edit title 08                  :
    chead(25) = ''                 !! 25 : TIME  : I10   :   Total time                     : 0
    chead(26) = ''                 !! 26 : UTIM  : A16   :   Time unit                      :
!   chead(27) = ''                 !! 27 : DATE  : A16   :   Time (yyyymmddhhmmss)          :
    chead(28) = ''                 !! 28 : TDUR  : I10   :   Data representative time       : 0
!   chead(29) = ''                 !! 29 : AITM1 : A16   :   Axis grid identification name  :
!   chead(30) = ''                 !! 30 : ASTR1 : I10   :   Start of axis grid number      : 1
!   chead(31) = ''                 !! 31 : AEND1 : I10   :   End of axis grid number        :
!   chead(32) = ''                 !! 32 : AITM2 : A16   :   Axis grid identification name  :
!   chead(33) = ''                 !! 33 : ASTR2 : I10   :   Start of axis grid number      : 1
!   chead(34) = ''                 !! 34 : AEND2 : I10   :   End of axis grid number        :
!   chead(35) = ''                 !! 35 : AITM3 : A16   :   Axis grid identification name  :
!   chead(36) = ''                 !! 36 : ASTR3 : I10   :   Start of axis grid number      : 1
!   chead(37) = ''                 !! 37 : AEND3 : I10   :   End of axis grid number        :
!   chead(38) = ''                 !! 38 : DFMT  : A16   :   Data format                    : UR4
!   chead(39) = ''                 !! 39 : MISS  : E15.7 :   Missing value                  : -999.
!   chead(40) = ''                 !! 40 : DMIN  : E15.7 :   Range (minimum)                : -999.
!   chead(41) = ''                 !! 41 : DMAX  : E15.7 :   Range (maximum)                : -999.
!   chead(42) = ''                 !! 42 : DIVS  : E15.7 :   Interval (small)               : -999.
!   chead(43) = ''                 !! 43 : DIVL  : E15.7 :   Interval (Large)               : -999.
!   chead(44) = ''                 !! 44 : STYP  : I10   :   Scaling type                   : 1
    chead(45) = ''                 !! 45 : COPTN : A16   :   Option (character)             :
    chead(46) = ''                 !! 46 : IOPTN : I16   :   Option (integer)               : 0
    chead(47) = ''                 !! 47 : ROPTN : E16.7 :   Option (real)                  : 0
    chead(48) = ''                 !! 48 : DATE1 : I16   :   Time (auxiliary)               :
    chead(49) = ''                 !! 49 : DATE2 : A16   :   Time unit (auxiliary)          :
!   chead(50) = ''                 !! 50 : MEMO1 : A16   :   Memo O1                        :
    chead(51) = ''                 !! 51 : MEMO2 : A16   :   Memo O2                        :
    chead(52) = ''                 !! 52 : MEMO3 : A16   :   Memo O3                        :
    chead(53) = ''                 !! 53 : MEMO4 : A16   :   Memo O4                        :
    chead(54) = ''                 !! 54 : MEMO5 : A16   :   Memo O5                        :
    chead(55) = ''                 !! 55 : MEMO6 : A16   :   Memo O6                        :
    chead(56) = ''                 !! 56 : MEMO7 : A16   :   Memo O7                        :
    chead(57) = ''                 !! 57 : MEMO8 : A16   :   Memo O8                        :
    chead(58) = ''                 !! 58 : MEMO9 : A16   :   Memo O9                        :
    chead(59) = ''                 !! 59 : MEM1O : A16   :   Memo 1O                        :
    chead(60) = ''                 !! 60 : CDATE : A16   :   Data creation date             :
    chead(61) = ''                 !! 61 : CSIGN : A16   :   Data creator                   :
    chead(62) = ''                 !! 62 : MDATE : A16   :   Data change date               :
    chead(63) = ''                 !! 63 : MSIGN : A16   :   Data modified by               :
    chead(64) = ''                 !! 64 : SIZE  : I10   :   Array size                     : 0
    write(chead(39), '(E15.7)') -999.0
    write(chead(40), '(E15.7)') -999.0
    write(chead(41), '(E15.7)') -999.0
    write(chead(42), '(E15.7)') -999.0
    write(chead(43), '(E15.7)') -999.0
    write(chead(44), '(I10)') 1
    write(chead(29), '(A16)') 'TX'
    write(chead(32), '(A16)') 'TY'
    write(chead(35), '(A16)') 'TZ'
    write(chead(30), '(I10)') 1
    write(chead(31), '(I10)') nx
    write(chead(33), '(I10)') 1
    write(chead(34), '(I10)') ny
    write(chead(36), '(I10)') 1
    write(chead(37), '(I10)') nz
    write(chead(64), '(I10)') nx*ny*nz
   !chead(27) = '09990101 000000 '
   !chead(50) = '  09990101000000'
    chead(27) = chead_orig(27, t)
    chead(50) = chead_orig(50, t)
    chead(38) = 'UR8'

    write(unit) chead
   !write(unit) (((REAL(var(i,j,k)),i=1,nx),j=1,ny),k=1,nz)
   !write(unit) var(:, t)
    write(unit) vartmp(:, t)

    end do

    close(unit)

  end subroutine write_gtool_beta

end module gtool_utilities
