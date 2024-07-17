module gtool_utilities
    implicit none

contains

  subroutine read_gtool(filename, nsize, var)
    implicit none

    character(len=*), intent(in) :: filename
   !character(len=*), intent(in) :: varname
    integer,intent(in) :: nsize
    double precision, allocatable, intent(out) :: var(:)
    character(len=16) :: chead(64)
    integer :: unit, iostat

    ! Variable declarations
    open(newunit=unit, file=filename, form='unformatted', status='old', iostat=iostat)
    WRITE(*,*)'iostat' , iostat
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
    double precision, intent(in) :: var(nx, ny, nz)
    character(len=16) :: chead(64)
    character(len=16) :: citem
    character(len=16) :: ctitle
    character(len=16) :: cunit
    integer :: unit, i, j, k
    double precision :: missing_value
    character(len=16) :: tmpstr

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
    WRITE(tmpstr, '(I10)') nx
    chead(31) = trim(adjustl(tmpstr))
    chead(32) = 'TY'
    chead(33) = '1'
    WRITE(tmpstr, '(I10)') ny
    chead(34) = trim(adjustl(tmpstr))
    chead(35) = 'TZ'
    chead(36) = '1'
    WRITE(tmpstr, '(I10)') nz
    chead(37) = trim(adjustl(tmpstr))
    chead(38) = 'UR8'
    chead(50) = '09990101000000'
    WRITE(tmpstr, '(E15.7)') missing_value
    chead(39:43) = trim(adjustl(tmpstr))
    chead(44) = '1'
    chead(46) = '0'
    chead(47) = '0.'
    chead(64) = '0'

    ! Open file for writing
    open(newunit=unit, file=filename, form='unformatted', status='replace', action='write')

    ! Write header and data
    WRITE(unit) chead
   !WRITE(unit) (((real(var(i,j,k)),i=1,nx),j=1,ny),k=1,nz)
    WRITE(unit) var

    ! Close the file
    close(unit)

  end subroutine write_gtool

  subroutine write_gtoolfile( &
 &           YY,    DAT3,  CFILE, CTITLE,  CUNIT, &
 &           NX,    NY,    NZ     )
    implicit none

    integer :: I, J, K
    integer :: NX, NY, NZ
    integer :: YY
    double precision :: DAT3(NX,NY,NZ)
    character(len=128) :: CFILE, CTITLE, CUNIT
    character(len=16) :: CHEAD(64)
    character(len=8) :: YEAR

    DATA CHEAD / 64*'' /
    WRITE(YEAR, 900) YY
900 FORMAT(I8.8)

      CHEAD(1) = '            9010'  !!  1 : IDFM  : I10   :   GTOOL3 version                 : 9010
      CHEAD(2) = ''                  !!  2 : DSET  : A16   :   Dataset name                   :
      CHEAD(3) = CTITLE              !!  3 : ITEM  : A16   :   Identification (variable) name :
      CHEAD(4) = ''                  !!  4 : EDIT1 : A16   :   Edit abbreviation 01           :
      CHEAD(5) = ''                  !!  5 : EDIT2 : A16   :   Edit abbreviation 02           :
      CHEAD(6) = ''                  !!  6 : EDIT3 : A16   :   Edit abbreviation 03           :
      CHEAD(7) = ''                  !!  7 : EDIT4 : A16   :   Edit abbreviation 04           :
      CHEAD(8) = ''                  !!  8 : EDIT5 : A16   :   Edit abbreviation 05           :
      CHEAD(9) = ''                  !!  9 : EDIT6 : A16   :   Edit abbreviation 06           :
      CHEAD(10) = ''                 !! 10 : EDIT7 : A16   :   Edit abbreviation 07           :
      CHEAD(11) = ''                 !! 11 : EDIT8 : A16   :   Edit abbreviation 08           :
      CHEAD(12) = ''                 !! 12 : FNUM  : I10   :   File number                    : 0
      CHEAD(13) = ''                 !! 13 : DNUM  : I10   :   Data number                    : 0
      CHEAD(14) = CTITLE             !! 14 : TITL1 : A16   :   Title 1                        :
      CHEAD(15) = ''                 !! 15 : TITL2 : A16   :   Title 2                        :
      CHEAD(16) = CUNIT              !! 16 : UNIT  : A16   :   Unit                           :
      CHEAD(17) = ''                 !! 17 : ETTL1 : A16   :   Edit title 01                  :
      CHEAD(18) = ''                 !! 18 : ETTL2 : A16   :   Edit title 02                  :
      CHEAD(19) = ''                 !! 19 : ETTL3 : A16   :   Edit title 03                  :
      CHEAD(20) = ''                 !! 20 : ETTL4 : A16   :   Edit title 04                  :
      CHEAD(21) = ''                 !! 21 : ETTL5 : A16   :   Edit title 05                  :
      CHEAD(22) = ''                 !! 22 : ETTL6 : A16   :   Edit title 06                  :
      CHEAD(23) = ''                 !! 23 : ETTL7 : A16   :   Edit title 07                  :
      CHEAD(24) = ''                 !! 24 : ETTL8 : A16   :   Edit title 08                  :
      CHEAD(25) = ''                 !! 25 : TIME  : I10   :   Total time                     : 0
      CHEAD(26) = ''                 !! 26 : UTIM  : A16   :   Time unit                      :
!     CHEAD(27) = ''                 !! 27 : DATE  : A16   :   Time (yyyymmddhhmmss)          :
      CHEAD(28) = ''                 !! 28 : TDUR  : I10   :   Data representative time       : 0
!     CHEAD(29) = ''                 !! 29 : AITM1 : A16   :   Axis grid identification name  :
!     CHEAD(30) = ''                 !! 30 : ASTR1 : I10   :   Start of axis grid number      : 1
!     CHEAD(31) = ''                 !! 31 : AEND1 : I10   :   End of axis grid number        :
!     CHEAD(32) = ''                 !! 32 : AITM2 : A16   :   Axis grid identification name  :
!     CHEAD(33) = ''                 !! 33 : ASTR2 : I10   :   Start of axis grid number      : 1
!     CHEAD(34) = ''                 !! 34 : AEND2 : I10   :   End of axis grid number        :
!     CHEAD(35) = ''                 !! 35 : AITM3 : A16   :   Axis grid identification name  :
!     CHEAD(36) = ''                 !! 36 : ASTR3 : I10   :   Start of axis grid number      : 1
!     CHEAD(37) = ''                 !! 37 : AEND3 : I10   :   End of axis grid number        :
!     CHEAD(38) = ''                 !! 38 : DFMT  : A16   :   Data format                    : UR4
!     CHEAD(39) = ''                 !! 39 : MISS  : E15.7 :   Missing value                  : -999.
!     CHEAD(40) = ''                 !! 40 : DMIN  : E15.7 :   Range (minimum)                : -999.
!     CHEAD(41) = ''                 !! 41 : DMAX  : E15.7 :   Range (maximum)                : -999.
!     CHEAD(42) = ''                 !! 42 : DIVS  : E15.7 :   Interval (small)               : -999.
!     CHEAD(43) = ''                 !! 43 : DIVL  : E15.7 :   Interval (Large)               : -999.
!     CHEAD(44) = ''                 !! 44 : STYP  : I10   :   Scaling type                   : 1
      CHEAD(45) = ''                 !! 45 : COPTN : A16   :   Option (character)             :
      CHEAD(46) = ''                 !! 46 : IOPTN : I16   :   Option (integer)               : 0
      CHEAD(47) = ''                 !! 47 : ROPTN : E16.7 :   Option (real)                  : 0
      CHEAD(48) = ''                 !! 48 : DATE1 : I16   :   Time (auxiliary)               :
      CHEAD(49) = ''                 !! 49 : DATE2 : A16   :   Time unit (auxiliary)          :
!     CHEAD(50) = ''                 !! 50 : MEMO1 : A16   :   Memo O1                        :
      CHEAD(51) = ''                 !! 51 : MEMO2 : A16   :   Memo O2                        :
      CHEAD(52) = ''                 !! 52 : MEMO3 : A16   :   Memo O3                        :
      CHEAD(53) = ''                 !! 53 : MEMO4 : A16   :   Memo O4                        :
      CHEAD(54) = ''                 !! 54 : MEMO5 : A16   :   Memo O5                        :
      CHEAD(55) = ''                 !! 55 : MEMO6 : A16   :   Memo O6                        :
      CHEAD(56) = ''                 !! 56 : MEMO7 : A16   :   Memo O7                        :
      CHEAD(57) = ''                 !! 57 : MEMO8 : A16   :   Memo O8                        :
      CHEAD(58) = ''                 !! 58 : MEMO9 : A16   :   Memo O9                        :
      CHEAD(59) = ''                 !! 59 : MEM1O : A16   :   Memo 1O                        :
      CHEAD(60) = ''                 !! 60 : CDATE : A16   :   Data creation date             :
      CHEAD(61) = ''                 !! 61 : CSIGN : A16   :   Data creator                   :
      CHEAD(62) = ''                 !! 62 : MDATE : A16   :   Data change date               :
      CHEAD(63) = ''                 !! 63 : MSIGN : A16   :   Data modified by               :
      CHEAD(64) = ''                 !! 64 : SIZE  : I10   :   Array size                     : 0
      WRITE(CHEAD(39), '(E15.7)') -999.0
      WRITE(CHEAD(40), '(E15.7)') -999.0
      WRITE(CHEAD(41), '(E15.7)') -999.0
      WRITE(CHEAD(42), '(E15.7)') -999.0
      WRITE(CHEAD(43), '(E15.7)') -999.0
      WRITE(CHEAD(44), '(I10)') 1
      WRITE(CHEAD(29), '(A16)') 'TX'
      WRITE(CHEAD(32), '(A16)') 'TY'
      WRITE(CHEAD(35), '(A16)') 'TZ'
      WRITE(CHEAD(30), '(I10)') 1
      WRITE(CHEAD(31), '(I10)') NX
      WRITE(CHEAD(33), '(I10)') 1
      WRITE(CHEAD(34), '(I10)') NY
      WRITE(CHEAD(36), '(I10)') 1
      WRITE(CHEAD(37), '(I10)') NZ
      WRITE(CHEAD(64), '(I10)') NX*NY*NZ
      CHEAD(27) = '09990101 000000 '
      CHEAD(50) = '  09990101000000'
      CHEAD(27) = YEAR//' 000000 '
      CHEAD(50) = '  '//YEAR//'000000'
      CHEAD(38) = 'UR4'

      OPEN(13, FILE=CFILE, FORM='UNFORMATTED', ACCESS='APPEND')
      WRITE(13) CHEAD
      WRITE(13) (((REAL(DAT3(I,J,K)),I=1,NX),J=1,NY),K=1,NZ)
      CLOSE(13)

  end subroutine write_gtoolfile

end module gtool_utilities
