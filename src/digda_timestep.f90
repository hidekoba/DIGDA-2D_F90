!######################################################################
module module_digda_timestep
  use module_zcfile
  implicit none

  ! Time step settings
!!! 1 year
!!integer :: iutstp = 1
!!double precision :: sdtmstp = 1

!!! 1/2 year
!!integer :: iutstp = 2
!!double precision :: sdtmstp = 6

!!! 50 day
!!integer :: iutstp = 3
!!double precision :: sdtmstp = 50

!!! 30 day
!!integer :: iutstp = 3
!!double precision :: sdtmstp = 30

!!! 10 day
!!integer :: iutstp = 3
!!double precision :: sdtmstp = 10

!!! 5 day
!!integer :: iutstp = 3
!!double precision :: sdtmstp = 5

!!! 2 day
!!integer :: iutstp = 3
!!double precision :: sdtmstp = 2

  ! 1 day
  double precision :: sdtmstp = 1
  integer :: iutstp = 3

!!! 1/10 day
!!double precision :: sdtmstp = 0.1
!!integer :: iutstp = 3

!!! 1/25 day
!!double precision :: sdtmstp = 0.04
!!integer :: iutstp = 3

!!! 1/50 day
!!double precision :: sdtmstp = 0.02
!!integer :: iutstp = 3

  double precision :: tssec !! [s]
  double precision :: sdtsyr !! [yr]
  character(len=6), parameter :: calunit(6) = (/ 'year  ', 'month ', 'day   ', 'hour  ', 'minute', 'sec   ' /)

  ! Conversion parameter
  integer :: imonyr = 12
  integer :: idaymo = 30
  integer :: ihrday = 24
  integer :: iminhr = 60
  integer :: isecmn = 60
  integer :: isecyr

  ! Start and end year
  integer :: tstr = 1
  integer :: sdyrstr = 0
  integer :: sdyrend = 500
  integer :: tend
  logical :: wotstp = .true.  !! write out time step ?
  namelist /nmsdtm/ sdtmstp, iutstp, sdyrstr, sdyrend, wotstp

  ! Interval for outputting history files
  integer :: sdyrhst = 100
  integer :: thst
  logical :: wohist = .true.  !! write out history file ?
  logical :: hstavg = .true.  !! output the average file as history file ?
  namelist /nmsdhi/ sdyrhst, wohist, hstavg

  contains

  subroutine initialize_timestep()

    integer :: ios

    ! Open and read the namelist file
    rewind(nfparm)
    read(nfparm, NML=nmsdtm, IOSTAT=ios)
    write(nfstdo, nmsdtm)
    rewind(nfparm)
    read(nfparm, NML=nmsdhi, IOSTAT=ios)
    write(nfstdo, nmsdhi)

    ! Calculate derived values
    isecyr = isecmn * iminhr * ihrday * idaymo * imonyr
    if (iutstp == 6) then
      tssec = sdtmstp
    else if (iutstp == 5) then
      tssec = sdtmstp * isecmn
    else if (iutstp == 4) then
      tssec = sdtmstp * isecmn * iminhr
    else if (iutstp == 3) then
      tssec = sdtmstp * isecmn * iminhr * ihrday
    else if (iutstp == 2) then
      tssec = sdtmstp * isecmn * iminhr * ihrday * idaymo
    else if (iutstp == 1) then
      tssec = sdtmstp * isecmn * iminhr * ihrday * idaymo * imonyr
    end if
    tend = (sdyrend - sdyrstr) * isecyr / tssec !! Converts years to time steps
    thst = sdyrhst * isecyr / tssec !! Converts years to time steps
    sdtsyr = tssec / isecmn / iminhr / ihrday / idaymo / imonyr !! Converts seconds to years

    write(nfstdo, *) "'--------------------------------------------------'"
    write(nfstdo, *) "'-- Time Control Parameters                      --'"
    write(nfstdo, *)
    write(nfstdo, *) '        Start year :', sdyrstr
    write(nfstdo, *) '          End year :', sdyrend
    write(nfstdo, *) '      History year :', sdyrhst
    write(nfstdo, *) '         Time step :', sdtmstp, ' ', calunit(iutstp)
    write(nfstdo, *) '                   :', tssec, ' sec'
    write(nfstdo, *) '     End time step :', tend, ' steps'
    write(nfstdo, *) ' History time step :', thst, ' steps'
    write(nfstdo, *) "'--------------------------------------------------'"

  end subroutine initialize_timestep

end module module_digda_timestep
