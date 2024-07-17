!######################################################################
module module_digda_mpi
  use mpi
  implicit none

  ! Declaration of parameters and variables
  integer, parameter :: IROOT = 0
  integer :: IERR, NPROCS, MYRANK
  integer :: IREQ, ISTATUS
 !integer, parameter :: NCPU = 12
  integer, parameter :: NCPU = 16
 !integer, parameter :: NCPU = 24
  integer :: IRANK, ISIZE1, ISIZE2, ISIZE3, IIST, IIEN, IJKG
  integer :: ISEND1, ISEND2, ISEND3, ISEND4, ISEND5
  integer :: IRECV1, IRECV2, IRECV3, IRECV4, IRECV5
  integer :: IREQ1, IREQ2, IREQ3, IREQ4, IREQ5
  integer :: IREQ11, IREQ12, IREQ13, IREQ14, IREQ15
  integer :: IREQ21, IREQ22, IREQ23, IREQ24, IREQ25
  integer :: JJSTR(0:NCPU-1), JJEND(0:NCPU-1)

  ! Parallel digda_configuration
  logical, parameter :: use_mpi = .true.

contains

  subroutine initialize_digda_mpi()

    IF (use_mpi) THEN
      call MPI_INIT(IERR)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, NPROCS, IERR)
      call MPI_COMM_RANK(MPI_COMM_WORLD, MYRANK, IERR)

      IF ( NPROCS /= NCPU ) THEN
        IF ( MYRANK == 0 ) THEN
          PRINT *, '=== MPI ERROR ==='
          call MPI_FINALIZE(IERR)
          STOP
        END IF
      END IF

!     DO IRANK = 0, NPROCS-1
!       call set_parallel_range( &
!  &         1, NXYDIM, NPROCS, IRANK, &
!  &         JJSTR(IRANK), JJEND(IRANK) &
!  &    )
!     END DO
    END IF

  end subroutine initialize_digda_mpi

  subroutine set_parallel_range( &
 &  N1, N2, NPR, IRA, &
 &  IST, IEN )

    integer, intent(in) :: N1, N2
    integer, intent(in) :: NPR, IRA
    integer, intent(out) :: IST, IEN
    integer :: IWORK1, IWORK2

    IWORK1 = ( N2 - N1 + 1 ) / NPR
    IWORK2 = MOD( N2 - N1 + 1, NPR )
    IST = IRA*IWORK1 + N1 + MIN( IRA, IWORK2 )
    IEN = IST + IWORK1 - 1
    IF ( IWORK2 > IRA ) THEN
      IEN = IEN + 1
    END IF

  end subroutine set_parallel_range

end module module_digda_mpi
