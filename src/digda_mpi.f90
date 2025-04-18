!######################################################################
module module_digda_mpi
  use mpi
  implicit none
#include "mpif.h"

  ! Declaration of parameters and variables
  integer, parameter :: iroot = 0
  integer :: ierr, nprocs, myrank
  integer :: ireq, istatus
 !integer, parameter :: ncpu = 12
 !integer, parameter :: ncpu = 16
 !integer, parameter :: ncpu = 24
 !integer, parameter :: ncpu = 45
  integer, parameter :: ncpu = 48 !! 256x192
  integer, parameter :: ncpu = 64
 !integer, parameter :: ncpu = 90 !! 2160x1080
 !integer, parameter :: ncpu = 96
 !integer, parameter :: ncpu = 180
 !integer, parameter :: ncpu = 192
  integer :: nsndct2d, nsndct3d, nsndcw3d
  integer :: irank, isize1, isize2, isize3, ijkg
  integer :: isend1, isend2, isend3
  integer :: irecv1, irecv2, irecv3
  integer :: jjstr(0:ncpu-1), jjend(0:ncpu-1)

  ! Parallel digda_configuration
  logical, parameter :: use_mpi = .true.

contains

  subroutine initialize_digda_mpi()

    call mpi_init(ierr)
    call mpi_comm_size(mpi_comm_world, nprocs, ierr)
    call mpi_comm_rank(mpi_comm_world, myrank, ierr)

    if ( nprocs /= ncpu ) then
      if ( myrank == iroot ) then
        write(*, *) "***** MPI ERROR *****"
        call mpi_finalize(ierr)
        stop
      end if
    end if

  end subroutine initialize_digda_mpi

  subroutine set_parallel_range( n1, n2, npr, ira, & !! IN
                               & ist, ien          ) !! OUT

    integer, intent(in)  :: n1, n2
    integer, intent(in)  :: npr, ira
    integer, intent(out) :: ist, ien
    integer :: iwork1, iwork2

    iwork1 = ( n2 - n1 + 1 ) / npr
    iwork2 = mod( n2 - n1 + 1, npr )
    ist = ira * iwork1 + n1 + min( ira, iwork2 )
    ien = ist + iwork1 - 1
    if ( iwork2 > ira ) then
      ien = ien + 1
    end if

  end subroutine set_parallel_range

  subroutine finalize_digda_mpi()

    call mpi_finalize(ierr)

  end subroutine finalize_digda_mpi

end module module_digda_mpi
