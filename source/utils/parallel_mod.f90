module parallel_mod

use constants_mod
use mpi_f08

character(MPI_MAX_PROCESSOR_NAME) hostname

# if defined(__GFORTRAN__)
interface
    function MPI_Comm_f2c(fcomm) result(ccomm) bind(c,name="MPI_Comm_f2c")
        use iso_c_binding, only: c_ptr, c_int

        integer(c_int), intent(in), value   :: fcomm
        type(c_ptr)                         :: ccomm

    end function
end interface
# endif

interface
    !! int sched_getcpu(void);
    function sched_getcpu() result(cpuid) bind(c,name="sched_getcpu")
        use iso_c_binding, only: c_int
        integer(c_int) :: cpuid
    end function
end interface

type(MPI_COMM), save :: COMM

integer :: nranks,ierr

contains

subroutine hook_init_100

    use share_mod, only: rt => runtime

    integer :: ierr,provided

# ifdef _OPENMP
    call MPI_Init_thread(MPI_THREAD_FUNNELED, provided, ierr)
# else
    call MPI_Init_thread(MPI_THREAD_SINGLE, provided, ierr)
# endif

    rt%mpi%comm = MPI_COMM_WORLD

    call MPI_COMM_RANK(rt%mpi%comm, rt%mpi%rank, ierr)
    call MPI_COMM_SIZE(rt%mpi%comm, rt%mpi%size, ierr)

    rt%mpi%root = rt%mpi%rank == 0
    rt%mpi%info = MPI_INFO_NULL
    rt%mpi%ccomm = MPI_Comm_f2c(rt%mpi%comm%mpi_val)

    rt%ncpus = rt%mpi%size

    COMM = rt%mpi%comm
    nranks = rt%mpi%size

end subroutine

subroutine hook_nuke_900

    call MPI_FINALIZE(ierr)

end subroutine

subroutine parallel_ping(dest,tag)

    integer, intent(in)             :: dest
    integer, intent(in), optional   :: tag
    integer                         :: ierr

    if (present(tag)) then
        call MPI_SEND(0, 1, MPI_INTEGER, dest, tag, COMM, ierr)
    else
        call MPI_SEND(0, 1, MPI_INTEGER, dest, 0, COMM, ierr)
    endif

end subroutine

subroutine parallel_pong(src,tag)

    integer, intent(in)             :: src
    integer, intent(in), optional   :: tag

    integer                         :: ierr
    type(MPI_STATUS)                :: stat

    if (present(tag)) then
        call MPI_Probe(src, tag, COMM, stat, ierr)
    else
        call MPI_Probe(src, 0, COMM, stat, ierr)
    end if

end subroutine

subroutine barrier(ierr)

    integer, intent(out), optional :: ierr

    integer :: ierr2

    if (nranks < 2) return

    if (present(ierr)) then
         call MPI_BARRIER(comm, ierr)
    else
         call MPI_BARRIER(comm, ierr2)
    end if

end subroutine

end module
