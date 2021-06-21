module buffer_mod

use mpi_f08, only: MPI_Win
use mpi_f08, only: MPI_Comm
use mpi_f08, only: MPI_ADDRESS_KIND
use iso_c_binding, only: c_ptr

type :: buffer_shared_t

    logical                         :: active = .false.

    integer(kind=MPI_ADDRESS_KIND)  :: size !! in bytes
    integer(kind=8)                 :: elem_size
    integer(kind=8)                 :: elem_count

    integer                         :: mpirank
    type(MPI_Comm)                  :: mpicomm
    type(MPI_Win)                   :: window
    type(c_ptr)                     :: baseptr

    integer                         :: disp_unit

end type

type :: buffer_serial_t

    logical                         :: active = .false.

    integer(kind=8)                 :: size !! in bytes
    integer(kind=8)                 :: elem_size
    integer(kind=8)                 :: elem_count

    character(len=1), allocatable   :: buffer(:)
    type(c_ptr)                     :: baseptr

end type

interface buffer_allocate
    procedure :: buffer_serial_allocate
    procedure :: buffer_shared_allocate
end interface

interface buffer_free
    procedure :: buffer_serial_free
    procedure :: buffer_shared_free
end interface

interface buffer_move
    procedure :: buffer_serial_move
    procedure :: buffer_shared_move
end interface

contains

subroutine buffer_serial_allocate(elem_size,elem_count,handle)

    use iso_c_binding, only: c_loc
    
    integer(kind=8), intent(in)         :: elem_size !! in bytes
    integer(kind=8), intent(in)         :: elem_count
    type(buffer_serial_t), target, intent(out)  :: handle

    handle%elem_size = elem_size
    handle%elem_count = elem_count
    handle%size = elem_size * elem_count

    allocate(handle%buffer(handle%size))

    handle%baseptr = c_loc(handle%buffer)

    handle%active = .true.

end subroutine

subroutine buffer_serial_free(handle)
    
    use iso_c_binding, only: c_null_ptr

    type(buffer_serial_t), intent(inout)   :: handle

    if (allocated(handle%buffer)) then
        handle%active = .false.
        deallocate(handle%buffer)
        handle%baseptr = c_null_ptr

        handle%elem_size = 0
        handle%elem_count = 0
        handle%size = 0

    end if

end subroutine

subroutine buffer_serial_move(src,snk)

    use iso_c_binding, only: c_null_ptr, c_loc
    
    type(buffer_serial_t), target, intent(inout) :: src,snk

    snk%active = .true.
    snk%elem_size = src%elem_size
    snk%elem_count = src%elem_count
    snk%size = src%size

    call move_alloc(src%buffer,snk%buffer)
    snk%baseptr = c_loc(snk%buffer)

    !! -- !!

    src%active = .false.
    src%baseptr = c_null_ptr
    src%elem_size = 0
    src%elem_count = 0
    src%size = 0

end subroutine

subroutine buffer_shared_allocate(elem_size,elem_count,comm,handle)
    
    use MPI_f08, only: MPI_Comm
    use mpi_f08, only: MPI_Comm_rank
    use mpi_f08, only: MPI_INFO_NULL
    use mpi_f08, only: MPI_Address_Kind
    use mpi_f08, only: MPI_Win_allocate_shared

    integer(kind=8), intent(in)                 :: elem_size !! in bytes
    integer(kind=8), intent(in)                 :: elem_count !! in bytes
    type(MPI_Comm), intent(in)                  :: comm
    type(buffer_shared_t), intent(out)          :: handle

    integer :: ierror

    call MPI_Comm_rank(comm,handle%mpirank,ierror)

    handle%elem_size    = elem_size
    handle%elem_count   = elem_count
    handle%mpicomm      = comm

    if (handle%mpirank == 0) then
        handle%size = elem_size * elem_count
        call MPI_Win_allocate_shared(handle%size, 1, MPI_INFO_NULL, comm, handle%baseptr, handle%window, ierror)
    else
        handle%size = 0
        call MPI_Win_allocate_shared(handle%size, 1, MPI_INFO_NULL, comm, handle%baseptr, handle%window, ierror)
        call MPI_Win_shared_query(handle%window, 0, handle%size, handle%disp_unit, handle%baseptr, ierror)
    end if

    handle%active = .true.

end subroutine

subroutine buffer_shared_free(handle)
    
    use iso_c_binding, only: c_null_ptr

    type(buffer_shared_t), intent(inout) :: handle

    integer :: ierror, assert

    if (handle%active) then
        handle%active = .false.
        call MPI_Win_free(handle%window, ierror)
        handle%baseptr = c_null_ptr

        handle%elem_size = 0
        handle%elem_count = 0
        handle%size = 0
    end if

end subroutine

subroutine buffer_shared_move(src,snk)

    use iso_c_binding, only: c_null_ptr
    use mpi_f08, only: MPI_COMM_NULL, MPI_WIN_NULL
    
    type(buffer_shared_t), intent(inout) :: src,snk

    snk = src

    src%active = .false.
    src%window = MPI_WIN_NULL
    src%mpicomm = MPI_COMM_NULL
    src%baseptr = c_null_ptr
    src%disp_unit = -1

    src%elem_size = 0
    src%elem_count = 0
    src%size = 0

end subroutine

end module
