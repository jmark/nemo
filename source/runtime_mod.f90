module runtime_mod

use constants_mod
use mpi_f08, only: MPI_COMM,MPI_INFO
use iso_c_binding, only: c_ptr, c_int

type :: runtime_mpi_t

    type(MPI_COMM) :: comm

# if defined(__INTEL_COMPILER)
    integer(c_int) :: ccomm
# else
    type(c_ptr) :: ccomm
# endif
    integer :: rank = 0
    integer :: size = 1

    integer :: erno = 0
    type(MPI_INFO) :: info

    logical :: root = .false.

end type

integer, parameter :: STATUS_OK     = 0
integer, parameter :: STATUS_ERROR  = 1
integer, parameter :: STARTUP_ERROR = -1

type status_t

    logical                 :: ok       = .true.
    integer                 :: errcode  = STATUS_OK
    character(len=4096)     :: message  = ''
    
end type

type :: runtime_t

    real(dp) :: simtime
    real(dp) :: timestep
    real(dp) :: maxtimestep

    real(dp) :: rktimestep
    real(dp) :: rksimtime
    integer  :: rkstage

    integer :: chkptcount  = 0
    integer :: timesteps   = 0

    integer :: ncpus = 1

    real(dp) :: nextdump

    type(runtime_mpi_t) :: mpi

    type(status_t) :: stat

end type

end module
