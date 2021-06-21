# define TRY(CALL) ok = CALL
# define TERMINATE 715 continue

program nemo_prog

use share_mod, only: traceback
use share_mod, only: mesh, rt => runtime

use parallel_mod, only: barrier

use hook_init_mod, only: hook_init
use hook_nuke_mod, only: hook_nuke

use parallel_mod, only: sched_getcpu
use parallel_mod, only: parallel_ping
use parallel_mod, only: parallel_pong

use driver_mod, only: mainloop

# if PP_CHECKPOINT
use hook_checkpoint_mod, only: hook_checkpoint
# endif

use timer_mod

use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
                                          stdout=>output_unit, &
                                          stderr=>error_unit

# ifdef _OPENMP
USE OMP_LIB, only: omp_get_thread_num
USE OMP_LIB, only: omp_get_num_threads
# endif

logical :: ok
integer :: ierr

real(dp) :: elapsed
type(timer_t) :: timer

ok = .true.

!! ========================================================================= !!
!! Define defaults.

call hook_init()

if (rt%mpi%root) then
    write (*,*) ''
    write (*,*) '  _   _ ______ __  __  ____             .                 '
    write (*,*) ' | \ | |  ____|  \/  |/ __ \           ":"                '
    write (*,*) ' |  \| | |__  | \  / | |  | |        ___:____     |"\/"|  '
    write (*,*) " | . ` |  __| | |\/| | |  | |      ,'        `.    \  /   "
    write (*,*) ' | |\  | |____| |  | | |__| |      |  O        \___/  |   '
    write (*,*) ' |_| \_|______|_|  |_|\____/    -~^~^~^~^~^~^~^~^~^~^~^~^~'
    write (*,*) '                                  -~^~^ ~^~ ^~^~^ ~^~^-   '
    write (*,*) '                                                          '      
    write (*,*) '         ## v0.9 ##             Mobilis adversus flumine. '
    write (*,*) ''
    !write (*,*) '## compiled at:  ' // ${time.strftime("'%Y-%m-%d %H:%M:%S'")}$
    write (*,*)
    
    if (rt%mpi%size > 1) then
        write (*,*) '## MPI:    number of ranks  ', rt%mpi%size
    end if

# ifdef _OPENMP
    !$OMP Parallel default(shared)
    if (omp_get_thread_num() == 0) then
        write (*,*) '## OpenMP: number of threads', omp_get_num_threads()
    end if
    !$OMP end Parallel
# endif
    write (*,*)
endif

call flush(stdout)
call barrier()

# ifdef _OPENMP
! if (rt%mpi%rank > 0) call parallel_pong(rt%mpi%rank-1)

!$OMP Parallel default(shared)
if (omp_get_thread_num() == 0) then
    rt%ncpus = rt%ncpus * omp_get_num_threads()
end if
!write (*,*) 'mpirank / threadid / cpuid', rt%mpi%rank, omp_get_thread_num(), sched_getcpu()
!$OMP end Parallel

!call flush(stdout)
!call barrier()

! if (rt%mpi%rank + 1 < rt%mpi%size) call parallel_ping(rt%mpi%rank+1)
! call barrier()

if (rt%mpi%root) write (*,*)
# endif

!! ========================================================================= !!

if (rt%mpi%root) then
    call timer_reset(timer)
    call timer_start(timer)
end if

if (rt%stat%ok) then
    !! ========================================================================= !!
    !! Load initial state from checkpoint if given.
    !! TODO: Not implemented yet.

    !if (len_trim(rt%chkptfile) > 0) then
    !    if (rt%mpi%root%world) write (*,*) "Loading checkpoint file " // trim(rt%chkptfile)
    !    !TRY(load_checkpoint_file(runtime))
    !end if

    call mainloop()
    ! call hook_checkpoint()
end if

if (rt%mpi%root) then
    call timer_stop(timer)
    elapsed = timer_read(timer)
end if

if (.not.rt%stat%ok) then
    if (rt%mpi%rank > 0) then
        call parallel_pong(rt%mpi%rank-1)
    end if

    if (.not. rt%stat%ok) then
        write (*,*)
        write (*,11) ' ERROR in rank: ', rt%mpi%rank
        write (*,*) 'errcode:', rt%stat%errcode
        write (*,*) 'message: ', trim(rt%stat%message)
        write (*,*)
    end if

    !if (.not. traceback%ok) then
    !    write (*,*)
    !    write (*,11) ' ERROR in rank: ', rt%mpi%rank
    !    write (*,*) 'errcode:', traceback%errcode
    !    write (*,*) 'message: ', trim(traceback%message)
    !    write (*,*)
    !end if

    if (rt%mpi%rank + 1 < rt%mpi%size) then
        call parallel_ping(rt%mpi%rank+1)
    end if

    if (rt%mpi%root) then
        write (*,*)
        write (*,2) '[FAILURE]'
    end if
else
    if (rt%mpi%root) write(*,*)
    if (rt%mpi%root) write(*,2) '[SUCCESS] total runtime:', elapsed, 's'
end if

call hook_nuke()

2 FORMAT(a25,5x,F24.5,a4)
11 FORMAT(a15,1x,i5)
1001 format ('[CHECKPOINT]', i5, 4x, i6, 2x, ES12.3, 2x, ES12.4,2x,ES12.4)

end program
