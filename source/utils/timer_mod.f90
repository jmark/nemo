module timer_mod

use constants_mod

type :: timer_t

    real(dp) :: tic,toc

end type

contains

function now()

    use mpi, only: mpi_wtime

    real(dp) :: now

    now = mpi_wtime()

end function

subroutine timer_reset(timer)

    type(timer_t), intent(inout) :: timer

    timer%tic = 0.0_dp
    timer%toc = 0.0_dp

end subroutine

subroutine timer_start(timer)

    type(timer_t), intent(inout) :: timer

    timer%tic = now()

end subroutine

subroutine timer_stop(timer)

    type(timer_t), intent(inout) :: timer

    timer%toc = now()

end subroutine

pure function timer_read(timer) result(timespan)

    type(timer_t), intent(in) :: timer
    real(dp) :: timespan
   
    timespan = timer%toc - timer%tic

end function

function timer_elapsed(timer) result(elapsed)

    type(timer_t), intent(in) :: timer
    real(dp) :: elapsed
   
    elapsed = now() - timer%tic

end function

end module
