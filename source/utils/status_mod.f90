module status_mod

!use status_types_mod

contains

subroutine raise_error(msg)

    use share_mod, only: rt => runtime

    character(len=*), intent(in) :: msg

    !$omp critical
    rt%stat%ok = .false.
    rt%stat%errcode = -1
    rt%stat%message = trim(rt%stat%message) // trim(msg) // NEW_LINE('A')

    !write (*,*) rt%stat%message
    !$omp end critical

end subroutine

function sighandler(erno) result(retval)

    use share_mod, only: rt => runtime

    integer, intent(in) :: erno
    integer :: retval

    write (*,*) 'got signal', rt%mpi%rank, erno

    retval = 0

end function

end module
