!# define DO_PROFILING 0

module profiler_mod

use constants_mod
use timer_mod

character(len=*), parameter :: reportfile = 'profiling.dat'

integer, parameter :: nslots = 32

integer, save :: nactive = 0

type :: profiler_t

    integer             :: id
    character(len=256)  :: name
    type(timer_t)       :: timer

    real(dp)            :: total = 0
    integer             :: count = 0

end type

type(profiler_t), save :: slots(nslots)

contains

function profiler_init(name) result(id)
    
    use status_mod, only: raise_error

    character(len=*), intent(in)    :: name
    integer                         :: id

    if (nactive >= nslots) then
        call raise_error('Increase number of available profiler slots!')
        id = 0
        return
    end if

    nactive = nactive + 1

    id = nactive

    slots(id)%name = trim(name)
    slots(id)%id = id

    slots(id)%total = 0.0_dp
    slots(id)%count = 0

end function

subroutine profiler_start(id)

    integer, intent(in) :: id 

    call timer_start(slots(id)%timer)

end subroutine

subroutine profiler_stop(id)

    integer, intent(in) :: id 

    call timer_stop(slots(id)%timer)

    slots(id)%total = slots(id)%total + timer_read(slots(id)%timer)
    slots(id)%count = slots(id)%count + 1

end subroutine

subroutine profiler_report
    
    use share_mod, only: rt => runtime
    use utils_mod, only: utils_get_newunit

    integer :: id
    integer :: ioUnit
    integer :: openStat

    if (rt%mpi%root) then
        OPEN(UNIT     = utils_get_newunit(ioUnit), &
             FILE     = TRIM(reportfile)   , &
             FORM     = 'FORMATTED'        , &
             STATUS   = 'UNKNOWN'          , &
             POSITION = 'APPEND'           , &
             RECL     = 50000              , &
             IOSTAT   = openStat             )

        IF (openStat.NE.0) THEN
            write (*,*) 'ERROR: cannot open '// TRIM(reportfile)
            return
        END IF

        do id = 1,nactive
            write(ioUnit, '((A32),(I9),(ES24.16),(I9),(ES24.16))') slots(id)%name, rt%timesteps, rt%simtime, slots(id)%count, slots(id)%total
        end do

        close(ioUnit)
    end if

end subroutine

end module
