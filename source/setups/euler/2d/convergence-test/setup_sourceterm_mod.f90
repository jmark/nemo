!# define MODULE_SOURCETERM 1

module setup_sourceterm_mod

use constants_mod

contains

pure subroutine calc_sourceterm(quad,rhs)

    use equations_const_mod
    use mesh_types_mod, only: quad_t

    type(quad_t), intent(inout) :: quad
    real(dp), intent(out)       :: rhs(N_NODES,N_NODES,N_VARS)

    ! rhs = 0.0_dp
    rhs = quad%aux%manufacsol%st

end subroutine

function setup_sourceterm_dt(quad) result(dt)

    use mesh_types_mod, only: quad_t

    type(quad_t), intent(in)    :: quad
    real(dp)                    :: dt

    dt = HUGE(1.0_dp)

end function

end module
