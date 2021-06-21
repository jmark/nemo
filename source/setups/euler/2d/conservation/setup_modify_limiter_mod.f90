!# define PP_SETUP_MODIFY_LIMITER 1

module setup_modify_limiter_mod

use constants_mod

contains

subroutine setup_modify_limiter_cb(quad)

    use setup_config_mod
    use equations_const_mod
    use kernel_const_mod

    use mesh_types_mod, only: quad_t

    type(quad_t), intent(inout) :: quad

    integer :: i

# if 0
    quad%aux%hydro%subcell = .false.

    quad%aux%hydro%qmin = 3
    quad%aux%hydro%qmax = 3

    quad%aux%hydro%blend = 1.0_dp

    return
# endif

# if PP_KERNEL_DGFV
    call RANDOM_NUMBER(quad%aux%hydro%blend)

    quad%aux%hydro%subcell = .true.

    quad%aux%hydro%qmin = 1
    quad%aux%hydro%qmax = 3

    ! quad%aux%hydro%blend = 1.0_dp
    ! quad%aux%hydro%subcell = .false.
    ! quad%aux%hydro%qmin = DGFV_QMAX
    ! quad%aux%hydro%qmax = DGFV_QMAX


    ! quad%aux%hydro%blend = 1.0_dp
    ! quad%aux%hydro%subcell = .true.
    ! quad%aux%hydro%subcell = .false.
    ! quad%aux%hydro%qmin = 2
    ! quad%aux%hydro%qmax = 2
# endif

end subroutine

end module
