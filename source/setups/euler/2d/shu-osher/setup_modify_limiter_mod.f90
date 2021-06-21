!# define PP_SETUP_MODIFY_LIMITER 1

module setup_modify_limiter_mod

use constants_mod

contains

subroutine setup_modify_limiter_cb(quad)

    use setup_config_mod
    use equations_const_mod

    use mesh_types_mod, only: quad_t
    use quad_utils_mod, only: quad_get_center

    type(quad_t), intent(inout) :: quad

    real(dp) :: center(N_DIMS)

    return

    !! center = quad_get_center(quad)

    !! if (&
    !!     (abs(xlim(1)-center(1)) < quad%delta(1)) .or. &
    !!     (abs(xlim(2)-center(1)) < quad%delta(1)) .or. &
    !!     (abs(ylim(1)-center(2)) < quad%delta(2)) .or. &
    !!     (abs(ylim(2)-center(2)) < quad%delta(2))) then

    !!     quad%aux%hydro%blend = 0.0_dp
    !!     quad%aux%hydro%subcell = .true.
    !! 
    !!     quad%aux%hydro%qmin = 1
    !!     quad%aux%hydro%qmax = 0
    !! end if

end subroutine

end module
