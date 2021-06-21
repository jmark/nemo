!! ========================================================================== !!
!!
!!   5 stages / 4th - order / SSP / Runge-Kutta Scheme
!!
!! ========================================================================== !!
!!
!! Raymond J. Spiteri and Steven J. Ruuth
!! A NEW CLASS OF OPTIMAL HIGH-ORDER STRONG-STABILITY-PRESERVING TIME DISCRETIZATION METHODS
!! 2002
!! SIAM J. NUMER. ANAL.
!! Vol. 40, No. 2, pp. 469-491
!!
!! ========================================================================== !!

module timedisc_hooks_mod

use constants_mod

contains

subroutine rkstep_cb(quad)

    use mesh_types_mod, only: quad_t
    use share_mod, only: rt => runtime

    use timedisc_const_mod

    use kernel_evolve_mod, only: kernel_evolve

    type(quad_t), intent(inout) :: quad

# if PP_N_DIMS == 2
    real(dp) :: flux(N_NODES,N_NODES,N_VARS)
    real(dp) :: state(N_NODES,N_NODES,N_VARS)
# elif PP_N_DIMS == 3
    real(dp) :: flux(N_NODES,N_NODES,N_NODES,N_VARS)
    real(dp) :: state(N_NODES,N_NODES,N_NODES,N_VARS)
# endif

    call kernel_evolve(quad,flux,state)

    associate(dt => rt%timestep, &
              uu => quad%pld%hydro%state, &
              r0 => quad%aux%timedisc%register0, &
              r1 => quad%aux%timedisc%register1)

    select case(rt%rkstage)

        case(1)
            r0 = uu

        case(2)
            r1 = E4 * state

        case(3)
            r1 = r1 + E5 * state

        case(4)
            r1 = r1 + E6 * dt * flux

        case(5)
            !! Do nothing.

    end select
    end associate

    quad%pld%hydro%state = state

end subroutine

subroutine hook_mid_rk_timestep_500

    use share_mod, only: mesh
    use mesh_utils_mod, only: mesh_iterate_quads

    call mesh_iterate_quads(mesh,rkstep_cb)

end subroutine

end module
