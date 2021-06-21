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
    use equations_mod, only: cons2entr

    use timedisc_const_mod

    use kernel_evolve_mod, only: kernel_evolve

    type(quad_t), intent(inout) :: quad

    integer :: i,j

# if PP_N_DIMS == 2
    real(dp) :: flux(N_NODES,N_NODES,N_VARS)
    real(dp) :: state(N_NODES,N_NODES,N_VARS)
# elif PP_N_DIMS == 3
    real(dp) :: flux(N_NODES,N_NODES,N_NODES,N_VARS)
    real(dp) :: state(N_NODES,N_NODES,N_NODES,N_VARS)
# endif

    call kernel_evolve(quad,flux,state)

    if (rt%rkstage == 1) then
        quad%aux%timedisc%register = quad%pld%hydro%state
    end if

    ! quad%aux%hydro%entropy = quad%aux%hydro%entropy + (cons2entr(state) - cons2entr(quad%pld%hydro%state))/rt%timestep
    ! quad%aux%hydro%entropy = quad%aux%hydro%entropy + (cons2entr(state) - cons2entr(quad%pld%hydro%state))/rt%timestep

    quad%pld%hydro%state = state

end subroutine

subroutine hook_mid_rk_timestep_500

    use share_mod, only: mesh
    use mesh_utils_mod, only: mesh_iterate_quads

    call mesh_iterate_quads(mesh,rkstep_cb)

end subroutine

end module
