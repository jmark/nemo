!! ========================================================================== !!
!!
!! 4th - order / 5 stages / low-storage / Runge-Kutta Scheme
!! 
!! ========================================================================== !!
!!
!! Kennedy, Christopher A., Mark H. Carpenter, and R. Michael Lewis:
!! "Low-storage, explicit Runge–Kutta schemes for the compressible Navier–Stokes
!! equations." Applied numerical mathematics 35.3 (2000): 177-219.
!!
!! ========================================================================== !!

module timedisc_hooks_mod

use constants_mod

contains

subroutine rkstep_cb(quad)

    use mesh_types_mod, only: quad_t
    use share_mod, only: rt => runtime

    use timedisc_const_mod, only: a_rk
    use timedisc_const_mod, only: c_rk

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

    if (rt%rkstage == 1) then
        quad%aux%timedisc%register0 = flux
    else
        quad%aux%timedisc%register0 = a_rk(rt%rkstage) * quad%aux%timedisc%register0 + flux
    end if

    quad%pld%hydro%state = state

end subroutine

subroutine hook_mid_rk_timestep_500

    use share_mod, only: mesh
    use mesh_utils_mod, only: mesh_iterate_quads

    call mesh_iterate_quads(mesh,rkstep_cb)

end subroutine

end module
