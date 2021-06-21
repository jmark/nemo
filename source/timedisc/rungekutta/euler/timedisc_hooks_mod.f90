module timedisc_hooks_mod

use constants_mod

contains

subroutine rkstep_cb(quad)

    use mesh_types_mod, only: quad_t
    use share_mod, only: rt => runtime

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

    quad%pld%hydro%state = state

end subroutine

subroutine hook_mid_rk_timestep_500

    use share_mod, only: mesh
    use mesh_utils_mod, only: mesh_iterate_quads

    call mesh_iterate_quads(mesh,rkstep_cb)

end subroutine

end module
