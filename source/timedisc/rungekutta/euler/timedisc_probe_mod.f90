module timedisc_probe_mod

use constants_mod

contains

pure subroutine timedisc_probe(quad,flux,state)

    use mesh_types_mod, only: quad_t
    use share_mod, only: rt => runtime

    type(quad_t), intent(in)    :: quad

# if PP_N_DIMS == 2
    real(dp), intent(in)        :: flux(N_NODES,N_NODES,N_VARS)
    real(dp), intent(out)       :: state(N_NODES,N_NODES,N_VARS)
# elif PP_N_DIMS == 3
    real(dp), intent(in)        :: flux(N_NODES,N_NODES,N_VARS)
    real(dp), intent(out)       :: state(N_NODES,N_NODES,N_VARS)
# endif

    state = quad%pld%hydro%state + rt%timestep * flux

end subroutine

end module
