module timedisc_probe_mod

use constants_mod

contains

pure subroutine timedisc_probe(quad,flux,state)

    use mesh_types_mod, only: quad_t
    use share_mod, only: rt => runtime

    use timedisc_const_mod, only: a_rk
    use timedisc_const_mod, only: c_rk

    type(quad_t), intent(in)    :: quad

# if PP_N_DIMS == 2
    real(dp), intent(in)        :: flux(N_NODES,N_NODES,N_VARS)
    real(dp), intent(out)       :: state(N_NODES,N_NODES,N_VARS)
# elif PP_N_DIMS == 3
    real(dp), intent(in)        :: flux(N_NODES,N_NODES,N_NODES,N_VARS)
    real(dp), intent(out)       :: state(N_NODES,N_NODES,N_NODES,N_VARS)
# endif

# if PP_N_DIMS == 2
    real(dp) :: register(N_NODES,N_NODES,N_VARS)
# elif PP_N_DIMS == 3
    real(dp) :: register(N_NODES,N_NODES,N_NODES,N_VARS)
# endif

    if (rt%rkstage == 1) then
        register = flux
    else
        register = a_rk(rt%rkstage) * quad%aux%timedisc%register0 + flux
    end if

    state = quad%pld%hydro%state + c_rk(rt%rkstage) * rt%timestep * register

end subroutine

end module
