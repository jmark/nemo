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
    real(dp), intent(in)        :: flux(N_NODES,N_NODES,N_NODES,N_VARS)
    real(dp), intent(out)       :: state(N_NODES,N_NODES,N_NODES,N_VARS)
# endif

    associate(dt => rt%timestep, &
              uu => quad%pld%hydro%state, &
              r0 => quad%aux%timedisc%register)

    if (rt%rkstage == 1) then
        state = uu + 2.0/3.0_dp*dt*flux
        !state = quad%pld%hydro%state + 2.0/3.0_dp * rt%timestep * flux
    else
        !state = quad%pld%hydro%state + rt%timestep * (-5.0_dp/12.0_dp * quad%aux%timedisc%register + 0.75_dp * flux)
        state = 5.0/8.0*r0 + 3.0/8.0*uu + 3.0/4.0*dt*flux
    end if
    end associate

end subroutine

end module
