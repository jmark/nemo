module timedisc_probe_mod

use constants_mod

contains

pure subroutine timedisc_probe(quad,flux,state)

    use mesh_types_mod, only: quad_t
    use share_mod, only: rt => runtime

    use timedisc_const_mod

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

    select case(rt%rkstage)

        case(1)
            state = uu + dt * flux

        case(2)
            state = 3.0_dp/4.0_dp * r0 + 1.0_dp/4.0_dp * uu + 1.0_dp/4.0_dp * dt * flux

        case(3)
            state = 1.0_dp/3.0_dp * r0 + 2.0_dp/3.0_dp * uu + 2.0_dp/3.0_dp * dt * flux

    end select
    end associate

end subroutine

end module
