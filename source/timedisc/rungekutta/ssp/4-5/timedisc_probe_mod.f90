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
              r0 => quad%aux%timedisc%register0, &
              r1 => quad%aux%timedisc%register1)

    select case(rt%rkstage)

        case(1)
            state =                uu + A3 * dt * flux

        case(2)
            state = B1 * r0 + B2 * uu + B3 * dt * flux

        case(3)
            state = C1 * r0 + C2 * uu + C3 * dt * flux

        case(4)
            state = D1 * r0 + D2 * uu + D3 * dt * flux

        case(5)
            state = E1 * r0 + E2 * uu + E3 * dt * flux + r1

    end select
    end associate

end subroutine

end module
