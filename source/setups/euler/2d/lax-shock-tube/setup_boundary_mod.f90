module setup_boundary_mod

use constants_mod

contains

subroutine setup_side_boundary(side,inner,outer)

    use setup_config_mod
    use equations_const_mod
    use mesh_types_mod, only: side_t
    use quad_utils_mod, only: side_get_coords
    use equations_mod, only: pressure

    type(side_t), intent(in)    :: side
    real(dp), intent(in)        :: inner(N_NODES,N_VARS)
    real(dp), intent(out)       :: outer(N_NODES,N_VARS)

    real(dp) :: coords(N_NODES,2), center(2)

    real(dp) :: dens(N_NODES)
    real(dp) :: velx(N_NODES)
    real(dp) :: vely(N_NODES)
    real(dp) :: pres(N_NODES)

    ! outer = inner

    select case(side%boundary)

        case (NOR)
            outer(:,DENS_VAR) = densL
            outer(:,MOMX_VAR) = densL*velxL
            outer(:,MOMY_VAR) = densL*velyL
            outer(:,ENER_VAR) = presL/(kappa-1.0_dp) + 0.5_dp*densL*(velxL**2+velyL**2)

        case (SOU)
            outer(:,DENS_VAR) = densR
            outer(:,MOMX_VAR) = densR*velxR
            outer(:,MOMY_VAR) = densR*velyR
            outer(:,ENER_VAR) = presR/(kappa-1.0_dp) + 0.5_dp*densR*(velxR**2+velyR**2)

        case default
            outer(:,DENS_VAR) =  inner(:,DENS_VAR)
            outer(:,MOMX_VAR) =  inner(:,MOMX_VAR)
            outer(:,MOMY_VAR) = -inner(:,MOMY_VAR)
            outer(:,ENER_VAR) =  inner(:,ENER_VAR)

    end select

end subroutine

end module
