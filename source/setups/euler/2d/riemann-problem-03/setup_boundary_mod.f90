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

    outer = inner

# if 0
    pres = pressure(inner)

    outer(:,DENS_VAR) = inner(:,DENS_VAR)
    outer(:,MOMX_VAR) = 0.0_dp ! inner(:,MOMX_VAR)
    outer(:,MOMY_VAR) = 0.0_dp ! inner(:,MOMY_VAR)
    outer(:,ENER_VAR) = pres/(kappa-1.0_dp)

    return

    coords = side_get_coords(side)
    center = sum(coords,dim=1) / REAL(N_NODES,dp)

    if ((side%boundary == NOR .and. center(2) < 0.5_dp) .or. (side%boundary == WES .and. center(1) < 0.5_dp)) then
        outer(:,DENS_VAR) = dens_UL
        outer(:,MOMX_VAR) = velx_UL * dens_UL
        outer(:,MOMY_VAR) = vely_UL * dens_UL
        outer(:,ENER_VAR) = pres_UL/(kappa-1.0_dp) + 0.5_dp*dens_UL*(velx_UL**2 + vely_UL**2)

    else if ((side%boundary == NOR .and. center(2) > 0.5_dp) .or. (side%boundary == EAS .and. center(1) < 0.5_dp)) then
        outer(:,DENS_VAR) = dens_UR
        outer(:,MOMX_VAR) = velx_UR * dens_UR
        outer(:,MOMY_VAR) = vely_UR * dens_UR
        outer(:,ENER_VAR) = pres_UR/(kappa-1.0_dp) + 0.5_dp*dens_UR*(velx_UR**2 + vely_UR**2)

    else if ((side%boundary == SOU .and. center(2) < 0.5_dp) .or. (side%boundary == WES .and. center(1) > 0.5_dp)) then
        outer(:,DENS_VAR) = dens_LL
        outer(:,MOMX_VAR) = velx_LL * dens_LL
        outer(:,MOMY_VAR) = vely_LL * dens_LL
        outer(:,ENER_VAR) = pres_LL/(kappa-1.0_dp) + 0.5_dp*dens_LL*(velx_LL**2 + vely_LL**2)

    else
        outer(:,DENS_VAR) = dens_LR
        outer(:,MOMX_VAR) = velx_LR * dens_LR
        outer(:,MOMY_VAR) = vely_LR * dens_LR
        outer(:,ENER_VAR) = pres_LR/(kappa-1.0_dp) + 0.5_dp*dens_LR*(velx_LR**2 + vely_LR**2)
    end if
# endif

end subroutine

end module
