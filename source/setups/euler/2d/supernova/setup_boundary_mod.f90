module setup_boundary_mod

use constants_mod

contains

subroutine setup_side_boundary(side,inner,outer)

    use setup_config_mod
    use equations_const_mod
    use mesh_types_mod, only: side_t

    type(side_t), intent(in)    :: side
    real(dp), intent(in)        :: inner(N_NODES,N_VARS)
    real(dp), intent(out)       :: outer(N_NODES,N_VARS)

    if (side%boundary == NOR) then
        outer(:,DENS_VAR) =  inner(:,DENS_VAR)
        outer(:,MOMX_VAR) = -inner(:,MOMX_VAR)
        outer(:,MOMY_VAR) =  inner(:,MOMY_VAR)
        outer(:,ENER_VAR) =  inner(:,ENER_VAR)

    else if (side%boundary == WES) then
        outer(:,DENS_VAR) =  inner(:,DENS_VAR)
        outer(:,MOMX_VAR) =  inner(:,MOMX_VAR)
        outer(:,MOMY_VAR) = -inner(:,MOMY_VAR)
        outer(:,ENER_VAR) =  inner(:,ENER_VAR)

    else
        outer(:,DENS_VAR) = dens0
        outer(:,MOMX_VAR) = dens0*velx0
        outer(:,MOMY_VAR) = dens0*vely0
        outer(:,ENER_VAR) = pres0/(kappa-1.0_dp) + 0.5_dp*dens0*(velx0**2 + vely0**2)

    end if

end subroutine

end module
