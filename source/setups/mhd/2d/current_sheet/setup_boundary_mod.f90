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

    !! Dirichlet as "outflow" condition

    outer(:,DENS_VAR) = dens0

    outer(:,MOMX_VAR) = dens0*velx0
    outer(:,MOMY_VAR) = dens0*vely0
    outer(:,MOMZ_VAR) = dens0*velz0

    outer(:,ENER_VAR) = pres0/(kappa-1.0_dp) &
        + 0.5_dp*dens0*(velx0**2+vely0**2+velz0**2) &
        + 0.5_dp/mu0*(magx0**2+magy0**2+magz0**2) &
        + 0.5_dp/mu0*glmp0**2

    outer(:,MAGX_VAR) = magx0
    outer(:,MAGY_VAR) = magy0
    outer(:,MAGZ_VAR) = magz0

    outer(:,GLMP_VAR) = glmp0

end subroutine

end module
