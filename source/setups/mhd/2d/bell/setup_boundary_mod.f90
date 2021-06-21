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

# if 1
    !! Dirichlet

    outer(:,DENS_VAR) = dens0

    outer(:,MOMX_VAR) = dens0*velx0
    outer(:,MOMY_VAR) = dens0*vely0
    outer(:,MOMZ_VAR) = dens0*velz0

    outer(:,ENER_VAR) = ener0

    outer(:,MAGX_VAR) = magx0
    outer(:,MAGY_VAR) = magy0
    outer(:,MAGZ_VAR) = magz0

    outer(:,GLMP_VAR) = glmp0

# else

    !! reflecting wall
    select case (side%direction)

        case (X_DIR)

            outer(:,DENS_VAR) =  inner(:,DENS_VAR)
            outer(:,MOMX_VAR) = -inner(:,MOMX_VAR)
            outer(:,MOMY_VAR) =  inner(:,MOMY_VAR)
            outer(:,ENER_VAR) =  inner(:,ENER_VAR)

        case (Y_DIR)

            outer(:,DENS_VAR) =  inner(:,DENS_VAR)
            outer(:,MOMX_VAR) =  inner(:,MOMX_VAR)
            outer(:,MOMY_VAR) = -inner(:,MOMY_VAR)
            outer(:,ENER_VAR) =  inner(:,ENER_VAR)

    end select
# endif

end subroutine

end module
