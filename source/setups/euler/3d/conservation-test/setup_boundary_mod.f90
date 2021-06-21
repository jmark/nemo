module setup_boundary_mod

use constants_mod

contains

subroutine setup_side_boundary(side,inner,outer)

    use setup_config_mod
    use equations_const_mod
    use mesh_types_mod, only: side_t
    use quad_utils_mod, only: side_get_coords

    type(side_t), intent(in)    :: side
    real(dp), intent(in)        :: inner(N_NODES,N_NODES,N_VARS)
    real(dp), intent(out)       :: outer(N_NODES,N_NODES,N_VARS)

    real(dp) :: pres(N_NODES,N_NODES)
    real(dp) :: ener(N_NODES,N_NODES)

    real(dp) :: velx(N_NODES,N_NODES)
    real(dp) :: vely(N_NODES,N_NODES)
    real(dp) :: velz(N_NODES,N_NODES)

    real(dp) :: coords(N_NODES,N_NODES,N_DIMS)

# if 0
    outer = inner

    !return

    !! coords = side_get_coords(side)

    !! velx = angl0 * (-coords(:,2))
    !! vely = angl0 *   coords(:,1)

    !! velx = inner(:,MOMX_VAR)/inner(:,DENS_VAR)
    !! vely = inner(:,MOMY_VAR)/inner(:,DENS_VAR)
    !! velz = inner(:,MOMZ_VAR)/inner(:,DENS_VAR)

    !! outer(:,DENS_VAR) = dens0
    !! outer(:,MOMX_VAR) = dens0*velx
    !! outer(:,MOMY_VAR) = dens0*vely
    !! outer(:,ENER_VAR) = 0.5*dens0*(velx**2 + vely**2) + K_isothermal*dens0/(kappa-1)

    !! outer(:,ENER_VAR) = 0.5*dens0*(velx**2 + vely**2) + dens0*R_specific*temp_min*(1 + (dens0/dens_critical)**(kappa-1))
    !! outer(:,ENER_VAR) = 0.5*dens0*(velx**2 + vely**2) + pres0/(kappa-1)

# else

    !! reflecting wall
    select case (side%direction)

        case (X_DIR)

            outer(:,:,DENS_VAR) =  inner(:,:,DENS_VAR)
            outer(:,:,MOMX_VAR) = -inner(:,:,MOMX_VAR)
            outer(:,:,MOMY_VAR) =  inner(:,:,MOMY_VAR)
            outer(:,:,MOMZ_VAR) =  inner(:,:,MOMZ_VAR)
            outer(:,:,ENER_VAR) =  inner(:,:,ENER_VAR)

        case (Y_DIR)

            outer(:,:,DENS_VAR) =  inner(:,:,DENS_VAR)
            outer(:,:,MOMX_VAR) =  inner(:,:,MOMX_VAR)
            outer(:,:,MOMY_VAR) = -inner(:,:,MOMY_VAR)
            outer(:,:,MOMZ_VAR) =  inner(:,:,MOMZ_VAR)
            outer(:,:,ENER_VAR) =  inner(:,:,ENER_VAR)

        case (Z_DIR)

            outer(:,:,DENS_VAR) =  inner(:,:,DENS_VAR)
            outer(:,:,MOMX_VAR) =  inner(:,:,MOMX_VAR)
            outer(:,:,MOMY_VAR) =  inner(:,:,MOMY_VAR)
            outer(:,:,MOMZ_VAR) = -inner(:,:,MOMZ_VAR)
            outer(:,:,ENER_VAR) =  inner(:,:,ENER_VAR)

    end select
# endif

end subroutine

end module
