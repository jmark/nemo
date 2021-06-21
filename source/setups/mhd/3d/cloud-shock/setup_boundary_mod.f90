module setup_boundary_mod

use constants_mod

contains

subroutine setup_side_boundary(side,inner,outer)

    use setup_config_mod
    use equations_const_mod
    use mesh_types_mod, only: side_t
    use quad_utils_mod, only: side_get_coords

    use share_mod, only: rt => runtime

    type(side_t), intent(in)    :: side
    real(dp), intent(in)        :: inner(N_NODES,N_NODES,N_VARS)
    real(dp), intent(out)       :: outer(N_NODES,N_NODES,N_VARS)

    real(dp) :: dens(N_NODES,N_NODES)

    real(dp) :: velx(N_NODES,N_NODES)
    real(dp) :: vely(N_NODES,N_NODES)
    real(dp) :: velz(N_NODES,N_NODES)

    real(dp) :: momx(N_NODES,N_NODES)
    real(dp) :: momy(N_NODES,N_NODES)
    real(dp) :: momz(N_NODES,N_NODES)

    real(dp) :: magx(N_NODES,N_NODES)
    real(dp) :: magy(N_NODES,N_NODES)
    real(dp) :: magz(N_NODES,N_NODES)
                     
    real(dp) :: pres(N_NODES,N_NODES)
    real(dp) :: ener(N_NODES,N_NODES)
    real(dp) :: glmp(N_NODES,N_NODES)

    real(dp) :: coords(N_NODES,N_NODES,N_DIMS)

    real(dp) :: x_shock

# if 1
    select case(side%boundary)

        case (NOR)
            !! inflow

            dens =  3.86859_dp

            velx =  11.2536_dp
            vely =  0.0_dp
            velz =  0.0_dp
            
            pres =  167.345_dp

            magx =  0.0_dp
            magy =  2.1826182_dp
            magz = -2.1826182_dp

        case (SOU)

            !! outflow

            dens = 1.0_dp

            velx = 0.0_dp
            vely = 0.0_dp
            velz = 0.0_dp

            pres = 1.0_dp

            magx = 0.0_dp
            magy = 0.56418958_dp
            magz = 0.56418958_dp

        case default

            outer = inner
            outer(:,:,GLMP_VAR) = 0.0_dp
            return
                  
    end select
# endif

# if 0
    ! x_shock = 0.05_dp + 15.176642330901245_dp * (rt%simtime + 0.5_dp*rt%timestep)
    x_shock = 0.05_dp + 15.176642330901245_dp * rt%simtime
    coords = side_get_coords(side)

    where (coords(:,:,X_DIR) < x_shock)

        dens =  3.86859_dp

        velx =  11.2536_dp
        vely =  0.0_dp
        velz =  0.0_dp
        
        pres =  167.345_dp

        magx =  0.0_dp
        magy =  2.1826182_dp
        magz = -2.1826182_dp

    else where

        dens = 1.0_dp

        velx = 0.0_dp
        vely = 0.0_dp
        velz = 0.0_dp

        pres = 1.0_dp

        magx = 0.0_dp
        magy = 0.56418958_dp
        magz = 0.56418958_dp

    end where
# endif

    glmp = 0.0_dp

    momx = dens * velx
    momy = dens * vely
    momz = dens * velz

    ener = pres/(kappa-1.0_dp) &
        + 0.5_dp*dens*(velx**2+vely**2+velz**2) &
        + 0.5_dp/mu0*(magx**2+magy**2+magz**2) &
        + 0.5_dp/mu0*glmp**2

    outer(:,:,DENS_VAR) = dens

    outer(:,:,MOMX_VAR) = momx
    outer(:,:,MOMY_VAR) = momy
    outer(:,:,MOMZ_VAR) = momz

    outer(:,:,ENER_VAR) = ener

    outer(:,:,MAGX_VAR) = magx
    outer(:,:,MAGY_VAR) = magy
    outer(:,:,MAGZ_VAR) = magz

    outer(:,:,GLMP_VAR) = glmp

end subroutine

end module
