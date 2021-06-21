module setup_init_mod

use constants_mod

contains

subroutine fill(quad)

    use setup_config_mod
    use equations_const_mod

    use mesh_types_mod, only: quad_t
    use quad_utils_mod, only: quad_get_coords

    type(quad_t), intent(inout) :: quad

    real(dp) :: dens(N_NODES,N_NODES,N_NODES)

    real(dp) :: velx(N_NODES,N_NODES,N_NODES)
    real(dp) :: vely(N_NODES,N_NODES,N_NODES)
    real(dp) :: velz(N_NODES,N_NODES,N_NODES)

    real(dp) :: momx(N_NODES,N_NODES,N_NODES)
    real(dp) :: momy(N_NODES,N_NODES,N_NODES)
    real(dp) :: momz(N_NODES,N_NODES,N_NODES)

    real(dp) :: magx(N_NODES,N_NODES,N_NODES)
    real(dp) :: magy(N_NODES,N_NODES,N_NODES)
    real(dp) :: magz(N_NODES,N_NODES,N_NODES)
                             
    real(dp) :: pres(N_NODES,N_NODES,N_NODES)
    real(dp) :: ener(N_NODES,N_NODES,N_NODES)
    real(dp) :: glmp(N_NODES,N_NODES,N_NODES)

    real(dp) :: coords(N_NODES,N_NODES,N_NODES,N_DIMS)
    real(dp) :: radius(N_NODES,N_NODES,N_NODES)

    coords = quad_get_coords(quad)

    where (coords(:,:,:,X_DIR) < 0.05_dp)

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

    glmp = 0.0_dp

    radius = SQRT(&
        (coords(:,:,:,X_DIR) - 0.25_dp)**2 + &
        (coords(:,:,:,Y_DIR) - 0.5_dp)**2 + &
        (coords(:,:,:,Z_DIR) - 0.5_dp)**2)

    where (radius <= 0.15_dp)
        dens = 10.0_dp
    end where

    momx = dens * velx
    momy = dens * vely
    momz = dens * velz

    ener = pres/(kappa-1.0_dp) &
        + 0.5_dp*dens*(velx**2+vely**2+velz**2) &
        + 0.5_dp/mu0*(magx**2+magy**2+magz**2) &
        + 0.5_dp/mu0*glmp**2

    quad%pld%hydro%state(:,:,:,DENS_VAR) = dens

    quad%pld%hydro%state(:,:,:,MOMX_VAR) = momx
    quad%pld%hydro%state(:,:,:,MOMY_VAR) = momy
    quad%pld%hydro%state(:,:,:,MOMZ_VAR) = momz

    quad%pld%hydro%state(:,:,:,ENER_VAR) = ener

    quad%pld%hydro%state(:,:,:,MAGX_VAR) = magx
    quad%pld%hydro%state(:,:,:,MAGY_VAR) = magy
    quad%pld%hydro%state(:,:,:,MAGZ_VAR) = magz
                           
    quad%pld%hydro%state(:,:,:,GLMP_VAR) = glmp

end subroutine

subroutine hook_init_810

    use share_mod, only: mesh

    use mesh_utils_mod, only: mesh_iterate_quads
    use setup_amr_mod, only: setup_amr_adapt

    integer :: level
    
    do level = mesh%minlevel,mesh%maxlevel-1
        call mesh_iterate_quads(mesh,fill)
        call setup_amr_adapt()
    end do
    call mesh_iterate_quads(mesh,fill)

end subroutine

end module
