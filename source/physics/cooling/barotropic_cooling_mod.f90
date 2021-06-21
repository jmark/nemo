module barotropic_cooling_mod

use constants_mod

# if DO_PROFILING
    integer :: profid
# endif

contains

subroutine hook_init_914

# if DO_PROFILING
    use profiler_mod, only: profiler_init
    profid = profiler_init('barotropic_cooling')
# endif

end subroutine

subroutine hook_post_rk_timestep_200

    use share_mod, only: mesh
    use mesh_utils_mod, only: mesh_iterate_quads

# if DO_PROFILING
    use profiler_mod, only: profiler_start,profiler_stop

    call profiler_start(profid)
# endif

    call mesh_iterate_quads(mesh,barotropic_cooling_cb)

# if DO_PROFILING
    call profiler_stop(profid)
# endif

end subroutine

subroutine barotropic_cooling_cb(quad)

    use setup_config_mod
    use equations_const_mod
    use mesh_types_mod, only: quad_t

    type(quad_t), intent(inout) :: quad

    where (quad%pld%hydro%state(:,:,:,DENS_VAR) <= dens_critical)
        quad%pld%hydro%state(:,:,:,ENER_VAR) = 0.5/quad%pld%hydro%state(:,:,:,DENS_VAR)*(quad%pld%hydro%state(:,:,:,MOMX_VAR)**2 + quad%pld%hydro%state(:,:,:,MOMY_VAR)**2 + quad%pld%hydro%state(:,:,:,MOMZ_VAR)**2) + (K_isothermal/(kappa-1.0_dp))*quad%pld%hydro%state(:,:,:,DENS_VAR)
    else where                                                                                                       
        quad%pld%hydro%state(:,:,:,ENER_VAR) = 0.5/quad%pld%hydro%state(:,:,:,DENS_VAR)*(quad%pld%hydro%state(:,:,:,MOMX_VAR)**2 + quad%pld%hydro%state(:,:,:,MOMY_VAR)**2 + quad%pld%hydro%state(:,:,:,MOMZ_VAR)**2) + (K_adiabatic /(kappa-1.0_dp))*quad%pld%hydro%state(:,:,:,DENS_VAR)**kappa
    end where

end subroutine

end module
