module setup_init_mod

use constants_mod

contains

subroutine fill(quad)

    use setup_config_mod
    use equations_const_mod

    use mesh_types_mod, only: quad_t
    use quad_utils_mod, only: quad_get_coords
    use utils_mod, only: vmax

    type(quad_t), intent(inout) :: quad

    real(dp) :: dens(N_NODES,N_NODES,N_NODES)
    real(dp) :: velx(N_NODES,N_NODES,N_NODES)
    real(dp) :: vely(N_NODES,N_NODES,N_NODES)
    real(dp) :: velz(N_NODES,N_NODES,N_NODES)
    real(dp) :: momx(N_NODES,N_NODES,N_NODES)
    real(dp) :: momy(N_NODES,N_NODES,N_NODES)
    real(dp) :: momz(N_NODES,N_NODES,N_NODES)
    real(dp) :: pres(N_NODES,N_NODES,N_NODES)
    real(dp) :: ener(N_NODES,N_NODES,N_NODES)

    real(dp) :: coords(N_NODES,N_NODES,N_NODES,N_DIMS)
    real(dp) :: radius(N_NODES,N_NODES,N_NODES)

    coords = quad_get_coords(quad)
    !radius = SQRT(SUM(coords**2,dim=N_DIMS+1))
    radius = SQRT(coords(:,:,:,1)**2 + coords(:,:,:,2)**2 + coords(:,:,:,3)**2)

    dens = dens0
    velx = velx0
    vely = vely0
    velz = velz0
    pres = pres0

    dens = vmax(dens,merge(dc,dc*(radius/rc)**(-7),radius < rc))

    momx = dens * velx
    momy = dens * vely
    momz = dens * velz
    ener = pres/(kappa-1.0_dp) + 0.5_dp*dens*(velx**2+vely**2+velz**2)

    ener = ener + blast_ener_normalized * exp(-0.5_dp*(radius/blast_sigma)**2)

    quad%pld%hydro%state(:,:,:,DENS_VAR) = dens
    quad%pld%hydro%state(:,:,:,MOMX_VAR) = momx
    quad%pld%hydro%state(:,:,:,MOMY_VAR) = momy
    quad%pld%hydro%state(:,:,:,MOMZ_VAR) = momz
    quad%pld%hydro%state(:,:,:,ENER_VAR) = ener

end subroutine

function setup_init_probe_refine(quad) result(ok)

    use setup_config_mod
    use mesh_types_mod, only: quad_t
    use mesh_const_mod, only: MESH_CORNERS
    use quad_utils_mod, only: quad_get_vertices

    type(quad_t), intent(in)        :: quad
    logical                         :: ok

    real(dp) :: verts(N_DIMS,MESH_CORNERS)

    verts = quad_get_vertices(quad)

    if (any(abs(SQRT(SUM((verts)**2,dim=1))) < 10*blast_sigma)) then
        ok = .true.
    else
        ok = .false.
    end if

end function

subroutine hook_init_810

    use share_mod, only: mesh

    use mesh_utils_mod, only: mesh_iterate_quads
    use mesh_utils_mod, only: mesh_refine

    integer :: level
    
    do level = mesh%minlevel,mesh%maxlevel-1
        call mesh_refine(mesh,setup_init_probe_refine)
    end do
    call mesh_iterate_quads(mesh,fill)

end subroutine

end module
