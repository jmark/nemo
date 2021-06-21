module setup_init_mod

use constants_mod

contains

subroutine fill(quad)

    use setup_config_mod
    use equations_const_mod

    use mesh_types_mod, only: quad_t
    use quad_utils_mod, only: quad_get_coords

    type(quad_t), intent(inout) :: quad

    real(dp) :: dens(N_NODES,N_NODES)
    real(dp) :: velx(N_NODES,N_NODES)
    real(dp) :: vely(N_NODES,N_NODES)
    real(dp) :: momx(N_NODES,N_NODES)
    real(dp) :: momy(N_NODES,N_NODES)
    real(dp) :: pres(N_NODES,N_NODES)
    real(dp) :: ener(N_NODES,N_NODES)

    real(dp) :: coords(N_NODES,N_NODES,N_DIMS)
    real(dp) :: radius(N_NODES,N_NODES)

    coords = quad_get_coords(quad)
    radius = SQRT(SUM((coords)**2,dim=N_DIMS+1))

    dens = dens0

    !dens = dens + blast_ener_normalized * exp(-0.5_dp*(radius/blast_sigma)**2)
    !dens = dens + exp(-0.5_dp*(radius/0.01_dp)**2)

    velx = velx0
    vely = vely0
    pres = pres0

    momx = dens * velx
    momy = dens * vely
    ener = pres/(kappa-1.0_dp) + 0.5_dp*dens*(velx**2+vely**2)

    ener = ener + blast_ener_normalized * exp(-0.5_dp*(radius/blast_sigma)**2)

    quad%pld%hydro%state(:,:,DENS_VAR) = dens
    quad%pld%hydro%state(:,:,MOMX_VAR) = momx
    quad%pld%hydro%state(:,:,MOMY_VAR) = momy
    quad%pld%hydro%state(:,:,ENER_VAR) = ener

end subroutine

function setup_init_probe_refine(quad) result(ok)

    use setup_config_mod
    use mesh_types_mod, only: quad_t
    use mesh_const_mod, only: MESH_CORNERS
    use quad_utils_mod, only: quad_get_center
    use quad_utils_mod, only: quad_get_vertices

    type(quad_t), intent(in)        :: quad
    logical                         :: ok

    real(dp) :: verts(N_DIMS,MESH_CORNERS)
    real(dp) :: center(N_DIMS)

    verts = quad_get_vertices(quad)
    ! center = quad_get_center(quad)

    ok = .false.

    if (any(abs(SQRT(SUM((verts)**2,dim=1))) < init_amr_min_diameter)) then
        ok = .true.
    end if

end function

function setup_init_probe_coarsen(quads) result(ok)

    use mesh_types_mod, only: quad_t
    use mesh_const_mod, only: MESH_CHILDREN
    use setup_config_mod, only: setup_amr_threshold_coarsen

    type(quad_t), intent(in)    :: quads(MESH_CHILDREN)
    logical                     :: ok

    real(dp)                    :: sensor(MESH_CHILDREN)

    integer :: i

    ok = .false.

end function

subroutine hook_init_810

    use share_mod, only: mesh

    use mesh_utils_mod, only: mesh_iterate_quads
    use mesh_utils_mod, only: mesh_refine

    integer :: level
    
    do level = mesh%minlevel,mesh%maxlevel-1
        call mesh_refine(mesh,setup_init_probe_refine)
        call mesh_iterate_quads(mesh,fill)
    end do
    call mesh_iterate_quads(mesh,fill)

end subroutine

end module
