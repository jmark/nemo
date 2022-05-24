module setup_init_mod

use constants_mod

contains

subroutine fill(quad)

    use setup_config_mod
    use equations_const_mod

    use mesh_types_mod, only: quad_t
    use quad_utils_mod, only: quad_get_coords
    use quad_utils_mod, only: quad_get_center

    type(quad_t), intent(inout) :: quad

    real(dp) :: dens(N_NODES,N_NODES)
    real(dp) :: velx(N_NODES,N_NODES)
    real(dp) :: vely(N_NODES,N_NODES)
    real(dp) :: momx(N_NODES,N_NODES)
    real(dp) :: momy(N_NODES,N_NODES)
    real(dp) :: pres(N_NODES,N_NODES)
    real(dp) :: ener(N_NODES,N_NODES)

    real(dp) :: center(N_DIMS)
    real(dp) :: coords(N_NODES,N_NODES,N_DIMS)
    real(dp) :: radius(N_NODES,N_NODES)

    center = quad_get_center(quad)
    coords = quad_get_coords(quad)
    radius = SQRT(SUM(coords**2,dim=N_DIMS+1))

    dens = 1.0_dp + 0.5_dp * exp(-0.5_dp*(radius/0.1_dp)**2)
    velx = 0.5_dp
    vely = 0.5_dp
    pres = 1.0_dp

    associate(x => coords(:,:,X_DIR), y => coords(:,:,Y_DIR))
    !! Primitive variables.

    dens = 1_dp*x ! + y

    ! dens = 1.0_dp

    ! dens = 10_dp*y

    ! dens = dens0 + dens1 * 0.5_dp*(1.0_dp+(tanh(slope*(y + 0.5_dp)) - (tanh(slope*(y - 0.5_dp)) + 1.0_dp)))
    ! velx = velx0 * (tanh(slope*(y + 0.5_dp)) - (tanh(slope*(y - 0.5_dp)) + 1.0_dp))
    ! vely = vely0 * sin(2.0_dp*pi*x)
    ! pres = pres0
    end associate

    !! Conservative variables.
    momx = dens * velx
    momy = dens * vely
    ener = pres/(kappa-1.0_dp) + 0.5_dp*dens*(velx**2+vely**2)

    ! quad%pld%hydro%state = 42.0_dp

    quad%pld%hydro%state(:,:,DENS_VAR) = dens
    quad%pld%hydro%state(:,:,MOMX_VAR) = momx
    quad%pld%hydro%state(:,:,MOMY_VAR) = momy
    quad%pld%hydro%state(:,:,MOMZ_VAR) = 0.0_dp
    quad%pld%hydro%state(:,:,ENER_VAR) = ener

    ! quad%pld%hydro%state(:,:,DENS_VAR) = 5 - 10*(coords(:,:,1) + coords(:,:,2))**2
    ! quad%pld%hydro%state(:,:,DENS_VAR) = (coords(:,:,1) + coords(:,:,2))**5

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

    center = quad_get_center(quad)

    ok = .false.

    ! return
    !if (rand() < 0.2) then
    if (any(abs(SQRT(SUM((verts)**2,dim=1))) <= refine_radius)) then
    !if (center(1) < 0) then
    !if (center(2) < 0) then
    !if (all(center < 0)) then
        ok = .true.
    end if

end function

function setup_init_probe_coarsen(quads) result(ok)

    use mesh_types_mod, only: quad_t
    use mesh_const_mod, only: MESH_CHILDREN

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
    
    ! do level = mesh%minlevel,mesh%maxlevel-1
    !     call mesh_refine(mesh,setup_init_probe_refine)
    ! end do
    call mesh_iterate_quads(mesh,fill)

end subroutine

end module
