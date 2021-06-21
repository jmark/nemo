module setup_init_mod

use constants_mod

contains

subroutine fill(quad)

    use setup_config_mod
    use equations_const_mod

    use mesh_types_mod, only: quad_t
    use quad_utils_mod, only: quad_get_coords

    use math_mod, only: sigmoid

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
    real(dp) :: angles(N_NODES,N_NODES,N_NODES)

    real(dp) :: xyradi(N_NODES,N_NODES,N_NODES)

    integer :: i,j,k,l,r,m

    coords = quad_get_coords(quad)
    radius = SQRT(SUM(coords**2,dim=N_DIMS+1))
    xyradi = SQRT(coords(:,:,:,1)**2 + coords(:,:,:,2)**2)
    angles = atan2(coords(:,:,:,2),coords(:,:,:,1))

    !! original
    !! dens = merge(dens_cloud*(1+0.1*cos(2*angles)),dens0,radius < cloudradius)

    !! dens = dens0 + dens_cloud &
    !!     * sigmoid(-25*(radius/cloudradius - 1.10)) &
    !!     * sigmoid( 25*(xyradi/cloudradius - 0.40)) &
    !!     * (1+0.5*cos(2*angles))

    dens = dens0
    !! velx = velx0
    !! vely = vely0
    !! velz = velz0
    !! pres = 1.0_dp

    dens = dens + dens_cloud &
        * sigmoid(-25*(radius/cloudradius - 1.10))

    radius = SQRT((coords(:,:,:,1)-xofs)**2 + (coords(:,:,:,2))**2 + (coords(:,:,:,3))**2)
    dens = dens + dens_cloud * 1.5 * exp(-0.5*(radius/bell_sigma)**2)

    radius = SQRT((coords(:,:,:,1)+xofs)**2 + (coords(:,:,:,2))**2 + (coords(:,:,:,3))**2)
    dens = dens + dens_cloud * 1.5 * exp(-0.5*(radius/bell_sigma)**2)

    radius = SQRT(SUM(coords**2,dim=N_DIMS+1))
    velx = angl0 * (-coords(:,:,:,2))* sigmoid(-25*(radius/cloudradius - 1.2))
    vely = angl0 *   coords(:,:,:,1) * sigmoid(-25*(radius/cloudradius - 1.2))
    velz = velz0

    where (dens <= dens_critical)
        pres = K_isothermal * dens
    else where                                                                                                       
        pres = K_adiabatic * dens**kappa
    end where

    momx = dens * velx
    momy = dens * vely
    momz = dens * velz

    ener = pres/(kappa-1.0_dp) + 0.5*dens*(velx**2+vely**2+velz**2)

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
    use quad_utils_mod, only: quad_get_center
    use quad_utils_mod, only: quad_get_vertices

    type(quad_t), intent(in)        :: quad
    logical                         :: ok

    real(dp) :: verts(N_DIMS,MESH_CORNERS)
    real(dp) :: center(N_DIMS)

    verts = quad_get_vertices(quad)
    center = quad_get_center(quad)

    ok = .false.

    if (any(abs(SQRT(SUM((verts)**2,dim=1))) < 1.25*cloudradius)) then
    !if (SQRT(SUM(quad%center**2)) < 2.5*cloudradius) then
    !if (any(abs(SQRT(SUM((verts)**2,dim=1))) < 0.25 * 0.5*diameter)) then
    !if (center(1) + center(2) + center(3) < 0) then
    !if (all(center < 0)) then
        ok = .true.
    end if

end function

function setup_init_probe_refine2(quad) result(ok)

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

    if (any(abs(SQRT(SUM((verts)**2,dim=1))) < 0.02)) then
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

    use share_mod, only: rt => runtime
    use share_mod, only: mesh

    use mesh_utils_mod, only: mesh_refine
    use mesh_utils_mod, only: mesh_iterate_quads

    integer :: level

    ! do level = mesh%minlevel,mesh%maxlevel-2
    !     call mesh_refine(mesh,setup_init_probe_refine)
    ! end do
    ! call mesh_refine(mesh,setup_init_probe_refine2)

    call mesh_iterate_quads(mesh,fill)

end subroutine

end module
