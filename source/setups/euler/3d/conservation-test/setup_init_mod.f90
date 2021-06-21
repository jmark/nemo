module setup_init_mod

use constants_mod

real(dp), save :: nodes(N_NODES),ws(N_NODES)

contains

subroutine hook_init_809

    use share_mod, only: mesh

    use interpol_mod, only: LegendreGaussNodesAndWeights

    integer :: level

    call LegendreGaussNodesAndWeights(N_NODES-1,nodes,ws)
    nodes = 0.5_dp*nodes
    ws = 0.5_dp*ws

end subroutine

subroutine fill(quad)

    use setup_config_mod
    use equations_const_mod
    use setup_manufacsol_mod

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
    real(dp) :: pres(N_NODES,N_NODES,N_NODES)
    real(dp) :: ener(N_NODES,N_NODES,N_NODES)

    real(dp) :: coords(N_NODES,N_NODES,N_NODES,N_DIMS)

    real(dp) :: u(N_NODES,N_NODES,N_NODES,N_VARS)
    real(dp) :: r(N_NODES,N_NODES,N_NODES,N_VARS)

    real(dp) :: xs(N_NODES)
    real(dp) :: ys(N_NODES)
    real(dp) :: zs(N_NODES)

    integer :: i,j,k, ii,jj,kk

    coords = quad_get_coords(quad)

    u = 0.0_dp
    r = 0.0_dp

    do k = 1,N_NODES; do j = 1,N_NODES; do i = 1,N_NODES;
        xs = coords(i,j,k,1) + quad%delta(1)/REAL(N_NODES,dp)*nodes
        ys = coords(i,j,k,2) + quad%delta(2)/REAL(N_NODES,dp)*nodes
        zs = coords(i,j,k,3) + quad%delta(3)/REAL(N_NODES,dp)*nodes

        do kk = 1,N_NODES; do jj = 1,N_NODES; do ii = 1,N_NODES;
            u(i,j,k,:) = u(i,j,k,:) + ws(ii)*ws(jj)*ws(kk)*init_manufacsol(xs(ii),ys(jj),zs(kk))
            r(i,j,k,:) = r(i,j,k,:) + ws(ii)*ws(jj)*ws(kk)*calc_manufacsol(xs(ii),ys(jj),zs(kk))
        end do; end do; end do
    end do; end do; end do

    dens = u(:,:,:,DENS_VAR)
    velx = u(:,:,:,VELX_VAR)
    vely = u(:,:,:,VELY_VAR)
    velz = u(:,:,:,VELZ_VAR)
    pres = u(:,:,:,PRES_VAR)

    momx = dens * velx
    momy = dens * vely
    momz = dens * velz
    ener = pres/(kappa-1) + 0.5*dens*(velx**2+vely**2+velz**2)

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
    ! center = quad_get_center(quad)

    ok = .false.
    if (any(abs(SQRT(SUM((verts)**2,dim=1))) < 0.125*diameter)) then
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
    
    call srand(86456)

    do level = mesh%minlevel,mesh%maxlevel-1
        call mesh_refine(mesh,setup_init_probe_refine)
    end do
    call mesh_iterate_quads(mesh,fill)

end subroutine

end module
