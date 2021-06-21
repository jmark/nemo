!# define MODULE_SOURCETERM 1

module setup_sourceterm_mod

use constants_mod

real(dp), save :: nodes(N_NODES),ws(N_NODES)

contains

subroutine hook_init_820

    use share_mod, only: mesh
    use interpol_mod, only: LegendreGaussNodesAndWeights

    integer :: level

    call LegendreGaussNodesAndWeights(N_NODES-1,nodes,ws)
    nodes = 0.5_dp*nodes
    ws = 0.5_dp*ws

end subroutine

pure subroutine calc_sourceterm(quad,rhs)

    use equations_const_mod
    use mesh_types_mod, only: quad_t
    use quad_utils_mod, only: quad_get_coords
    use setup_manufacsol_mod

    use share_mod, only: rt => runtime

    type(quad_t), intent(inout) :: quad
    real(dp), intent(out)       :: rhs(N_NODES,N_NODES,N_NODES,N_VARS)

    real(dp) :: coords(N_NODES,N_NODES,N_NODES,N_DIMS)

    real(dp) :: xs(N_NODES)
    real(dp) :: ys(N_NODES)
    real(dp) :: zs(N_NODES)

    real(dp) :: time

    integer :: i,j,k, ii,jj,kk

    ! rhs = quad%aux%manufacsol%st
    ! return

    coords = quad_get_coords(quad)

    time = rt%rksimtime
    ! time = 0.0_dp

    rhs = 0.0_dp
    do k = 1,N_NODES; do j = 1,N_NODES; do i = 1,N_NODES;
        xs = coords(i,j,k,1) + quad%delta(1)/REAL(N_NODES,dp)*nodes
        ys = coords(i,j,k,2) + quad%delta(2)/REAL(N_NODES,dp)*nodes
        zs = coords(i,j,k,3) + quad%delta(3)/REAL(N_NODES,dp)*nodes

        do kk = 1,N_NODES; do jj = 1,N_NODES; do ii = 1,N_NODES;
            rhs(i,j,k,:) = rhs(i,j,k,:) + ws(ii)*ws(jj)*ws(kk)*calc_manufacsol(time,xs(ii),ys(jj),zs(kk))
        end do; end do; end do
    end do; end do; end do

end subroutine

function setup_sourceterm_dt(quad) result(dt)

    use mesh_types_mod, only: quad_t

    type(quad_t), intent(in)    :: quad
    real(dp)                    :: dt

    dt = HUGE(1.0_dp)

end function

end module
