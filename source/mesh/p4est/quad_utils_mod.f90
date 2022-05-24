# if PP_N_DIMS == 3
# define p4est_constants_mod p8est_constants_mod
# define P4EST_MAXLEVEL P8EST_MAXLEVEL
# endif

module quad_utils_mod

use constants_mod
use mesh_types_mod
use mesh_const_mod
use kernel_share_mod, only: nodes, weights

public

! real(dp) :: nodes(N_NODES), weights(N_NODES)

# if PP_N_DIMS == 2
real(dp) :: weightsD(N_NODES,N_NODES)
# elif PP_N_DIMS == 3
real(dp) :: weightsD(N_NODES,N_NODES,N_NODES)
# endif

contains

subroutine hook_init_510

    ! use interpol_mod, only: EquidistantNodesAndWeights

    integer :: n1,n2,n3

    !! The prefactor is necessary since the weights
    !! are defined for the intervall: [-1,1]
    ! real(dp), parameter :: prefactor = 0.5**N_DIMS

    ! call EquidistantNodesAndWeights(N_NODES-1, nodes, weights)

# if PP_N_DIMS == 2
    do n2 = 1,N_NODES; do n1 = 1,N_NODES;
        weightsD(n1,n2) = weights(n1)*weights(n2)
    end do; end do;
# elif PP_N_DIMS == 3
    do n3 = 1,N_NODES; do n2 = 1,N_NODES; do n1 = 1,N_NODES;
        weightsD(n1,n2,n3) = weights(n1)*weights(n2)*weights(n3)
    end do; end do; end do;
# endif

end subroutine

pure function quad_get_center(quad) result(center)

    type(quad_t), intent(in)    :: quad
    real(dp)                    :: verts(N_DIMS,MESH_CORNERS)
    real(dp)                    :: center(N_DIMS)

    verts = quad_get_vertices(quad)

# if PP_N_DIMS == 2
    center = 0.5*(verts(:,NWF) + verts(:,SEF))
# elif PP_N_DIMS == 3
    center = 0.5*(verts(:,NWF) + verts(:,SEB))
# endif

end function

# if PP_N_DIMS == 2
pure function quad_get_vertices(quad) result(verts)

    use setup_config_mod
    use p4est_interfaces_mod, only: p4est_qcoord_to_vertex

    type(quad_t), intent(in)    :: quad
    real(dp)                    :: verts(N_DIMS,MESH_CORNERS)

    real(dp), parameter :: upleft(N_DIMS) = (/xlim(1),ylim(1)/)
    real(dp), parameter :: cublen(N_DIMS) = (/xlim(2)-xlim(1),ylim(2)-ylim(1)/)

    real(dp) :: difvec(N_DIMS)
    real(dp) :: coords(3)

    difvec = 1.0_dp / 2**quad%level

    call p4est_qcoord_to_vertex(quad%mesh%conn,0,quad%morton(1),quad%morton(2),coords)

    verts(:,NWF) = upleft + cublen * (coords(1:2) + (/0.0_dp,0.0_dp/) * difvec)
    verts(:,NEF) = upleft + cublen * (coords(1:2) + (/0.0_dp,1.0_dp/) * difvec)
    verts(:,SWF) = upleft + cublen * (coords(1:2) + (/1.0_dp,0.0_dp/) * difvec)
    verts(:,SEF) = upleft + cublen * (coords(1:2) + (/1.0_dp,1.0_dp/) * difvec)

end function

# elif PP_N_DIMS == 3

pure function quad_get_vertices(quad) result(verts)

    use setup_config_mod, only: xlim,ylim,zlim
    use p8est_interfaces_mod, only: p8est_qcoord_to_vertex

    type(quad_t), intent(in)    :: quad
    real(dp)                    :: verts(N_DIMS,MESH_CORNERS)

    real(dp), parameter :: upleft(N_DIMS) = (/xlim(1),ylim(1),zlim(1)/)
    real(dp), parameter :: cublen(N_DIMS) = (/xlim(2)-xlim(1),ylim(2)-ylim(1),zlim(2)-zlim(1)/)

    real(dp) :: difvec(N_DIMS)
    real(dp) :: coords(3)

    difvec = 1.0_dp / REAL(2**quad%level,dp)

    call p8est_qcoord_to_vertex(quad%mesh%conn,0,quad%morton(1),quad%morton(2),quad%morton(3),coords)

    verts(:,NWF) = upleft + cublen * (coords + (/0.0_dp,0.0_dp,0.0_dp/) * difvec)
    verts(:,NEF) = upleft + cublen * (coords + (/0.0_dp,1.0_dp,0.0_dp/) * difvec)
    verts(:,SWF) = upleft + cublen * (coords + (/1.0_dp,0.0_dp,0.0_dp/) * difvec)
    verts(:,SEF) = upleft + cublen * (coords + (/1.0_dp,1.0_dp,0.0_dp/) * difvec)

    verts(:,NWB) = upleft + cublen * (coords + (/0.0_dp,0.0_dp,1.0_dp/) * difvec)
    verts(:,NEB) = upleft + cublen * (coords + (/0.0_dp,1.0_dp,1.0_dp/) * difvec)
    verts(:,SWB) = upleft + cublen * (coords + (/1.0_dp,0.0_dp,1.0_dp/) * difvec)
    verts(:,SEB) = upleft + cublen * (coords + (/1.0_dp,1.0_dp,1.0_dp/) * difvec)

end function
# endif

pure function quad_get_delta(quad) result(delta)

    type(quad_t), intent(in)    :: quad
    real(dp)                    :: delta(N_DIMS)
    real(dp)                    :: verts(N_DIMS,MESH_CORNERS)

    verts = quad_get_vertices(quad)

    delta(1) = abs(verts(1,SWF)-verts(1,NWF))
    delta(2) = abs(verts(2,NEF)-verts(2,NWF))

# if PP_N_DIMS == 3
    delta(3) = abs(verts(3,NEB)-verts(3,NWF))
# endif

end function

pure function quad_get_nodes(quad) result(ns)

    type(quad_t), intent(in)    :: quad
    real(dp)                    :: ns(N_NODES)

    ns = nodes

end function

pure function quad_get_weights(quad) result(ws)

    type(quad_t), intent(in)    :: quad
# if PP_N_DIMS == 2
    real(dp)                    :: ws(N_NODES,N_NODES)
# elif PP_N_DIMS == 3
    real(dp)                    :: ws(N_NODES,N_NODES,N_NODES)
# endif

    ws = weightsD

end function

pure function quad_get_volumes(quad) result(volumes)

    type(quad_t), intent(in)    :: quad
# if PP_N_DIMS == 2
    real(dp)                    :: volumes(N_NODES,N_NODES)
# elif PP_N_DIMS == 3
    real(dp)                    :: volumes(N_NODES,N_NODES,N_NODES)
# endif

    volumes = product(quad%delta)*weightsD

end function

pure function quad_get_coords(quad) result(coords)

    type(quad_t), intent(in)    :: quad
# if PP_N_DIMS == 2
    real(dp)                    :: coords(N_NODES,N_NODES,N_DIMS)
# elif PP_N_DIMS == 3
    real(dp)                    :: coords(N_NODES,N_NODES,N_NODES,N_DIMS)
# endif

    real(dp) :: verts(N_DIMS,MESH_CORNERS),delta(N_DIMS)

    integer :: n1,n2,n3

    verts = quad_get_vertices(quad)

# if PP_N_DIMS == 2
    delta = abs(verts(:,SEF)-verts(:,NWF))
    do n2 = 1,N_NODES; do n1 = 1,N_NODES;
        coords(n1,n2,1) = verts(1,NWF) + (nodes(n1)+0.5_dp)*delta(1)
        coords(n1,n2,2) = verts(2,NWF) + (nodes(n2)+0.5_dp)*delta(2)
    end do; end do;
# elif PP_N_DIMS == 3
    delta = abs(verts(:,SEB)-verts(:,NWF))
    do n3 = 1,N_NODES; do n2 = 1,N_NODES; do n1 = 1,N_NODES;
        coords(n1,n2,n3,1) = verts(1,NWF) + (nodes(n1)+0.5_dp)*delta(1)
        coords(n1,n2,n3,2) = verts(2,NWF) + (nodes(n2)+0.5_dp)*delta(2)
        coords(n1,n2,n3,3) = verts(3,NWF) + (nodes(n3)+0.5_dp)*delta(3)
    end do; end do; end do;
# endif

end function

pure function quad_get_face_coords(quad,faceid) result(coords)

    type(quad_t), intent(in)    :: quad
    integer, intent(in)         :: faceid
# if PP_N_DIMS == 2
    real(dp)                    :: coords(N_NODES,N_DIMS)
# elif PP_N_DIMS == 3
    real(dp)                    :: coords(N_NODES,N_NODES,N_DIMS)
# endif

    real(dp)    :: verts(N_DIMS,MESH_CORNERS),delta(N_DIMS)
    integer     :: n1,n2

    verts = quad_get_vertices(quad)

# if PP_N_DIMS == 2
    delta = abs(verts(:,SEF)-verts(:,NWF))

    select case(faceid)
        case(NOR)
            do n1 = 1,N_NODES
                coords(n1,1) = verts(1,NWF)
                coords(n1,2) = verts(2,NWF) + (nodes(n1)+0.5_dp)*delta(2)
            end do;

        case(SOU)
            do n1 = 1,N_NODES
                coords(n1,1) = verts(1,SWF)
                coords(n1,2) = verts(2,SWF) + (nodes(n1)+0.5_dp)*delta(2)
            end do;

        case(WES)
            do n1 = 1,N_NODES
                coords(n1,1) = verts(1,NWF) + (nodes(n1)+0.5_dp)*delta(1)
                coords(n1,2) = verts(2,NWF)
            end do;

        case(EAS)
            do n1 = 1,N_NODES
                coords(n1,1) = verts(1,NWF) + (nodes(n1)+0.5_dp)*delta(1)
                coords(n1,2) = verts(2,NEF)
            end do;
    end select

# elif PP_N_DIMS == 3
    delta = abs(verts(:,SEB)-verts(:,NWF))

    select case(faceid)
        case(NOR)
            do n2 = 1,N_NODES; do n1 = 1,N_NODES
                coords(n1,n2,1) = verts(1,NWF)
                coords(n1,n2,2) = verts(2,NWF) + (nodes(n1)+0.5_dp)*delta(2)
                coords(n1,n2,3) = verts(3,NWF) + (nodes(n2)+0.5_dp)*delta(3)
            end do; end do;

        case(SOU)
            do n2 = 1,N_NODES; do n1 = 1,N_NODES
                coords(n1,n2,1) = verts(1,SWF)
                coords(n1,n2,2) = verts(2,SWF) + (nodes(n1)+0.5_dp)*delta(2)
                coords(n1,n2,3) = verts(3,SWF) + (nodes(n2)+0.5_dp)*delta(3)
            end do; end do;

        case(WES)
            do n2 = 1,N_NODES; do n1 = 1,N_NODES
                coords(n1,n2,1) = verts(1,NWF) + (nodes(n1)+0.5_dp)*delta(1)
                coords(n1,n2,2) = verts(2,NWF)
                coords(n1,n2,3) = verts(3,NWF) + (nodes(n2)+0.5_dp)*delta(3)
            end do; end do;

        case(EAS)
            do n2 = 1,N_NODES; do n1 = 1,N_NODES
                coords(n1,n2,1) = verts(1,NEF) + (nodes(n1)+0.5_dp)*delta(1)
                coords(n1,n2,2) = verts(2,NEF)
                coords(n1,n2,3) = verts(3,NEF) + (nodes(n2)+0.5_dp)*delta(3)
            end do; end do;

        case(FRO)
            do n2 = 1,N_NODES; do n1 = 1,N_NODES
                coords(n1,n2,1) = verts(1,NWF) + (nodes(n1)+0.5_dp)*delta(1)
                coords(n1,n2,2) = verts(2,NWF) + (nodes(n2)+0.5_dp)*delta(2)
                coords(n1,n2,3) = verts(3,NWF)
            end do; end do;

        case(BAC)
            do n2 = 1,N_NODES; do n1 = 1,N_NODES
                coords(n1,n2,1) = verts(1,NWB) + (nodes(n1)+0.5_dp)*delta(1)
                coords(n1,n2,2) = verts(2,NWB) + (nodes(n2)+0.5_dp)*delta(2)
                coords(n1,n2,3) = verts(3,NWB)
            end do; end do;
    end select
# endif

end function

!! ------------------------------------------------------------------------- !!

pure function quad_is_equal(a,b) result(yes)

    type(quad_t), intent(in) :: a,b
    logical                  :: yes

    yes = (a%level == b%level) .and. all(a%morton == b%morton)

end function

!! ------------------------------------------------------------------------- !!
!! ------------------------------------------------------------------------- !!

pure function vert_get_dist(a,b) result(dist)

    type(vert_t), intent(in)    :: a,b
    real(dp)                    :: dist
    real(dp)                    :: v(N_DIMS)

    v = vert_get_coords(b) - vert_get_coords(a)
    dist = SQRT(SUM(v*v))

end function

pure function vert_get_coords(vert) result(coords)

    use mesh_const_mod, only: map_to_opposite_vertid

    type(vert_t), intent(in)    :: vert
    real(dp)                    :: coords(N_DIMS)
    real(dp)                    :: verts(N_DIMS,MESH_CORNERS)

    verts = quad_get_vertices(vert%mesh%quads(vert%iquads(vert%refquad)))
    coords = verts(:,map_to_opposite_vertid(vert%refquad))

end function

pure function side_get_coords(side) result(coords)

    type(side_t), intent(in)    :: side

# if PP_N_DIMS == 2
    real(dp)                    :: coords(N_NODES,N_DIMS)
# elif PP_N_DIMS == 3
    real(dp)                    :: coords(N_NODES,N_NODES,N_DIMS)
# endif

    coords = quad_get_face_coords(side%quads(1,side%inner)%p,side%faceid(side%inner))

end function

function quad_is_family(family) result(ok)

    use mesh_const_mod, only: MESH_CHILDREN

    use p4est_constants_mod, only: P4EST_MAXLEVEL
    ! use p4est_interfaces_mod, only: p4est_quadrant_is_family

    type(quad_t), intent(in) :: family(MESH_CHILDREN)
    logical                  :: ok

    integer :: i

    integer(kind=4) :: inc

    ok = .false.

    if (family(1)%level == 0) return

    do i = 2,MESH_CHILDREN
        if (family(1)%level /= family(i)%level) then
            ok = .false.
            return
        end if
    end do

    !! #define P4EST_QUADRANT_LEN(l) ((p4est_qcoord_t) 1 << (P4EST_MAXLEVEL - (l)))
    !inc = LSHIFT(1,P4EST_MAXLEVEL-family(1)%level)
    inc = family(2)%morton(1) - family(1)%morton(1)

    ! write (*,*) 'first level!', inc, family(1)%morton(1), family(2)%morton(1), family(2)%morton(1) - family(1)%morton(1)

# if PP_N_DIMS == 2
    !! ((q0->x + inc == q1->x && q0->y       == q1->y) &&
    !!  (q0->x       == q2->x && q0->y + inc == q2->y) &&
    !!  (q1->x       == q3->x && q2->y       == q3->y))

    if (family(1)%morton(1) + inc == family(2)%morton(1) .and. family(1)%morton(2)       == family(2)%morton(2)) then
    if (family(1)%morton(1)       == family(3)%morton(1) .and. family(1)%morton(2) + inc == family(3)%morton(2)) then
    if (family(2)%morton(1)       == family(4)%morton(1) .and. family(3)%morton(2)       == family(4)%morton(2)) then
        ok = .true.
    end if; end if; end if

# elif PP_N_DIMS == 3
    !! ((q0->x + inc == q1->x && q0->y       == q1->y && q0->z       == q1->z) &&
    !!  (q0->x       == q2->x && q0->y + inc == q2->y && q0->z       == q2->z) &&
    !!  (q1->x       == q3->x && q2->y       == q3->y && q0->z       == q3->z) &&
    !!  (q0->x       == q4->x && q0->y       == q4->y && q0->z + inc == q4->z) &&
    !!  (q1->x       == q5->x && q1->y       == q5->y && q4->z       == q5->z) &&
    !!  (q2->x       == q6->x && q2->y       == q6->y && q4->z       == q6->z) &&
    !!  (q3->x       == q7->x && q3->y       == q7->y && q4->z       == q7->z))

    if (family(1)%morton(1) + inc == family(2)%morton(1) .and. family(1)%morton(2)       == family(2)%morton(2) .and. family(1)%morton(3)       == family(2)%morton(3)) then
        !write (*,*) 'check 1'
    if (family(1)%morton(1)       == family(3)%morton(1) .and. family(1)%morton(2) + inc == family(3)%morton(2) .and. family(1)%morton(3)       == family(3)%morton(3)) then
        !write (*,*) 'check 2'
    if (family(2)%morton(1)       == family(4)%morton(1) .and. family(3)%morton(2)       == family(4)%morton(2) .and. family(1)%morton(3)       == family(4)%morton(3)) then
        !write (*,*) 'check 3'
    if (family(1)%morton(1)       == family(5)%morton(1) .and. family(1)%morton(2)       == family(5)%morton(2) .and. family(1)%morton(3) + inc == family(5)%morton(3)) then
        !write (*,*) 'check 4'
    if (family(2)%morton(1)       == family(6)%morton(1) .and. family(2)%morton(2)       == family(6)%morton(2) .and. family(5)%morton(3)       == family(6)%morton(3)) then
        !write (*,*) 'check 5'
    if (family(3)%morton(1)       == family(7)%morton(1) .and. family(3)%morton(2)       == family(7)%morton(2) .and. family(5)%morton(3)       == family(7)%morton(3)) then
        !write (*,*) 'check 6'
    if (family(4)%morton(1)       == family(8)%morton(1) .and. family(4)%morton(2)       == family(8)%morton(2) .and. family(5)%morton(3)       == family(8)%morton(3)) then
        !write (*,*) 'check 7'
        ok = .true.
    end if; end if; end if; end if; end if; end if; end if

# endif

end function

end module
