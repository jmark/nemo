module morton_mod

use constants_mod

# if PP_N_DIMS == 2
integer(kind=4), parameter :: MORTON_MAXLEVEL = 30
integer(kind=4), parameter :: MORTON_QMAXLEVEL = 29
# else
integer(kind=4), parameter :: MORTON_MAXLEVEL = 19
integer(kind=4), parameter :: MORTON_QMAXLEVEL = 18
# endif
 
integer(kind=4), parameter :: MORTON_ROOT_LEN = SHIFTL(1,MORTON_MAXLEVEL)

integer, parameter :: MORTON_CHILDREN = 2**(N_DIMS)
contains

pure function morton_is_family(levels,mortons) result(ok)

    integer(kind=1), intent(in)     :: levels(MORTON_CHILDREN)
    integer(kind=4), intent(in)     :: mortons(N_DIMS,MORTON_CHILDREN)
    logical                         :: ok

    integer :: i

    integer(kind=4) :: inc

    ok = .false.

    if (levels(1) == 0) return

    do i = 2,MORTON_CHILDREN
        if (levels(1) /= levels(i)) then
            ok = .false.
            return
        end if
    end do

    inc = mortons(1,2) - mortons(1,1)

# if PP_N_DIMS == 2
    !! ((q0->x + inc == q1->x && q0->y       == q1->y) &&
    !!  (q0->x       == q2->x && q0->y + inc == q2->y) &&
    !!  (q1->x       == q3->x && q2->y       == q3->y))

    if (mortons(1,1) + inc == mortons(1,2) .and. mortons(2,1)       == mortons(2,2)) then
    if (mortons(1,1)       == mortons(1,3) .and. mortons(2,1) + inc == mortons(2,3)) then
    if (mortons(1,2)       == mortons(1,4) .and. mortons(2,3)       == mortons(2,4)) then
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

    if (mortons(1,1) + inc == mortons(1,2) .and. mortons(2,1)       == mortons(2,2) .and. mortons(3,1)       == mortons(3,2)) then
    if (mortons(1,1)       == mortons(1,3) .and. mortons(2,1) + inc == mortons(2,3) .and. mortons(3,1)       == mortons(3,3)) then
    if (mortons(1,2)       == mortons(1,4) .and. mortons(2,3)       == mortons(2,4) .and. mortons(3,1)       == mortons(3,4)) then
    if (mortons(1,1)       == mortons(1,5) .and. mortons(2,1)       == mortons(2,5) .and. mortons(3,1) + inc == mortons(3,5)) then
    if (mortons(1,2)       == mortons(1,6) .and. mortons(2,2)       == mortons(2,6) .and. mortons(3,5)       == mortons(3,6)) then
    if (mortons(1,3)       == mortons(1,7) .and. mortons(2,3)       == mortons(2,7) .and. mortons(3,5)       == mortons(3,7)) then
    if (mortons(1,4)       == mortons(1,8) .and. mortons(2,4)       == mortons(2,8) .and. mortons(3,5)       == mortons(3,8)) then
        ok = .true.
    end if; end if; end if; end if; end if; end if; end if
# endif

end function

subroutine morton_get_parent(clevel,cmorton,plevel,pmorton)

    !! child
    integer(kind=4), intent(in)     :: clevel
    integer(kind=4), intent(in)     :: cmorton(N_DIMS)

    !! parent
    integer(kind=4), intent(out)    :: plevel
    integer(kind=4), intent(out)    :: pmorton(N_DIMS)
  
    integer(kind=4) :: not_QUAD_LEN

    NOT_QUAD_LEN = not(SHIFTL(1,MORTON_MAXLEVEL-clevel))

    pmorton(1) = iand(cmorton(1),NOT_QUAD_LEN)
    pmorton(2) = iand(cmorton(2),NOT_QUAD_LEN)
# if PP_N_DIMS > 2
    pmorton(3) = iand(cmorton(3),NOT_QUAD_LEN)
# endif

    plevel = max(0,clevel-1)

end subroutine

pure function morton_is_ancestor(alevel,amorton,dlevel,dmorton) result(ok)

    !! ancestor
    integer(kind=4), intent(in)     :: alevel
    integer(kind=4), intent(in)     :: amorton(N_DIMS)

    !! descendant
    integer(kind=4), intent(in)    :: dlevel
    integer(kind=4), intent(in)    :: dmorton(N_DIMS)

    logical :: ok
  
    integer(kind=4) :: exclorx, exclory, exclorz

    if (alevel >= dlevel) then
        ok = .false.
        return
    end if

    exclorx = SHIFTR(IEOR(amorton(1),dmorton(1)),MORTON_MAXLEVEL - alevel)
    exclory = SHIFTR(IEOR(amorton(2),dmorton(2)),MORTON_MAXLEVEL - alevel)

# if PP_N_DIMS > 2
    exclorz = SHIFTR(IEOR(amorton(3),dmorton(3)),MORTON_MAXLEVEL - alevel)
    ok = (exclorx == 0) .and. (exclory == 0) .and. (exclorz == 0)
#else
    ok = (exclorx == 0) .and. (exclory == 0)
#endif

end function

pure function morton_get_child_id(level,morton) result(id)
    
    integer(kind=4), intent(in) :: level
    integer(kind=4), intent(in) :: morton(N_DIMS)
    
    integer :: id

    integer(kind=4) :: QUAD_LEN

    id = 0

    if (level == 0) return

    QUAD_LEN = SHIFTL(1,MORTON_MAXLEVEL-level)

    id = ior(id,merge(1,0,iand(morton(1),QUAD_LEN) /= 0))
    id = ior(id,merge(2,0,iand(morton(2),QUAD_LEN) /= 0))
# if PP_N_DIMS > 2
    id = ior(id,merge(4,0,iand(morton(3),QUAD_LEN) /= 0))
# endif

    id = id + 1

end function

pure function morton_to_vertex(morton) result(vertex)

    integer(kind=4), intent(in) :: morton(N_DIMS)
    real(dp)                    :: vertex(N_DIMS)

    integer :: i,xi,yi,zi
    real(dp) :: wx(2),wy(2),wz(2)
    real(dp) :: xfactor,yfactor

# if PP_N_DIMS == 2
    real(dp), parameter :: UNIT_VERTICES(2,4) = reshape((/&
        0.0_dp, 0.0_dp, &
        1.0_dp, 0.0_dp, &
        0.0_dp, 1.0_dp, &
        1.0_dp, 1.0_dp /), (/2,4/))

# elif PP_N_DIMS == 3
    real(dp), parameter :: UNIT_VERTICES(3,8) = reshape((/&
        0.0_dp, 0.0_dp, 0.0_dp, &
        1.0_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 1.0_dp, 0.0_dp, &
        1.0_dp, 1.0_dp, 0.0_dp, &
        0.0_dp, 0.0_dp, 1.0_dp, &
        1.0_dp, 0.0_dp, 1.0_dp, &
        0.0_dp, 1.0_dp, 1.0_dp, &
        1.0_dp, 1.0_dp, 1.0_dp /), (/3,8/))
# endif

    wx(2) = REAL(morton(1),dp) / real(MORTON_ROOT_LEN,dp)
    wx(1) = 1.0_dp - wx(2)

    wy(2) = REAL(morton(2),dp) / real(MORTON_ROOT_LEN,dp)
    wy(1) = 1.0_dp - wy(2)

# if PP_N_DIMS > 2
    wz(2) = REAL(morton(3),dp) / real(MORTON_ROOT_LEN,dp)
    wz(1) = 1.0_dp - wz(2)
# endif

    vertex = 0.0_dp
    i = 0
# if PP_N_DIMS > 2
    do zi = 1,2
# endif
        do yi = 1,2
# if PP_N_DIMS > 2
            yfactor = wz(zi) * wy(yi)
# else
            yfactor = wy(yi)
# endif
            do xi = 1,2
                i = i + 1
                xfactor = yfactor * wx(xi)
                vertex = vertex + xfactor * UNIT_VERTICES(:,i)
            end do
        end do
# if PP_N_DIMS > 2
    end do
# endif

end function

end module
