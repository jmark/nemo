module kernel_flux_mod

use constants_mod

use mesh_types_mod, only: quad_t

!! ------------------------------------------------------------------------- !!

contains

subroutine kernel_fill_halo(quad)

    use mesh_types_mod, only: quad_t
    use mesh_types_mod, only: side_t
    use mesh_types_mod, only: vert_t

    use share_mod, only: mesh

    type(quad_t), intent(inout) :: quad

    type(side_t), pointer :: side
    type(vert_t), pointer :: vert
    type(quad_t), pointer :: neig

    integer :: i,j,h

    integer, parameter :: N = N_NODES

    associate(states => quad%pld%hydro%state, halo => quad%aux%hydro%halo)
    do h = 1,N_VARS; do i = 1,N_NODES; do j = 1,N_NODES;
        halo(h,i,j) = states(i,j,h)
    end do; end do; end do;

    !! north
    side => mesh%sides(quad%isides(NOR))
    neig => mesh%quads(side%iquads(1,side%plus))
    do h = 1,N_VARS; do i = 1,N; do j = 1,N;
        halo(h,i-N,j) = neig%pld%hydro%state(i,j,h)
    end do; end do; end do;

    !! south
    side => mesh%sides(quad%isides(SOU))
    neig => mesh%quads(side%iquads(1,side%minus))
    do h = 1,N_VARS; do i = 1,N; do j = 1,N;
        halo(h,i+N,j) = neig%pld%hydro%state(i,j,h)
    end do; end do; end do;

    !! west
    side => mesh%sides(quad%isides(WES))
    neig => mesh%quads(side%iquads(1,side%plus))
    do h = 1,N_VARS; do i = 1,N; do j = 1,N;
        halo(h,i,j-N) = neig%pld%hydro%state(i,j,h)
    end do; end do; end do;

    !! east
    side => mesh%sides(quad%isides(EAS))
    neig => mesh%quads(side%iquads(1,side%minus))
    do h = 1,N_VARS; do i = 1,N; do j = 1,N;
        halo(h,i,j+N) = neig%pld%hydro%state(i,j,h)
    end do; end do; end do;

    !! north west
    vert => mesh%verts(quad%iverts(NWF))
    neig => mesh%quads(vert%iquads(NWF))
    do h = 1,N_VARS; do i = 1,N; do j = 1,N;
        halo(h,i-N,j-N) = neig%pld%hydro%state(i,j,h)
    end do; end do; end do;

    !! north east
    vert => mesh%verts(quad%iverts(NEF))
    neig => mesh%quads(vert%iquads(NEF))
    do h = 1,N_VARS; do i = 1,N; do j = 1,N;
        halo(h,i-N,j+N) = neig%pld%hydro%state(i,j,h)
    end do; end do; end do;

    !! south west
    vert => mesh%verts(quad%iverts(SWF))
    neig => mesh%quads(vert%iquads(SWF))
    do h = 1,N_VARS; do i = 1,N; do j = 1,N;
        halo(h,i+N,j-N) = neig%pld%hydro%state(i,j,h)
    end do; end do; end do;

    !! south east
    vert => mesh%verts(quad%iverts(SEF))
    neig => mesh%quads(vert%iquads(SEF))
    do h = 1,N_VARS; do i = 1,N; do j = 1,N;
        halo(h,i+N,j+N) = neig%pld%hydro%state(i,j,h)
    end do; end do; end do;
    end associate

end subroutine

subroutine kernel_calc_flux(quad,fvrhs)

    use mesh_types_mod, only: quad_t
    use riemann_inner_mod, only: xriemann
    use riemann_inner_mod, only: yriemann

    use kernel_config_mod, only: lofactor
    use equations_mod, only: cons2prim, prim2cons
    use equations_const_mod, only: posvars
    use utils_mod, only: vmin,vmax
    use utils_mod, only: pp
    use math_mod, only: minmod3

    use equations_mod, only: cons2evec

    integer, parameter :: N = N_NODES
    integer, parameter :: V = N_VARS
    
# if 1
    !! third order central reconstruction
    integer, parameter :: Q = 2
    integer, parameter :: R = 3

    !! Gauss quadrature rule
    real(dp), parameter :: xs(Q) = 0.5_dp/SQRT(3.0_dp)*(/-1.0_dp,1.0_dp/)
    real(dp), parameter :: ws(Q) = 0.5_dp*(/1.0_dp,1.0_dp/)

    real(dp), parameter :: xq(Q,R) = transpose(reshape((/&
        & 1.0_dp, xs(1), xs(1)**2, &
        & 1.0_dp, xs(2), xs(2)**2 /),(/R,Q/)))

    real(dp), parameter :: recoMat(R,R) = transpose(reshape((/&
        &   -1.0_dp/24.0_dp, 13.0_dp/12.0_dp, -1.0_dp/24.0_dp, & 
        &   -1.0_dp/ 2.0_dp,          0.0_dp,  1.0_dp/ 2.0_dp, &
        &    1.0_dp/ 2.0_dp,         -1.0_dp,  1.0_dp/ 2.0_dp /), (/R,R/)))

    real(dp), parameter :: xM(R) = (/ 1.0_dp, -1.0/2.0_dp,  1.0_dp/4.0_dp /)
    real(dp), parameter :: xP(R) = (/ 1.0_dp,  1.0/2.0_dp,  1.0_dp/4.0_dp /)

    integer, parameter :: stencil = 1
# endif

# if 0
    !! fifth order central reconstruction
    integer, parameter :: Q = 3
    integer, parameter :: R = 5

    !! Gauss quadrature rule
    real(dp), parameter :: xs(Q) = 0.5_dp*SQRT(3.0_dp/5.0_dp)*(/-1.0_dp,0.0_dp,1.0_dp/)
    real(dp), parameter :: ws(Q) = 0.5_dp*(/5.0_dp/9.0_dp,8.0_dp/9.0_dp,5.0_dp/9.0_dp/)

    real(dp), parameter :: xq(Q,R) = transpose(reshape((/&
        & 1.0_dp, xs(1), xs(1)**2, xs(1)**3, xs(1)**4, &
        & 1.0_dp, xs(2), xs(2)**2, xs(2)**3, xs(2)**4, &
        & 1.0_dp, xs(3), xs(3)**2, xs(3)**3, xs(3)**4   /),(/R,Q/)))

    real(dp), parameter :: recoMat(R,R) = transpose(reshape((/&
        &    3.0_dp/640.0_dp,  -29.0_dp/480.0_dp, 1067.0_dp/ 960.0_dp,  -29.0_dp/480.0_dp,  3.0_dp/640.0_dp, & 
        &    5.0_dp/ 48.0_dp,  -17.0_dp/ 24.0_dp,    0.0_dp,             17.0_dp/ 24.0_dp, -5.0_dp/ 48.0_dp, &
        &   -1.0_dp/ 16.0_dp,    3.0_dp/  4.0_dp,  -11.0_dp/ 8.0_dp,      3.0_dp/  4.0_dp, -1.0_dp/ 16.0_dp, &
        &   -1.0_dp/ 12.0_dp,    1.0_dp/  6.0_dp,   0.0_dp,              -1.0_dp/  6.0_dp,  1.0_dp/ 12.0_dp, &
        &    1.0_dp/ 24.0_dp,   -1.0_dp/  6.0_dp,   1.0_dp /4.0_dp,      -1.0_dp/  6.0_dp,  1.0_dp/ 24.0_dp  &
        & /), (/R,R/)))

    real(dp), parameter :: xM(R) = (/ 1.0_dp, -1.0/2.0_dp,  1.0_dp/4.0_dp, -1.0/8.0_dp, 1.0/16.0_dp /)
    real(dp), parameter :: xP(R) = (/ 1.0_dp,  1.0/2.0_dp,  1.0_dp/4.0_dp,  1.0/8.0_dp, 1.0/16.0_dp /)

    integer, parameter :: stencil = 2
# endif

    real(dp), parameter :: recoTam(R,R) = transpose(recoMat)

    type(quad_t), intent(inout) :: quad
    real(dp), intent(out)       :: fvrhs(N_NODES,N_NODES,N_VARS)

    !! ---------------------------------------------------------------------- !!
    !! inner interfaces

    real(dp) :: halo(-N+1:2*N,-N+1:2*N, N_VARS)

    real(dp) :: reco(R,R,0:N+1,0:N+1,V)

    real(dp) :: xuM(Q,0:N,N,N_VARS)
    real(dp) :: xuP(Q,0:N,N,N_VARS)

    real(dp) :: yuM(Q,N,0:N,N_VARS)
    real(dp) :: yuP(Q,N,0:N,N_VARS)

    !! ---------------------------------------------------------------------- !!
    !! FV operator

    real(dp) :: xrmq(Q, 0:N_NODES, 1:N_NODES, N_VARS)
    real(dp) :: yrmq(Q, 1:N_NODES, 0:N_NODES, N_VARS)

    real(dp) :: xrm(0:N_NODES, 1:N_NODES, N_VARS)
    real(dp) :: yrm(1:N_NODES, 0:N_NODES, N_VARS)

    integer :: i,j,k,h,qq

    real(dp) :: sdx,sdy

    sdx = REAL(N_NODES,dp)/quad%delta(1) 
    sdy = REAL(N_NODES,dp)/quad%delta(2) 

    !! ---------------------------------------------------------------------- !!
    !! preparing inner boundaries

    !! Transform to primitive variables.
    do h = 1,V; do j = -N+1,2*N; do i = -N+1,2*N;
        halo(i,j,h) = quad%aux%hydro%halo(h,i,j)
    end do; end do; end do;

    !! Reconstruct second-degree polynomials.
    do h = 1,V; do j = 0,N+1; do i = 0,N+1;
        reco(:,:,i,j,h) = matmul(recoMat,matmul(halo(i-stencil:i+stencil,j-stencil:j+stencil,h),recoTam))
    end do; end do; end do;

    !! Interpolate to interfaces in x-direction.
    do h = 1,V; do j = 1,N; do i = 0,N;
        xuP(:,i,j,h) = matmul(xq,matmul(xP,reco(:,:,i  ,j,h)))
        xuM(:,i,j,h) = matmul(xq,matmul(xM,reco(:,:,i+1,j,h)))
    end do; end do; end do;

    !! Interpolate to interfaces in y-direction.
    do h = 1,V; do j = 0,N; do i = 1,N;
        yuP(:,i,j,h) = matmul(xq,matmul(reco(:,:,i,j  ,h),xP))
        yuM(:,i,j,h) = matmul(xq,matmul(reco(:,:,i,j+1,h),xM))
    end do; end do; end do;

!! # if 1
!!     !! Transform to primitive variables.
!!     do j = -N+1,2*N; do i = -N+1,2*N;
!!         halo(:,i,j) = cons2prim(quad%aux%hydro%halo(:,i,j))
!!     end do; end do;
!! 
!!     !! 2nd order monotonized central limiter.
!!     do j = 1,N; do i = 0,N; do h = 1,V;
!!         xuP(h,i,j) = halo(h,i ,j) + 0.25_dp*minmod3(&
!!             4.0_dp*(halo(h,i  ,j) - halo(h,i-1,j)), &
!!                    (halo(h,i+1,j) - halo(h,i-1,j)), &
!!             4.0_dp*(halo(h,i+1,j) - halo(h,i  ,j)))
!! 
!!         xuM(h,i,j) = halo(h,i+1,j) - 0.25_dp*minmod3(&
!!             4.0_dp*(halo(h,i+1,j) - halo(h,i  ,j)), &
!!                    (halo(h,i+2,j) - halo(h,i  ,j)), &
!!             4.0_dp*(halo(h,i+2,j) - halo(h,i+1,j)))
!!     end do; end do; end do;
!! 
!!     do j = 0,N; do i = 1,N; do h = 1,V;
!!         yuP(h,i,j) = halo(h,i, j) + 0.25_dp*minmod3(&
!!             4.0_dp*(halo(h,i,j  ) - halo(h,i,j-1)), &
!!                    (halo(h,i,j+1) - halo(h,i,j-1)), &
!!             4.0_dp*(halo(h,i,j+1) - halo(h,i,j  )))
!! 
!!         yuM(h,i,j) = halo(h,i,j+1) - 0.25_dp*minmod3(&
!!             4.0_dp*(halo(h,i,j+1) - halo(h,i,j  )), &
!!                    (halo(h,i,j+2) - halo(h,i,j  )), &
!!             4.0_dp*(halo(h,i,j+2) - halo(h,i,j+1)))
!!     end do; end do; end do;
!! 
!!     do j = 1,N; do i = 0,N;
!!         xuP(:,i,j) = prim2cons(xuP(:,i,j))
!!         xuM(:,i,j) = prim2cons(xuM(:,i,j))
!!     end do; end do;
!! 
!!     do j = 0,N; do i = 1,N;
!!         yuP(:,i,j) = prim2cons(yuP(:,i,j))
!!         yuM(:,i,j) = prim2cons(yuM(:,i,j))
!!     end do; end do;
!! # else
!!     !! Transform to primitive variables.
!!     do j = -N+1,2*N; do i = -N+1,2*N;
!!         halo(:,i,j) = quad%aux%hydro%halo(:,i,j)
!!     end do; end do;
!! 
!!     !! First order FV scheme
!!     do j = 1,N; do i = 0,N;
!!         xuP(:,i,j) = halo(:,i  ,j)
!!         xuM(:,i,j) = halo(:,i+1,j)
!!     end do; end do;
!! 
!!     do j = 0,N; do i = 1,N;
!!         yuP(:,i,j) = halo(:,i,j  )
!!         yuM(:,i,j) = halo(:,i,j+1)
!!     end do; end do;
!! # endif

    do qq = 1,Q; do j = 1,N; do i = 0,N;
        xrmq(qq,i,j,:) = xriemann(xuP(qq,i,j,:),xuM(qq,i,j,:))
    end do; end do; end do;

    do qq = 1,Q; do j = 0,N; do i = 1,N;
        yrmq(qq,i,j,:) = yriemann(yuP(qq,i,j,:),yuM(qq,i,j,:))
    end do; end do; end do;

    xrm = 0.0_dp
    yrm = 0.0_dp

    do j = 1,N; do i = 0,N; do qq = 1,Q; 
        xrm(i,j,:) = xrm(i,j,:) + ws(qq)*xrmq(qq,i,j,:)
    end do; end do; end do;

    do j = 0,N; do i = 1,N; do qq = 1,Q; 
        yrm(i,j,:) = yrm(i,j,:) + ws(qq)*yrmq(qq,i,j,:)
    end do; end do; end do;

    do h = 1,V; do j = 1,N; do i = 1,N;
        fvrhs(i,j,h) = sdx*(xrm(i-1,j,h)-xrm(i,j,h)) + sdy*(yrm(i,j-1,h)-yrm(i,j,h))
    end do; end do; end do;

end subroutine

end module
