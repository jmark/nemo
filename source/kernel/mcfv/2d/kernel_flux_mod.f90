module kernel_flux_mod

use constants_mod

use mesh_types_mod, only: quad_t

!! ------------------------------------------------------------------------- !!

contains

subroutine kernel_fill_facevars(quad)

    use mesh_types_mod, only: quad_t

    type(quad_t), intent(inout) :: quad

    integer :: i,h

    integer, parameter :: M = N_NODES-1
    integer, parameter :: N = N_NODES

    associate(fvstate => quad%pld%hydro%state, fvars => quad%aux%hydro%facevars)
    do h = 1,N_VARS; do i = 1,N_NODES;
        !! Convention: increaseing indices denotes
        !! further distance from interface.
        !! 2,1, [inner values], 1,2

        fvars(i,h,2,SOU) = fvstate(M,i,h)
        fvars(i,h,1,SOU) = fvstate(N,i,h)
                   
        fvars(i,h,1,NOR) = fvstate(1,i,h)
        fvars(i,h,2,NOR) = fvstate(2,i,h)
                   
        fvars(i,h,2,EAS) = fvstate(i,M,h)
        fvars(i,h,1,EAS) = fvstate(i,N,h)
                   
        fvars(i,h,1,WES) = fvstate(i,1,h)
        fvars(i,h,2,WES) = fvstate(i,2,h)
    end do; end do;
    end associate

end subroutine

subroutine kernel_fill_faceflux(side)

    use mesh_types_mod, only: quad_t, side_t
    use mesh_const_mod, only: MESH_HALF, MESH_FACES

    use riemann_inner_mod, only: xriemann_inner => xriemann
    use riemann_inner_mod, only: yriemann_inner => yriemann

# if !PP_MESH_PERIODIC
    use riemann_north_mod, only: riemann_north => xriemann
    use riemann_south_mod, only: riemann_south => xriemann

    use riemann_west_mod, only: riemann_west => yriemann
    use riemann_east_mod, only: riemann_east => yriemann
# endif

    use setup_boundary_mod, only: setup_side_boundary

    use kernel_config_mod, only: lofactor
    use kernel_share_mod, only: rfaceMat
    use kernel_share_mod, only: cfaceMat

    use equations_mod, only: cons2prim, prim2cons
    use equations_const_mod, only: posvars
    use utils_mod, only: vmin,vmax
    use math_mod, only: minmod3

# if PP_NON_CONSERVATIVE_FLUX
    use Non_Conservative_flux_mod, only: Non_Conservative_two_point_xFlux
    use Non_Conservative_flux_mod, only: Non_Conservative_two_point_yFlux
# endif

    integer, parameter :: N = N_NODES
    integer, parameter :: V = N_VARS

    type(side_t), intent(inout) :: side

    real(dp) :: recovars(N_NODES,N_VARS,2,MESH_HALF,2)
    real(dp) :: facevars(N_NODES,N_VARS,MESH_HALF,2)
    real(dp) :: flexvars(N_NODES,N_VARS,MESH_HALF,2)
    real(dp) :: fluxvars(N_NODES,N_VARS,MESH_HALF,2)

    real(dp) :: ocervars(N_NODES,N_VARS,2,MESH_HALF,2)

    integer :: i,j,h,f,p,s,o, iside, iquad, q, qq, nquads

    !! Reconstruct face data for involved orders for conforming and non-conforming faces.
    do iside = 1,2
        !! non-conforming face
        if (side%nquads(iside) > 0 .and. side%nquads(iside) < side%nquads(3-iside)) then
            do iquad = 1,MESH_HALF
                do j = 1,2; do h = 1,N_VARS;
                    recovars(:,h,j,iquad,iside) &
                        = matmul(rfaceMat(:,:,iquad),side%quads(1,iside)%p%aux%hydro%facevars(:,h,j,side%faceid(iside)))
                end do; end do;
                recovars(:,:,2,iquad,iside) = 2.0_dp/3.0_dp*recovars(:,:,1,iquad,iside) + 1.0_dp/3.0_dp*recovars(:,:,2,iquad,iside)
            end do

        !! conforming face
        else
            do iquad = 1,side%nquads(iside)
                recovars(:,:,:,iquad,iside) = side%quads(iquad,iside)%p%aux%hydro%facevars(:,:,:,side%faceid(iside))
            end do
        end if
    end do

    nquads = maxval(side%nquads)

    !! Mirror inner states at boundaries.
    if (side%boundary > 0) then
        recovars(:,:,:,:,side%outer) = recovars(:,:,:,:,side%inner)
    end if

    do iside = 1,2
        do iquad = 1,side%nquads(iside)
            side%quads(iquad,iside)%p%aux%hydro%halovars(:,:,side%faceid(iside)) = recovars(:,:,1,iquad,3-iside)
        end do
    end do

    ocervars = recovars

# if PP_KERNEL_DGFV_MCFV_RECO
    !! Transform to primitive variables.
    do iside = 1,2; do iquad = 1,nquads; do j = 1,2; do i = 1,N_NODES;
        recovars(i,:,j,iquad,iside) = cons2prim(recovars(i,:,j,iquad,iside))
    end do; end do; end do; end do;

    !! Reconstruct via 2nd order monotonized minmod limiter.
    do iquad = 1,nquads
        do h = 1,N_VARS
            facevars(:,h,iquad,side%plus) = recovars(:,h,1,iquad,side%plus) + 0.25_dp*minmod3(&
                4.0_dp*(recovars(:,h,1,iquad,side%plus )-recovars(:,h,2,iquad,side%plus)), &
                       (recovars(:,h,1,iquad,side%minus)-recovars(:,h,2,iquad,side%plus)), &
                4.0_dp*(recovars(:,h,1,iquad,side%minus)-recovars(:,h,1,iquad,side%plus)))
        end do
    end do

    do iquad = 1,nquads
        do h = 1,N_VARS
            facevars(:,h,iquad,side%minus) = recovars(:,h,1,iquad,side%minus) - 0.25_dp*minmod3(&
                4.0_dp*(recovars(:,h,1,iquad,side%minus)-recovars(:,h,1,iquad,side%plus)), &
                       (recovars(:,h,2,iquad,side%minus)-recovars(:,h,1,iquad,side%plus)), &
                4.0_dp*(recovars(:,h,2,iquad,side%minus)-recovars(:,h,1,iquad,side%minus)))
        end do
    end do

    !! Reconstruct via 2nd order monotonized minmod limiter.
    do iside = 1,2; do iquad = 1,nquads; do h = 1,size(posvars);
        facevars(:,posvars(h),iquad,iside) = vmax(lofactor*recovars(:,posvars(h),1,iquad,iside),facevars(:,posvars(h),iquad,iside))
    end do; end do; end do;

    !! Transform to conservative variables.
    do iside = 1,2; do iquad = 1,nquads; do i = 1,N_NODES;
        facevars(i,:,iquad,iside) = prim2cons(facevars(i,:,iquad,iside))
    end do; end do; end do;
# else
    facevars = recovars(:,:,1,:,:)
    facevars = recovars(:,:,1,:,:)
# endif

    !! Process boundary conditions.
    if (side%boundary > 0) then
        call setup_side_boundary(side,facevars(:,:,1,side%inner),facevars(:,:,1,side%outer))
    end if

    !! Calculate face fluxes.
    if (side%boundary > 0) then

# if !PP_MESH_PERIODIC
        select case(side%boundary)
            case (NOR)
                do iquad = 1,nquads; do i = 1,N_NODES;
                    flexvars(i,:,iquad,1) = riemann_north(facevars(i,:,iquad,side%plus),facevars(i,:,iquad,side%minus))
                end do; end do;

            case (SOU)
                do iquad = 1,nquads; do i = 1,N_NODES;
                    flexvars(i,:,iquad,1) = riemann_south(facevars(i,:,iquad,side%plus),facevars(i,:,iquad,side%minus))
                end do; end do;

            case (WES)
                do iquad = 1,nquads; do i = 1,N_NODES;
                    flexvars(i,:,iquad,1) = riemann_west(facevars(i,:,iquad,side%plus),facevars(i,:,iquad,side%minus))
                end do; end do;

            case (EAS)
                do iquad = 1,nquads; do i = 1,N_NODES;
                    flexvars(i,:,iquad,1) = riemann_east(facevars(i,:,iquad,side%plus),facevars(i,:,iquad,side%minus))
                end do; end do;
        end select
# endif

    else

        select case(side%direction)
            case (X_DIR)
                do iquad = 1,nquads; do i = 1,N_NODES;
                    flexvars(i,:,iquad,1) = xriemann_inner(facevars(i,:,iquad,side%plus),facevars(i,:,iquad,side%minus))
                end do; end do;

            case (Y_DIR)
                do iquad = 1,nquads; do i = 1,N_NODES;
                    flexvars(i,:,iquad,1) = yriemann_inner(facevars(i,:,iquad,side%plus),facevars(i,:,iquad,side%minus))
                end do; end do;
        end select

    end if

    flexvars(:,:,:,2) = flexvars(:,:,:,1)

# if PP_NON_CONSERVATIVE_FLUX
    select case(side%direction)
        case (X_DIR)
            do iside = 1,2; do iquad = 1,nquads; do i = 1,N_NODES;
                flexvars(i,:,iquad,iside) = flexvars(i,:,iquad,iside) + Non_Conservative_two_point_xFlux(&
                    ocervars(i,:,1,iquad,iside),facevars(i,:,iquad,side%plus),facevars(i,:,iquad,side%minus))
            end do; end do; end do;

        case (Y_DIR)
            do iside = 1,2; do iquad = 1,nquads; do i = 1,N_NODES;
                flexvars(i,:,iquad,iside) = flexvars(i,:,iquad,iside) + Non_Conservative_two_point_yFlux(&
                    ocervars(i,:,1,iquad,iside),facevars(i,:,iquad,side%plus),facevars(i,:,iquad,side%minus))
            end do; end do; end do;
    end select
# endif

    !! Project fluxes back to block faces.
    do iside = 1,2
        if (side%nquads(iside) > 0 .and. side%nquads(iside) < side%nquads(3-iside)) then
            do iquad = 1,MESH_HALF
                do h = 1,N_VARS;
                    fluxvars(:,h,iquad,iside) = matmul(cfaceMat(:,:,iquad),flexvars(:,h,iquad,iside))
                end do
            end do
                
            do iquad = 2,MESH_HALF
                fluxvars(:,:,1,iside) = fluxvars(:,:,1,iside) + fluxvars(:,:,iquad,iside)
            end do
        else
            do iquad = 1,side%nquads(iside)
                fluxvars(:,:,iquad,iside) = flexvars(:,:,iquad,iside)
            end do
        end if
    end do

    !! Copy facefluxes to respective quads.
    do iside = 1,2
        do iquad = 1,side%nquads(iside)
            side%quads(iquad,iside)%p%aux%hydro%faceflux(:,:,side%faceid(iside)) = fluxvars(:,:,iquad,iside)
        end do
    end do

end subroutine

subroutine kernel_calc_flux(quad,fvrhs)

    use mesh_types_mod, only: quad_t
    use riemann_inner_mod, only: xriemann
    use riemann_inner_mod, only: yriemann

    use kernel_config_mod, only: lofactor
    use equations_mod, only: cons2prim, prim2cons
    use equations_const_mod, only: posvars
    use utils_mod, only: vmin,vmax
    use math_mod, only: minmod3

    use equations_mod, only: cons2evec

# if PP_NON_CONSERVATIVE_FLUX
    use equations_share_mod, only: GLM_alpha
    use equations_const_mod, only: GLMP_VAR

    use Non_Conservative_flux_mod, only: Non_Conservative_two_point_xFlux
    use Non_Conservative_flux_mod, only: Non_Conservative_two_point_yFlux
# endif

    integer, parameter :: R = N_NODES-1
    integer, parameter :: N = N_NODES
    integer, parameter :: P = N_NODES+1
    integer, parameter :: V = N_VARS

    type(quad_t), intent(inout) :: quad
    real(dp), intent(out)       :: fvrhs(N_NODES,N_NODES,N_VARS)

    !! ---------------------------------------------------------------------- !!
    !! inner interfaces

    real(dp) :: xlob(0:N+1,N,size(posvars))
    real(dp) :: ylob(N,0:N+1,size(posvars))

    real(dp) :: xrec(0:N+1,N,N_VARS)
    real(dp) :: yrec(N,0:N+1,N_VARS)

    real(dp) :: xuM(R,N,N_VARS)
    real(dp) :: xuP(R,N,N_VARS)

    real(dp) :: yuM(N,R,N_VARS)
    real(dp) :: yuP(N,R,N_VARS)

    !! ---------------------------------------------------------------------- !!
    !! FV operator

    real(dp) :: xrm(0:N_NODES, 1:N_NODES, N_VARS)
    real(dp) :: yrm(1:N_NODES, 0:N_NODES, N_VARS)

    real(dp) :: xeM(R,N)
    real(dp) :: xeP(R,N)

    real(dp) :: yeM(N,R)
    real(dp) :: yeP(N,R)

    real(dp) :: fvepr(N_NODES,N_NODES)

    integer :: i,j,k,h

    real(dp) :: sdx,sdy

    sdx = REAL(N_NODES,dp)/quad%delta(1) 
    sdy = REAL(N_NODES,dp)/quad%delta(2) 

    !! ---------------------------------------------------------------------- !!
    !! preparing inner boundaries

# if PP_KERNEL_DGFV_MCFV_RECO
    do h = 1,V; do j = 1,N; do i = 1,N;
        xrec(i,j,h) = quad%pld%hydro%state(i,j,h)
    end do; end do; end do;

    do h = 1,V; do j = 1,N; do i = 1,N;
        yrec(i,j,h) = quad%pld%hydro%state(i,j,h)
    end do; end do; end do;

    xrec(0,:,:) = quad%aux%hydro%halovars(:,:,NOR)
    xrec(P,:,:) = quad%aux%hydro%halovars(:,:,SOU)
    yrec(:,0,:) = quad%aux%hydro%halovars(:,:,WES)
    yrec(:,P,:) = quad%aux%hydro%halovars(:,:,EAS)

    !! Transform to primitive variables.
    do j = 1,N; do i = 0,P;
        xrec(i,j,:) = cons2prim(xrec(i,j,:))
    end do; end do;

    do j = 0,P; do i = 1,N;
        yrec(i,j,:) = cons2prim(yrec(i,j,:))
    end do; end do;

    !! Reconstruct via 2nd order monotonized minmod limiter.
    do h = 1,V; do j = 1,N; do i = 1,R;
        xuP(i,j,h) = xrec(i,j,h) + 0.25_dp*minmod3(&
            4.0_dp*(xrec(i  ,j,h) - xrec(i-1,j,h)), &
                   (xrec(i+1,j,h) - xrec(i-1,j,h)), &
            4.0_dp*(xrec(i+1,j,h) - xrec(i  ,j,h)))

        xuM(i,j,h) = xrec(i+1,j,h) - 0.25_dp*minmod3(&
            4.0_dp*(xrec(i+1,j,h) - xrec(i  ,j,h)), &
                   (xrec(i+2,j,h) - xrec(i  ,j,h)), &
            4.0_dp*(xrec(i+2,j,h) - xrec(i+1,j,h)))
    end do; end do; end do;

    do h = 1,V; do j = 1,R; do i = 1,N;
        yuP(i,j,h) = yrec(i,j,h) + 0.25_dp*minmod3(&
            4.0_dp*(yrec(i,j  ,h) - yrec(i,j-1,h)), &
                   (yrec(i,j+1,h) - yrec(i,j-1,h)), &
            4.0_dp*(yrec(i,j+1,h) - yrec(i,j  ,h)))

        yuM(i,j,h) = yrec(i,j+1,h) - 0.25_dp*minmod3(&
            4.0_dp*(yrec(i,j+1,h) - yrec(i,j  ,h)), &
                   (yrec(i,j+2,h) - yrec(i,j  ,h)), &
            4.0_dp*(yrec(i,j+2,h) - yrec(i,j+1,h)))
    end do; end do; end do;

    do h = 1,size(posvars)
        xlob(:,:,h) = lofactor*xrec(:,:,posvars(h))
        ylob(:,:,h) = lofactor*yrec(:,:,posvars(h))
    end do

    do h = 1,size(posvars); do j = 1,N; do i = 1,R;
        xuP(i,j,posvars(h)) = vmax(xlob(i,  j,h),xuP(i,j,posvars(h)))
        xuM(i,j,posvars(h)) = vmax(xlob(i+1,j,h),xuM(i,j,posvars(h)))
    end do; end do; end do;

    do h = 1,size(posvars); do j = 1,R; do i = 1,N;
        yuP(i,j,posvars(h)) = vmax(ylob(i,j  ,h),yuP(i,j,posvars(h)))
        yuM(i,j,posvars(h)) = vmax(ylob(i,j+1,h),yuM(i,j,posvars(h)))
    end do; end do; end do;

    do j = 1,N; do i = 1,R;
        xuP(i,j,:) = prim2cons(xuP(i,j,:))
        xuM(i,j,:) = prim2cons(xuM(i,j,:))
    end do; end do;

    do j = 1,R; do i = 1,N;
        yuP(i,j,:) = prim2cons(yuP(i,j,:))
        yuM(i,j,:) = prim2cons(yuM(i,j,:))
    end do; end do;
# else
    !! First order FV scheme
    do h = 1,N_VARS; do j = 1,N; do i = 1,R;
        xuP(i,j,h) = quad%pld%hydro%state(i  ,j  ,h)
        xuM(i,j,h) = quad%pld%hydro%state(i+1,j  ,h)
    end do; end do; end do;

    do h = 1,N_VARS; do j = 1,R; do i = 1,N;
        yuP(i,j,h) = quad%pld%hydro%state(i  ,j  ,h)
        yuM(i,j,h) = quad%pld%hydro%state(i  ,j+1,h)
    end do; end do; end do;
# endif

    do j = 1,N; do i = 1,N-1;
        xrm(i,j,:) = xriemann(xuP(i,j,:),xuM(i,j,:))
    end do; end do;

    do j = 1,N-1; do i = 1,N;
        yrm(i,j,:) = yriemann(yuP(i,j,:),yuM(i,j,:))
    end do; end do;

    xrm(0,:,:) = quad%aux%hydro%faceflux(:,:,NOR)
    xrm(N,:,:) = quad%aux%hydro%faceflux(:,:,SOU)

    yrm(:,0,:) = quad%aux%hydro%faceflux(:,:,WES)
    yrm(:,N,:) = quad%aux%hydro%faceflux(:,:,EAS)

    ! do j = 1,N; do i = 1,N-1;
    !     xeM(i,j) = dot_product(cons2evec(xuP(i,j,:)),xrm(i,j,:))
    !     xeP(i,j) = dot_product(cons2evec(xuM(i,j,:)),xrm(i,j,:))
    ! end do; end do;

    ! do j = 1,N-1; do i = 1,N;
    !     yeM(i,j) = dot_product(cons2evec(yuP(i,j,:)),yrm(i,j,:))
    !     yeP(i,j) = dot_product(cons2evec(yuM(i,j,:)),yrm(i,j,:))
    ! end do; end do;

    ! block
    ! real, parameter :: tol = 1e-8
    ! if (any(xeM > tol)) then
    !     write (*,'(a10,99(ES12.4))') 'xeM', xeM
    ! endif

    ! if (any(xeP > tol)) then
    !     write (*,'(a10,99(ES12.4))') 'xeP', xeP
    ! endif

    ! if (any(yeM > tol)) then
    !     write (*,'(a10,99(ES12.4))') 'yeM', yeM
    ! endif

    ! if (any(yeP > tol)) then
    !     write (*,'(a10,99(ES12.4))') 'yeP', yeP
    ! endif
    ! end block

    do h = 1,V; do j = 1,N; do i = 1,N;
        fvrhs(i,j,h) = sdx*(xrm(i-1,j,h)-xrm(i,j,h)) + sdy*(yrm(i,j-1,h)-yrm(i,j,h))
    end do; end do; end do;

    ! do j = 1,N; do i = 1,N;
    !     fvepr(i,j) = dot_product(cons2evec(quad%pld%hydro%state(i,j,:)),fvrhs(i,j,:))
    ! end do; end do;

    ! quad%aux%hydro%entropy = fvepr
    ! if (any(fvepr > 1e-8)) then
    !     write (*,'(8(ES12.4))') fvepr
    !     write (*,*)
    ! endif

# if PP_NON_CONSERVATIVE_FLUX
    associate(u => quad%pld%hydro%state)
    do j = 1,N; do i = 1,R;
        fvrhs(i  ,j,:) = fvrhs(i  ,j,:) - sdx*Non_Conservative_two_point_xFlux(u(i  ,j,:),xuP(i,j,:),xuM(i,j,:))
        fvrhs(i+1,j,:) = fvrhs(i+1,j,:) + sdx*Non_Conservative_two_point_xFlux(u(i+1,j,:),xuP(i,j,:),xuM(i,j,:))
    end do; end do;

    do j = 1,R; do i = 1,N;
        fvrhs(i,j  ,:) = fvrhs(i,j  ,:) - sdy*Non_Conservative_two_point_yFlux(u(i,j  ,:),yuP(i,j,:),yuM(i,j,:))
        fvrhs(i,j+1,:) = fvrhs(i,j+1,:) + sdy*Non_Conservative_two_point_yFlux(u(i,j+1,:),yuP(i,j,:),yuM(i,j,:))
    end do; end do;
    end associate

    !! GLM damping term.
    fvrhs(:,:,GLMP_VAR) = fvrhs(:,:,GLMP_VAR) - GLM_alpha*quad%pld%hydro%state(:,:,GLMP_VAR)
# endif

end subroutine

end module
