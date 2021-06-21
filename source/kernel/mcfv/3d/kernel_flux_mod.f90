module kernel_flux_mod

use constants_mod

use mesh_types_mod, only: quad_t

!! ------------------------------------------------------------------------- !!

contains

subroutine kernel_fill_facevars(quad)

    use mesh_types_mod, only: quad_t

    type(quad_t), intent(inout) :: quad

    integer :: i,j,h

    integer, parameter :: M = N_NODES-1
    integer, parameter :: N = N_NODES

    !! Convention: increaseing indices denotes further distance from interface
    associate(fvstate => quad%pld%hydro%state, fvars => quad%aux%hydro%facevars)
    do h = 1,N_VARS; do j = 1,N_NODES; do i = 1,N_NODES;
        fvars(i,j,h,2,SOU) = fvstate(M,i,j,h)
        fvars(i,j,h,1,SOU) = fvstate(N,i,j,h)
                   
        fvars(i,j,h,1,NOR) = fvstate(1,i,j,h)
        fvars(i,j,h,2,NOR) = fvstate(2,i,j,h)
                   
        fvars(i,j,h,2,EAS) = fvstate(i,M,j,h)
        fvars(i,j,h,1,EAS) = fvstate(i,N,j,h)
                   
        fvars(i,j,h,1,WES) = fvstate(i,1,j,h)
        fvars(i,j,h,2,WES) = fvstate(i,2,j,h)

        fvars(i,j,h,2,BAC) = fvstate(i,j,M,h)
        fvars(i,j,h,1,BAC) = fvstate(i,j,N,h)
                   
        fvars(i,j,h,1,FRO) = fvstate(i,j,1,h)
        fvars(i,j,h,2,FRO) = fvstate(i,j,2,h)
    end do; end do; end do;
    end associate

end subroutine

subroutine kernel_fill_faceflux(side)

    use mesh_types_mod, only: quad_t, side_t
    use mesh_const_mod, only: MESH_HALF, MESH_FACES

    use riemann_inner_mod, only: xriemann_inner => xriemann
    use riemann_inner_mod, only: yriemann_inner => yriemann
    use riemann_inner_mod, only: zriemann_inner => zriemann

# if !PP_MESH_PERIODIC
    use riemann_north_mod, only: riemann_north => xriemann
    use riemann_south_mod, only: riemann_south => xriemann

    use riemann_west_mod, only: riemann_west => yriemann
    use riemann_east_mod, only: riemann_east => yriemann

    use riemann_front_mod, only: riemann_front => zriemann
    use riemann_back_mod, only: riemann_back => zriemann
# endif

    use setup_boundary_mod, only: setup_side_boundary

    use kernel_share_mod, only: rfaceMat,rfaceTam
    use kernel_share_mod, only: cfaceMat,cfaceTam

    use equations_mod, only: cons2prim, prim2cons
    use equations_const_mod, only: posvars
    use utils_mod, only: vmax

    use kernel_config_mod, only: lofactor
    use math_mod, only: minmod3

# if PP_NON_CONSERVATIVE_FLUX
    use Non_Conservative_flux_mod, only: Non_Conservative_two_point_xFlux
    use Non_Conservative_flux_mod, only: Non_Conservative_two_point_yFlux
    use Non_Conservative_flux_mod, only: Non_Conservative_two_point_zFlux
# endif

    integer, parameter :: N = N_NODES
    integer, parameter :: V = N_VARS

    type(side_t), intent(inout) :: side

    real(dp) :: recovars(N_NODES,N_NODES,N_VARS,2,MESH_HALF,2)
    real(dp) :: facevars(N_NODES,N_NODES,N_VARS,  MESH_HALF,2)
    real(dp) :: flexvars(N_NODES,N_NODES,N_VARS,  MESH_HALF,2)
    real(dp) :: fluxvars(N_NODES,N_NODES,N_VARS,  MESH_HALF,2)

    real(dp) :: ocervars(N_NODES,N_NODES,N_VARS,2,MESH_HALF,2)

    integer :: i,j,k,h,f, iside, iquad,nquads

    !! Reconstruct face data for involved orders for conforming and non-conforming faces.
    do iside = 1,2
        !! non-conforming face
        if (side%nquads(iside) > 0 .and. side%nquads(iside) < side%nquads(3-iside)) then
            do iquad = 1,MESH_HALF
                do j = 1,2; do h = 1,N_VARS;
                    recovars(:,:,h,j,iquad,iside) &
                        = matmul(rfaceMat(:,:,iquad),&
                          matmul(side%quads(1,iside)%p%aux%hydro%facevars(:,:,h,j,side%faceid(iside)),rfaceTam(:,:,iquad)))
                end do; end do;
                recovars(:,:,:,2,iquad,iside) = 2.0_dp/3.0_dp*recovars(:,:,:,1,iquad,iside) + 1.0_dp/3.0_dp*recovars(:,:,:,2,iquad,iside)
            end do

        !! conforming face
        else
            do iquad = 1,side%nquads(iside)
                recovars(:,:,:,:,iquad,iside) = side%quads(iquad,iside)%p%aux%hydro%facevars(:,:,:,:,side%faceid(iside))
            end do
        end if
    end do

    nquads = maxval(side%nquads)

    !! Mirror inner states at boundaries.
    if (side%boundary > 0) then
        recovars(:,:,:,:,:,side%outer) = recovars(:,:,:,:,:,side%inner)
    end if

    do iside = 1,2
        do iquad = 1,side%nquads(iside)
            side%quads(iquad,iside)%p%aux%hydro%halovars(:,:,:,side%faceid(iside)) = recovars(:,:,:,1,iquad,3-iside)
        end do
    end do

    ocervars = recovars

    !! Transform to primitive variables.
    do iside = 1,2; do iquad = 1,nquads; do k = 1,2; do j = 1,N_NODES; do i = 1,N_NODES;
        recovars(i,j,:,k,iquad,iside) = cons2prim(recovars(i,j,:,k,iquad,iside))
    end do; end do; end do; end do; end do;

    !! Reconstruct via 2nd order monotonized minmod limiter.
    do iquad = 1,nquads
        do h = 1,N_VARS
            facevars(:,:,h,iquad,side%plus) = recovars(:,:,h,1,iquad,side%plus) + 0.25_dp*minmod3(&
                4.0_dp*(recovars(:,:,h,1,iquad,side%plus )-recovars(:,:,h,2,iquad,side%plus)), &
                       (recovars(:,:,h,1,iquad,side%minus)-recovars(:,:,h,2,iquad,side%plus)), &
                4.0_dp*(recovars(:,:,h,1,iquad,side%minus)-recovars(:,:,h,1,iquad,side%plus)))
        end do
    end do

    do iquad = 1,nquads
        do h = 1,N_VARS
            facevars(:,:,h,iquad,side%minus) = recovars(:,:,h,1,iquad,side%minus) - 0.25_dp*minmod3(&
                4.0_dp*(recovars(:,:,h,1,iquad,side%minus)-recovars(:,:,h,1,iquad,side%plus)), &
                       (recovars(:,:,h,2,iquad,side%minus)-recovars(:,:,h,1,iquad,side%plus)), &
                4.0_dp*(recovars(:,:,h,2,iquad,side%minus)-recovars(:,:,h,1,iquad,side%minus)))
        end do
    end do

    !! Reconstruct via 2nd order monotonized minmod limiter.
    do iside = 1,2; do iquad = 1,nquads; do h = 1,size(posvars);
        facevars(:,:,posvars(h),iquad,iside) = vmax(lofactor*recovars(:,:,posvars(h),1,iquad,iside),facevars(:,:,posvars(h),iquad,iside))
    end do; end do; end do;

    !! Transform to conservative variables.
    do iside = 1,2; do iquad = 1,nquads; do j = 1,N_NODES; do i = 1,N_NODES;
        facevars(i,j,:,iquad,iside) = prim2cons(facevars(i,j,:,iquad,iside))
    end do; end do; end do; end do;

    !! Process boundary conditions.
    if (side%boundary > 0) then
        call setup_side_boundary(side,facevars(:,:,:,1,side%inner),facevars(:,:,:,1,side%outer))
    end if

    !! Calculate face fluxes.
    if (side%boundary > 0) then
# if !PP_MESH_PERIODIC
        select case(side%boundary)
            case (NOR)
                do iquad = 1,nquads; do j = 1,N_NODES; do i = 1,N_NODES;
                    flexvars(i,j,:,iquad,1) = riemann_north(facevars(i,j,:,iquad,side%plus),facevars(i,j,:,iquad,side%minus))
                end do; end do; end do;

            case (SOU)
                do iquad = 1,nquads; do j = 1,N_NODES; do i = 1,N_NODES;
                    flexvars(i,j,:,iquad,1) = riemann_south(facevars(i,j,:,iquad,side%plus),facevars(i,j,:,iquad,side%minus))
                end do; end do; end do;

            case (WES)
                do iquad = 1,nquads; do j = 1,N_NODES; do i = 1,N_NODES;
                    flexvars(i,j,:,iquad,1) = riemann_west(facevars(i,j,:,iquad,side%plus),facevars(i,j,:,iquad,side%minus))
                end do; end do; end do;

            case (EAS)
                do iquad = 1,nquads; do j = 1,N_NODES; do i = 1,N_NODES;
                    flexvars(i,j,:,iquad,1) = riemann_east(facevars(i,j,:,iquad,side%plus),facevars(i,j,:,iquad,side%minus))
                end do; end do; end do;

           case (FRO)
                do iquad = 1,nquads; do j = 1,N_NODES; do i = 1,N_NODES;
                    flexvars(i,j,:,iquad,1) = riemann_front(facevars(i,j,:,iquad,side%plus),facevars(i,j,:,iquad,side%minus))
                end do; end do; end do;

           case (BAC)
                do iquad = 1,nquads; do j = 1,N_NODES; do i = 1,N_NODES;
                    flexvars(i,j,:,iquad,1) = riemann_back(facevars(i,j,:,iquad,side%plus),facevars(i,j,:,iquad,side%minus))
                end do; end do; end do;
        end select
# endif
    else
        select case(side%direction)
            case (X_DIR)
                do iquad = 1,nquads; do j = 1,N_NODES; do i = 1,N_NODES;
                    flexvars(i,j,:,iquad,1) = xriemann_inner(facevars(i,j,:,iquad,side%plus),facevars(i,j,:,iquad,side%minus))
                end do; end do; end do;

            case (Y_DIR)
                do iquad = 1,nquads; do j = 1,N_NODES; do i = 1,N_NODES;
                    flexvars(i,j,:,iquad,1) = yriemann_inner(facevars(i,j,:,iquad,side%plus),facevars(i,j,:,iquad,side%minus))
                end do; end do; end do;

            case (Z_DIR)
                do iquad = 1,nquads; do j = 1,N_NODES; do i = 1,N_NODES;
                    flexvars(i,j,:,iquad,1) = zriemann_inner(facevars(i,j,:,iquad,side%plus),facevars(i,j,:,iquad,side%minus))
                end do; end do; end do;
        end select
    end if

    flexvars(:,:,:,:,2) = flexvars(:,:,:,:,1)

# if PP_NON_CONSERVATIVE_FLUX
    select case(side%direction)
        case (X_DIR)
            do iside = 1,2; do iquad = 1,nquads; do j = 1,N_NODES; do i = 1,N_NODES;
                flexvars(i,j,:,iquad,iside) = flexvars(i,j,:,iquad,iside) + Non_Conservative_two_point_xFlux(&
                    ocervars(i,j,:,1,iquad,iside),facevars(i,j,:,iquad,side%plus),facevars(i,j,:,iquad,side%minus))
            end do; end do; end do; end do;

        case (Y_DIR)
            do iside = 1,2; do iquad = 1,nquads; do j = 1,N_NODES; do i = 1,N_NODES;
                flexvars(i,j,:,iquad,iside) = flexvars(i,j,:,iquad,iside) + Non_Conservative_two_point_yFlux(&
                    ocervars(i,j,:,1,iquad,iside),facevars(i,j,:,iquad,side%plus),facevars(i,j,:,iquad,side%minus))
            end do; end do; end do; end do;

        case (Z_DIR)
            do iside = 1,2; do iquad = 1,nquads; do j = 1,N_NODES; do i = 1,N_NODES;
                flexvars(i,j,:,iquad,iside) = flexvars(i,j,:,iquad,iside) + Non_Conservative_two_point_zFlux(&
                    ocervars(i,j,:,1,iquad,iside),facevars(i,j,:,iquad,side%plus),facevars(i,j,:,iquad,side%minus))
            end do; end do; end do; end do;
    end select
# endif

    !! Project fluxes back to block faces.
    do iside = 1,2
        if (side%nquads(iside) > 0 .and. side%nquads(iside) < side%nquads(3-iside)) then
            do iquad = 1,MESH_HALF
                forall (h = 1:N_VARS)
                    fluxvars(:,:,h,iquad,iside) = matmul(cfaceMat(:,:,iquad),matmul(flexvars(:,:,h,iquad,iside),cfaceTam(:,:,iquad)))
                end forall
            end do
                
            do iquad = 2,MESH_HALF
                fluxvars(:,:,:,1,iside) = fluxvars(:,:,:,1,iside) + fluxvars(:,:,:,iquad,iside)
            end do
        else
            do iquad = 1,side%nquads(iside)
                fluxvars(:,:,:,iquad,iside) = flexvars(:,:,:,iquad,iside)
            end do
        end if
    end do

    !! Copy facefluxes to respective quads.
    do iside = 1,2
        do iquad = 1,side%nquads(iside)
            side%quads(iquad,iside)%p%aux%hydro%faceflux(:,:,:,side%faceid(iside)) = fluxvars(:,:,:,iquad,iside)
        end do
    end do

end subroutine

subroutine kernel_calc_flux(quad,fvrhs)

    use mesh_types_mod, only: quad_t
    use riemann_inner_mod, only: xriemann
    use riemann_inner_mod, only: yriemann
    use riemann_inner_mod, only: zriemann

    use equations_mod, only: cons2prim, prim2cons
    use equations_const_mod, only: posvars
    use utils_mod, only: vmax

    use kernel_config_mod, only: lofactor
    use math_mod, only: minmod3

# if PP_NON_CONSERVATIVE_FLUX
    use equations_share_mod, only: GLM_alpha
    use equations_const_mod, only: GLMP_VAR

    use Non_Conservative_flux_mod, only: Non_Conservative_two_point_xFlux
    use Non_Conservative_flux_mod, only: Non_Conservative_two_point_yFlux
    use Non_Conservative_flux_mod, only: Non_Conservative_two_point_zFlux
# endif

    integer, parameter :: R = N_NODES-1
    integer, parameter :: N = N_NODES
    integer, parameter :: P = N_NODES+1
    integer, parameter :: V = N_VARS

    type(quad_t), intent(inout) :: quad
    real(dp), intent(out)       :: fvrhs(N_NODES,N_NODES,N_NODES,N_VARS)

    !! lower thresholds
    real(dp) :: xlob(0:P,N,N,size(posvars))
    real(dp) :: ylob(N,0:P,N,size(posvars))
    real(dp) :: zlob(N,N,0:P,size(posvars))
    
    !! stencils for reconstruction
    real(dp) :: xrec(0:P,N,N,N_VARS)
    real(dp) :: yrec(N,0:P,N,N_VARS)
    real(dp) :: zrec(N,N,0:P,N_VARS)

    !! inner interfaces
    real(dp) :: xuM(R,N,N,N_VARS)
    real(dp) :: xuP(R,N,N,N_VARS)

    real(dp) :: yuM(N,R,N,N_VARS)
    real(dp) :: yuP(N,R,N,N_VARS)

    real(dp) :: zuM(N,N,R,N_VARS)
    real(dp) :: zuP(N,N,R,N_VARS)

    !! interface fluxes
    real(dp) :: xrm(0:N_NODES,1:N_NODES,1:N_NODES,N_VARS)
    real(dp) :: yrm(1:N_NODES,0:N_NODES,1:N_NODES,N_VARS)
    real(dp) :: zrm(1:N_NODES,1:N_NODES,0:N_NODES,N_VARS)

    integer :: i,j,k,h

    real(dp) :: sdx,sdy,sdz

    sdx = REAL(N_NODES,dp)/quad%delta(1) 
    sdy = REAL(N_NODES,dp)/quad%delta(2) 
    sdz = REAL(N_NODES,dp)/quad%delta(3) 

    !! ---------------------------------------------------------------------- !!
    !! Interpolating inner interfaces. Monotonized centrel reconstruction.

    do h = 1,V; do k = 1,N; do j = 1,N; do i = 1,N;
        xrec(i,j,k,h) = quad%pld%hydro%state(i,j,k,h)
        yrec(i,j,k,h) = quad%pld%hydro%state(i,j,k,h)
        zrec(i,j,k,h) = quad%pld%hydro%state(i,j,k,h)
    end do; end do; end do; end do;

    xrec(0,:,:,:) = quad%aux%hydro%halovars(:,:,:,NOR)
    xrec(P,:,:,:) = quad%aux%hydro%halovars(:,:,:,SOU)
    yrec(:,0,:,:) = quad%aux%hydro%halovars(:,:,:,WES)
    yrec(:,P,:,:) = quad%aux%hydro%halovars(:,:,:,EAS)
    zrec(:,:,0,:) = quad%aux%hydro%halovars(:,:,:,FRO)
    zrec(:,:,P,:) = quad%aux%hydro%halovars(:,:,:,BAC)

    !! Transform to primitive variables.
    do k = 1,N; do j = 1,N; do i = 0,P;
        xrec(i,j,k,:) = cons2prim(xrec(i,j,k,:))
    end do; end do; end do;

    do k = 1,N; do j = 0,P; do i = 1,N;
        yrec(i,j,k,:) = cons2prim(yrec(i,j,k,:))
    end do; end do; end do;

    do k = 0,P; do j = 1,N; do i = 1,N;
        zrec(i,j,k,:) = cons2prim(zrec(i,j,k,:))
    end do; end do; end do;

    !! Reconstruct interfaces.
    do h = 1,V; do k = 1,N; do j = 1,N; do i = 1,R;
        xuP(i,j,k,h) = xrec(i  ,j,k,h) + 0.25_dp*minmod3(&
            4.0_dp*(xrec(i  ,j,k,h) - xrec(i-1,j,k,h)), &
                   (xrec(i+1,j,k,h) - xrec(i-1,j,k,h)), &
            4.0_dp*(xrec(i+1,j,k,h) - xrec(i  ,j,k,h)))

        xuM(i,j,k,h) = xrec(i+1,j,k,h) - 0.25_dp*minmod3(&
            4.0_dp*(xrec(i+1,j,k,h) - xrec(i  ,j,k,h)), &
                   (xrec(i+2,j,k,h) - xrec(i  ,j,k,h)), &
            4.0_dp*(xrec(i+2,j,k,h) - xrec(i+1,j,k,h)))
    end do; end do; end do; end do;

    do h = 1,V; do k = 1,N; do j = 1,R; do i = 1,N;
        yuP(i,j,k,h) = yrec(i,j  ,k,h) + 0.25_dp*minmod3(&
            4.0_dp*(yrec(i,j  ,k,h) - yrec(i,j-1,k,h)), &
                   (yrec(i,j+1,k,h) - yrec(i,j-1,k,h)), &
            4.0_dp*(yrec(i,j+1,k,h) - yrec(i,j  ,k,h)))

        yuM(i,j,k,h) = yrec(i,j+1,k,h) - 0.25_dp*minmod3(&
            4.0_dp*(yrec(i,j+1,k,h) - yrec(i,j  ,k,h)), &
                   (yrec(i,j+2,k,h) - yrec(i,j  ,k,h)), &
            4.0_dp*(yrec(i,j+2,k,h) - yrec(i,j+1,k,h)))
    end do; end do; end do; end do;

    do h = 1,V; do k = 1,R; do j = 1,N; do i = 1,N;
        zuP(i,j,k,h) = zrec(i,j,k  ,h) + 0.25_dp*minmod3(&
            4.0_dp*(zrec(i,j,k  ,h) - zrec(i,j,k-1,h)), &
                   (zrec(i,j,k+1,h) - zrec(i,j,k-1,h)), &
            4.0_dp*(zrec(i,j,k+1,h) - zrec(i,j,k  ,h)))

        zuM(i,j,k,h) = zrec(i,j,k+1,h) - 0.25_dp*minmod3(&
            4.0_dp*(zrec(i,j,k+1,h) - zrec(i,j,k  ,h)), &
                   (zrec(i,j,k+2,h) - zrec(i,j,k  ,h)), &
            4.0_dp*(zrec(i,j,k+2,h) - zrec(i,j,k+1,h)))
    end do; end do; end do; end do;

    !! Limit interface states to lower threshold.
    do h = 1,size(posvars);
        xlob(:,:,:,h) = lofactor*xrec(:,:,:,posvars(h))
        ylob(:,:,:,h) = lofactor*yrec(:,:,:,posvars(h))
        zlob(:,:,:,h) = lofactor*zrec(:,:,:,posvars(h))
    end do

    forall (h = 1:size(posvars), k = 1:N, j = 1:N, i = 1:R)
        xuP(i,j,k,posvars(h)) = vmax(xlob(i,  j,k,h),xuP(i,j,k,posvars(h)))
        xuM(i,j,k,posvars(h)) = vmax(xlob(i+1,j,k,h),xuM(i,j,k,posvars(h)))
    end forall

    forall (h = 1:size(posvars), k = 1:N, j = 1:R, i = 1:N)
        yuP(i,j,k,posvars(h)) = vmax(ylob(i,j  ,k,h),yuP(i,j,k,posvars(h)))
        yuM(i,j,k,posvars(h)) = vmax(ylob(i,j+1,k,h),yuM(i,j,k,posvars(h)))
    end forall

    forall (h = 1:size(posvars), k = 1:R, j = 1:N, i = 1:N)
        zuP(i,j,k,posvars(h)) = vmax(zlob(i,j,k  ,h),zuP(i,j,k,posvars(h)))
        zuM(i,j,k,posvars(h)) = vmax(zlob(i,j,k+1,h),zuM(i,j,k,posvars(h)))
    end forall

    !! Convert back to conservative variables.
    forall (k = 1:N, j = 1:N, i = 1:R)
        xuP(i,j,k,:) = prim2cons(xuP(i,j,k,:))
        xuM(i,j,k,:) = prim2cons(xuM(i,j,k,:))
    end forall

    forall (k = 1:N, j = 1:R, i = 1:N)
        yuP(i,j,k,:) = prim2cons(yuP(i,j,k,:))
        yuM(i,j,k,:) = prim2cons(yuM(i,j,k,:))
    end forall

    forall (k = 1:R, j = 1:N, i = 1:N)
        zuP(i,j,k,:) = prim2cons(zuP(i,j,k,:))
        zuM(i,j,k,:) = prim2cons(zuM(i,j,k,:))
    end forall

    !! ---------------------------------------------------------------------- !!
    !! Solve Riemann problems.

    forall (k = 1:N, j = 1:N, i = 1:R) xrm(i,j,k,:) = xriemann(xuP(i,j,k,:),xuM(i,j,k,:))
    forall (k = 1:N, j = 1:R, i = 1:N) yrm(i,j,k,:) = yriemann(yuP(i,j,k,:),yuM(i,j,k,:))
    forall (k = 1:R, j = 1:N, i = 1:N) zrm(i,j,k,:) = zriemann(zuP(i,j,k,:),zuM(i,j,k,:))

    !! ---------------------------------------------------------------------- !!
    !! Append block face fluxes.

    xrm(0,:,:,:) = quad%aux%hydro%faceflux(:,:,:,NOR)
    xrm(N,:,:,:) = quad%aux%hydro%faceflux(:,:,:,SOU)

    yrm(:,0,:,:) = quad%aux%hydro%faceflux(:,:,:,WES)
    yrm(:,N,:,:) = quad%aux%hydro%faceflux(:,:,:,EAS)

    zrm(:,:,0,:) = quad%aux%hydro%faceflux(:,:,:,FRO)
    zrm(:,:,N,:) = quad%aux%hydro%faceflux(:,:,:,BAC)

    !! ---------------------------------------------------------------------- !!
    !! Apply FV operator.

    xrm = sdx*xrm
    yrm = sdy*yrm
    zrm = sdz*zrm

    forall (h = 1:V, k = 1:N, j = 1:N, i = 1:N)
        fvrhs(i,j,k,h) = xrm(i-1,j,k,h)-xrm(i,j,k,h) + yrm(i,j-1,k,h)-yrm(i,j,k,h) + zrm(i,j,k-1,h)-zrm(i,j,k,h)
    end forall

# if PP_NON_CONSERVATIVE_FLUX
    associate(u => quad%pld%hydro%state)
    do k = 1,N; do j = 1,N; do i = 1,R;
        fvrhs(i  ,j,k,:) = fvrhs(i  ,j,k,:) - sdx*Non_Conservative_two_point_xFlux(u(i  ,j,k,:),xuP(i,j,k,:),xuM(i,j,k,:))
        fvrhs(i+1,j,k,:) = fvrhs(i+1,j,k,:) + sdx*Non_Conservative_two_point_xFlux(u(i+1,j,k,:),xuP(i,j,k,:),xuM(i,j,k,:))
    end do; end do; end do;

    do k = 1,N; do j = 1,R; do i = 1,N;
        fvrhs(i,j  ,k,:) = fvrhs(i,j  ,k,:) - sdy*Non_Conservative_two_point_yFlux(u(i,j  ,k,:),yuP(i,j,k,:),yuM(i,j,k,:))
        fvrhs(i,j+1,k,:) = fvrhs(i,j+1,k,:) + sdy*Non_Conservative_two_point_yFlux(u(i,j+1,k,:),yuP(i,j,k,:),yuM(i,j,k,:))
    end do; end do; end do;

    do k = 1,R; do j = 1,N; do i = 1,N;
        fvrhs(i,j,k  ,:) = fvrhs(i,j,k  ,:) - sdz*Non_Conservative_two_point_zFlux(u(i,j,k  ,:),zuP(i,j,k,:),zuM(i,j,k,:))
        fvrhs(i,j,k+1,:) = fvrhs(i,j,k+1,:) + sdz*Non_Conservative_two_point_zFlux(u(i,j,k+1,:),zuP(i,j,k,:),zuM(i,j,k,:))
    end do; end do; end do;
    end associate

    !! GLM damping term.
    fvrhs(:,:,:,GLMP_VAR) = fvrhs(:,:,:,GLMP_VAR) - GLM_alpha*quad%pld%hydro%state(:,:,:,GLMP_VAR)
# endif

end subroutine

end module
