module kernel_flux_mod

use constants_mod

contains

subroutine kernel_fill_facevars_cb(quad)

    use mesh_types_mod, only: quad_t
    use kernel_share_mod, only: lagrItplVecP
    use kernel_share_mod, only: lagrItplVecM

    type(quad_t), intent(inout) :: quad

    integer :: h

    !! ------------------------------------------------------------------ !!

    do h = 1,N_VARS
        quad%aux%hydro%facevars(:,h,NOR) = matmul(lagrItplVecM,quad%pld%hydro%state(:,:,h))
        quad%aux%hydro%facevars(:,h,SOU) = matmul(lagrItplVecP,quad%pld%hydro%state(:,:,h))
        quad%aux%hydro%facevars(:,h,WES) = matmul(quad%pld%hydro%state(:,:,h),lagrItplVecM)
        quad%aux%hydro%facevars(:,h,EAS) = matmul(quad%pld%hydro%state(:,:,h),lagrItplVecP)
    end do

end subroutine

subroutine kernel_fill_faceflux_cb(side)

    use mesh_types_mod, only: quad_t, side_t
    use mesh_const_mod, only: MESH_HALF, MESH_FACES

    use kernel_share_mod, only: rfaceMat
    use kernel_share_mod, only: cfaceMat

    use riemann_inner_mod, only: xriemann_inner => xriemann
    use riemann_inner_mod, only: yriemann_inner => yriemann

# if !PP_MESH_PERIODIC
    use riemann_north_mod, only: riemann_north => xriemann
    use riemann_south_mod, only: riemann_south => xriemann

    use riemann_west_mod, only: riemann_west => yriemann
    use riemann_east_mod, only: riemann_east => yriemann
# endif

    use setup_boundary_mod, only: setup_side_boundary

    use utils_mod, only: pp

    integer, parameter :: N = N_NODES
    integer, parameter :: V = N_VARS

    type(side_t), intent(inout) :: side

    real(dp) :: facevars(N_NODES,N_VARS,MESH_HALF,2)
    real(dp) :: flexvars(N_NODES,N_VARS,MESH_HALF)
    real(dp) :: fluxvars(N_NODES,N_VARS,MESH_HALF,2)

    integer :: i,j,h,f,p,s,o, iside, iquad, nquads

    facevars = 0.0_dp

    nquads = maxval(side%nquads)

    do iside = 1,2
        if (side%nquads(iside) == 0) cycle

        !! non-conforming face
        if (side%nquads(iside) < side%nquads(3-iside)) then
            do iquad = 1,MESH_HALF
                do h = 1,N_VARS
                    facevars(:,h,iquad,iside) = matmul(rfaceMat(:,:,iquad),side%quads(1,iside)%p%aux%hydro%facevars(:,h,side%faceid(iside)))
                end do
            end do

        else !! conforming face
            do iquad = 1,side%nquads(iside)
                facevars(:,:,iquad,iside) = side%quads(iquad,iside)%p%aux%hydro%facevars(:,:,side%faceid(iside))
            end do
        end if
    end do

# if !PP_MESH_PERIODIC
    !! Process boundary conditions.
    if (side%boundary > 0) then
        call setup_side_boundary(side,facevars(:,:,1,side%inner),facevars(:,:,1,side%outer))

        select case(side%boundary)
            case (NOR)
                do iquad = 1,nquads; do i = 1,N_NODES;
                    flexvars(i,:,iquad) = riemann_north(facevars(i,:,iquad,side%plus),facevars(i,:,iquad,side%minus))
                end do; end do;

            case (SOU)
                do iquad = 1,nquads; do i = 1,N_NODES;
                    flexvars(i,:,iquad) = riemann_south(facevars(i,:,iquad,side%plus),facevars(i,:,iquad,side%minus))
                end do; end do;

            case (WES)
                do iquad = 1,nquads; do i = 1,N_NODES;
                    flexvars(i,:,iquad) = riemann_west(facevars(i,:,iquad,side%plus),facevars(i,:,iquad,side%minus))
                end do; end do;

            case (EAS)
                do iquad = 1,nquads; do i = 1,N_NODES;
                    flexvars(i,:,iquad) = riemann_east(facevars(i,:,iquad,side%plus),facevars(i,:,iquad,side%minus))
                end do; end do;
        end select

    else
# endif

        select case(side%direction)
            case (X_DIR)
                do iquad = 1,nquads; do i = 1,N_NODES;
                    flexvars(i,:,iquad) = xriemann_inner(facevars(i,:,iquad,side%plus),facevars(i,:,iquad,side%minus))
                end do; end do;

            case (Y_DIR)
                do iquad = 1,nquads; do i = 1,N_NODES;
                    flexvars(i,:,iquad) = yriemann_inner(facevars(i,:,iquad,side%plus),facevars(i,:,iquad,side%minus))
                end do; end do;
        end select

# if !PP_MESH_PERIODIC
    end if
# endif

    fluxvars = 0.0_dp

    !! Project from nodal values to mean values for each side.
    do iside = 1,2
        if (side%nquads(iside) == 0) cycle

        if (side%nquads(iside) < side%nquads(3-iside)) then
            do iquad = 1,MESH_HALF
                do h = 1,N_VARS
                    fluxvars(:,h,1,iside) = fluxvars(:,h,1,iside) + matmul(cfaceMat(:,:,iquad),flexvars(:,h,iquad))
                end do
            end do
        else
            do iquad = 1,side%nquads(iside)
                fluxvars(:,:,iquad,iside) = flexvars(:,:,iquad)
            end do
        end if
    end do

    !! Copy facefluxes to respective quads
    do iside = 1,2
        do iquad = 1,side%nquads(iside)
            side%quads(iquad,iside)%p%aux%hydro%faceflux(:,:,side%faceid(iside)) = fluxvars(:,:,iquad,iside)
        end do
    end do

# if PP_NON_CONSERVATIVE_FLUX
    block
    real(dp) :: tempavgs(N_NODES,N_VARS,MESH_HALF)
    real(dp) :: faceavgs(N_NODES,N_VARS,MESH_HALF,2)

    tempavgs = 0.5_dp*(facevars(:,:,:,1) + facevars(:,:,:,2))
    
    faceavgs = 0.0_dp

    !! Project faceaverages back to block faces.
    do iside = 1,2
        if (side%nquads(iside) == 0) cycle

        if (side%nquads(iside) < side%nquads(3-iside)) then
            do iquad = 1,MESH_HALF
                do h = 1,N_VARS;
                    faceavgs(:,h,iquad,iside) = matmul(cfaceMat(:,:,iquad),tempavgs(:,h,iquad))
                end do
            end do
                
            do iquad = 2,MESH_HALF
                faceavgs(:,:,1,iside) = faceavgs(:,:,1,iside) + faceavgs(:,:,iquad,iside)
            end do
        else
            do iquad = 1,side%nquads(iside)
                faceavgs(:,:,iquad,iside) = tempavgs(:,:,iquad)
            end do
        end if
    end do

    !! Copy faceavgs to respective quads.
    do iside = 1,2
        do iquad = 1,side%nquads(iside)
            side%quads(iquad,iside)%p%aux%hydro%faceavgs(:,:,side%faceid(iside)) = faceavgs(:,:,iquad,iside)
        end do
    end do
    end block
# endif /* PP_NON_CONSERVATIVE_FLUX */

end subroutine

subroutine kernel_calc_flux(quad,rhs)

    use mesh_types_mod, only: quad_t

    use equations_mod, only: xflux
    use equations_mod, only: yflux

    use kernel_share_mod, only: dgDiffMat, dgDiffTam
    use kernel_share_mod, only: dgsurfMat, dgsurfTam

# if PP_NON_CONSERVATIVE_FLUX
    use setup_config_mod, only: mu0
    use equations_share_mod, only: GLM_alpha

    use Non_Conservative_flux_mod, only: Non_Conservative_two_point_xFlux
    use Non_Conservative_flux_mod, only: Non_Conservative_two_point_yFlux
# endif

    use utils_mod, only: pp

    integer, parameter :: R = N_NODES-1
    integer, parameter :: N = N_NODES
    integer, parameter :: P = N_NODES+1
    integer, parameter :: V = N_VARS

    type(quad_t), intent(inout) :: quad
    real(dp), intent(out)       :: rhs(N_NODES,N_NODES,N_VARS)

    real(dp) :: xrm(1:N_NODES,1:N_NODES,N_VARS)
    real(dp) :: yrm(1:N_NODES,1:N_NODES,N_VARS)

    !! ---------------------------------------------------------------------- !!

    real(dp) :: xfx(N_NODES,N_NODES,N_VARS)
    real(dp) :: yfy(N_NODES,N_NODES,N_VARS)

    real(dp) :: sfl(N_NODES,N_NODES)
    real(dp) :: vfl(N_NODES,N_NODES)

    integer :: i,j,k,h,q, ii, m, s

    real(dp) :: sdx,sdy

    !! ---------------------------------------------------------------------- !!

    sdx = REAL(N_NODES,dp)/quad%delta(1) 
    sdy = REAL(N_NODES,dp)/quad%delta(2) 

    !! ---------------------------------------------------------------------- !!

    xrm = 0.0_dp
    yrm = 0.0_dp

    xrm(1,:,:) = quad%aux%hydro%faceflux(:,:,NOR)
    xrm(N,:,:) = quad%aux%hydro%faceflux(:,:,SOU)

    yrm(:,1,:) = quad%aux%hydro%faceflux(:,:,WES)
    yrm(:,N,:) = quad%aux%hydro%faceflux(:,:,EAS)

    xrm = sdx*xrm
    yrm = sdy*yrm

    ! call pp('xrm', xrm)
    ! call pp('yrm', xrm)

    !! ---------------------------------------------------------------------- !!

    do j = 1,N_NODES; do i = 1,N_NODES;
        xfx(i,j,:) = sdx*xflux(quad%pld%hydro%state(i,j,:))
        yfy(i,j,:) = sdy*yflux(quad%pld%hydro%state(i,j,:))
    end do; end do;

    ! call pp('xfx', xfx)
    ! call pp('yfy', yfy)

    do h = 1,N_VARS
        sfl = matmul(dgSurfMat,xrm(:,:,h)) + matmul(yrm(:,:,h),dgSurfTam)
        vfl = matmul(dgDiffMat,xfx(:,:,h)) + matmul(yfy(:,:,h),dgDiffTam)

        rhs(:,:,h) = sfl + vfl
    end do

!! # if PP_NON_CONSERVATIVE_FLUX
!!     block
!!     real(dp) :: dgvel(N_NODES,N_NODES,3)
!!     real(dp) ::  divB(N_NODES,N_NODES)
!!     real(dp) ::  divG(N_NODES,N_NODES)
!!     real(dp) ::    vB(N_NODES,N_NODES)
!! 
!!     dgvel(:,:,X_DIR) = dgstate(:,:,MOMX_VAR,q)/dgstate(:,:,DENS_VAR,q)
!!     dgvel(:,:,Y_DIR) = dgstate(:,:,MOMY_VAR,q)/dgstate(:,:,DENS_VAR,q)
!!     dgvel(:,:,Z_DIR) = dgstate(:,:,MOMZ_VAR,q)/dgstate(:,:,DENS_VAR,q)
!! 
!!     divB = sdx*matmul(dgdiffMat(:,:,q),dgstate(:,:,MAGX_VAR,q)) + sdy*matmul(dgstate(:,:,MAGY_VAR,q),dgdiffTam(:,:,q))
!!     divB = divB + matmul(surfMat(:,:,q),xrms(:,:,MAGX_VAR,q)) + matmul(yrms(:,:,MAGY_VAR,q),surfTam(:,:,q))
!! 
!!     divG = sdx*matmul(dgdiffMat(:,:,q),dgstate(:,:,GLMP_VAR,q))*dgvel(:,:,X_DIR) &
!!          + sdy*matmul(dgstate(:,:,GLMP_VAR,q),dgdiffTam(:,:,q))*dgvel(:,:,Y_DIR)
!! 
!!     divG = divG + matmul(surfMat(:,:,q),xrms(:,:,GLMP_VAR,q))*dgvel(:,:,X_DIR) &
!!          + matmul(yrms(:,:,GLMP_VAR,q),surfTam(:,:,q))*dgvel(:,:,Y_DIR)
!! 
!!     vB = dgstate(:,:,MAGX_VAR,q)*dgvel(:,:,X_DIR) + dgstate(:,:,MAGY_VAR,q)*dgvel(:,:,Y_DIR)
!! 
!!     !! x-momentum term
!!     vfl = divB*dgstate(:,:,MAGX_VAR,q)
!!     vfl = matmul(dg2fvMat(:,:,q),matmul(vfl,dg2fvTam(:,:,q)))
!!     tmprhs(:,:,MOMX_VAR) = tmprhs(:,:,MOMX_VAR) + vfl
!! 
!!     !! y-momentum term
!!     vfl = divB*dgstate(:,:,MAGY_VAR,q)
!!     vfl = matmul(dg2fvMat(:,:,q),matmul(vfl,dg2fvTam(:,:,q)))
!!     tmprhs(:,:,MOMY_VAR) = tmprhs(:,:,MOMY_VAR) + vfl
!! 
!!     !! z-momentum term
!!     vfl = divB*dgstate(:,:,MAGZ_VAR,q)
!!     vfl = matmul(dg2fvMat(:,:,q),matmul(vfl,dg2fvTam(:,:,q)))
!!     tmprhs(:,:,MOMZ_VAR) = tmprhs(:,:,MOMZ_VAR) + vfl
!! 
!!     !! energy term
!!     vfl = divB*vB/mu0
!!     vfl = matmul(dg2fvMat(:,:,q),matmul(vfl,dg2fvTam(:,:,q)))
!!     tmprhs(:,:,ENER_VAR) = tmprhs(:,:,ENER_VAR) + vfl
!! 
!!     !! x-mag term
!!     vfl = divB*dgvel(:,:,X_DIR)
!!     vfl = matmul(dg2fvMat(:,:,q),matmul(vfl,dg2fvTam(:,:,q)))
!!     tmprhs(:,:,MAGX_VAR) = tmprhs(:,:,MAGX_VAR) + vfl
!! 
!!     !! y-mag term
!!     vfl = divB*dgvel(:,:,Y_DIR)
!!     vfl = matmul(dg2fvMat(:,:,q),matmul(vfl,dg2fvTam(:,:,q)))
!!     tmprhs(:,:,MAGY_VAR) = tmprhs(:,:,MAGY_VAR) + vfl
!! 
!!     !! z-mag term
!!     vfl = divB*dgvel(:,:,Z_DIR)
!!     vfl = matmul(dg2fvMat(:,:,q),matmul(vfl,dg2fvTam(:,:,q)))
!!     tmprhs(:,:,MAGZ_VAR) = tmprhs(:,:,MAGZ_VAR) + vfl
!! 
!!     !! energy term
!!     vfl = divG*dgstate(:,:,GLMP_VAR,q)/mu0
!!     vfl = matmul(dg2fvMat(:,:,q),matmul(vfl,dg2fvTam(:,:,q)))
!!     tmprhs(:,:,ENER_VAR) = tmprhs(:,:,ENER_VAR) + vfl
!! 
!!     !! glmp term
!!     vfl = matmul(dg2fvMat(:,:,q),matmul(divG,dg2fvTam(:,:,q)))
!!     tmprhs(:,:,GLMP_VAR) = tmprhs(:,:,GLMP_VAR) + vfl
!!     end block
!! 
!!     !! GLM damping term.
!!     dgrhs(:,:,GLMP_VAR) = dgrhs(:,:,GLMP_VAR) - GLM_alpha*quad%pld%hydro%state(:,:,GLMP_VAR)
!! # endif

end subroutine

end module
