module setup_amr_mod

use constants_mod

real(dp), save :: fv8fv3Mat(3,8)
real(dp), save :: fv8fv3Tam(8,3)

contains

subroutine hook_pre_timestep_300

    use share_mod, only: mesh

    use mesh_utils_mod, only: mesh_adapt
    use mesh_utils_mod, only: mesh_partition

    call mesh_adapt(mesh,setup_amr_probe_refine,setup_amr_replace_refine,setup_amr_probe_coarsen,setup_amr_replace_coarsen)
    call mesh_partition(mesh)

end subroutine

subroutine hook_init_223

    use interpol_mod, only: LegendreGaussNodesAndWeights
    use subcellfv_mod, only: subcellfv_vanderMondeMatrix
    use linalg_mod, only: invert
    use utils_mod, only: pp

    real(dp) :: dgnodes(8), dgweights(8)

    real(dp) :: dg2fvMat8(8,8)
    real(dp) :: fv2dgMat8(8,8)

    call LegendreGaussNodesAndWeights(7,dgnodes,dgweights)
    call subcellfv_vanderMondeMatrix(8,8,dgnodes,dgweights,dg2fvMat8)
    fv2dgMat8 = invert(dg2fvMat8,8)

    call subcellfv_vanderMondeMatrix(8,3,dgnodes,dgweights,fv8fv3Mat)

    fv8fv3Mat = matmul(fv8fv3Mat,fv2dgMat8)
    fv8fv3Tam = transpose(fv8fv3Mat)

    ! call pp('fv8fv3Mat', fv8fv3Mat)
    ! call pp('fv8fv3Tam', fv8fv3Tam)

end subroutine

!! pure function calc_loehner(u)
!! 
!!     real(dp), intent(in)    :: u(N_NODES,N_NODES)
!!     real(dp)                :: calc_loehner
!! 
!!     real(dp), parameter :: eps = 0.25_dp
!! 
!!     integer :: i,j
!! 
!!     real(dp) :: xtemp(2:N_NODES-1,1:N_NODES)
!!     real(dp) :: ytemp(1:N_NODES,2:N_NODES-1)
!! 
!!     forall (j = 1:N_NODES, i = 2:N_NODES-1)
!!         xtemp(i,j) = abs(u(i-1,j) - 2.0_dp*u(i,j) + u(i+1,j)) / &
!!             (abs(u(i+1,j)-u(i,j)) + abs(u(i,j)-u(i-1,j)) + &
!!             eps*(abs(u(i+1,j)) + 2.0_dp*abs(u(i,j)) + abs(u(i-1,j))))
!!     end forall
!! 
!!     forall (j = 2:N_NODES-1, i = 1:N_NODES)
!!         ytemp(i,j) = abs(u(i,j-1) - 2.0_dp*u(i,j) + u(i,j+1)) / &
!!             (abs(u(i,j+1)-u(i,j)) + abs(u(i,j)-u(i,j-1)) + &
!!             eps*(abs(u(i,j+1)) + 2.0_dp*abs(u(i,j)) + abs(u(i,j-1))))
!!     end forall
!! 
!!     calc_loehner = max(maxval(xtemp),maxval(ytemp))
!! 
!! end function

pure function calc_loehner(u)

    real(dp), intent(in)    :: u(N_NODES,N_NODES)
    real(dp)                :: calc_loehner

    real(dp), parameter :: eps = 0.2_dp

    real(dp) :: stencil(-1:1,-1:1)
    real(dp) :: xtemp(-1:1)
    real(dp) :: ytemp(-1:1)

    integer :: i,j

    stencil = matmul(fv8fv3Mat,matmul(u,fv8fv3Tam))

    forall (j = -1:1)
        xtemp(j) = abs(stencil(-1,j) - 2.0_dp*stencil(0,j) + stencil(1,j)) / &
                 (abs(stencil(1,j)-stencil(0,j)) + abs(stencil(0,j)-stencil(-1,j)) + &
                 eps*(abs(stencil(1,j)) + 2.0_dp*abs(stencil(0,j)) + abs(stencil(-1,j))))
    end forall

    forall (i = -1:1)
         ytemp(i) = abs(stencil(i,-1) - 2.0_dp*stencil(i,0) + stencil(i,1)) / &
             (abs(stencil(i,1)-stencil(i,0)) + abs(stencil(i,0)-stencil(i,-1)) + &
             eps*(abs(stencil(i,1)) + 2.0_dp*abs(stencil(i,0)) + abs(stencil(i,-1))))
     end forall

    ! calc_loehner = max(maxval(xtemp),maxval(ytemp))
    calc_loehner = max(sum(xtemp)/3,sum(ytemp)/3)

end function


function setup_amr_probe_refine(quad) result(ok)

    use setup_config_mod
    use equations_const_mod
    use mesh_types_mod, only: quad_t
# if PP_KERNEL_DGFV
    use kernel_const_mod, only: DGFV_QMAX
# endif

    type(quad_t), intent(in)    :: quad
    logical                     :: ok

    real(dp) :: loehner

    real(dp), parameter :: reftol = 0.3_dp

# if PP_KERNEL_DGFV
    if (quad%aux%hydro%subcell) then
        ok = .true.
        return
    end if
# endif

    loehner = calc_loehner(quad%pld%hydro%state(:,:,DENS_VAR))

    ok = loehner > reftol

end function

function setup_amr_probe_coarsen(family) result(ok)

    use setup_config_mod
    use equations_const_mod
    use mesh_types_mod, only: quad_t
    use mesh_const_mod, only: MESH_CHILDREN
    use mesh_const_mod, only: MESH_CORNERS

# if PP_KERNEL_DGFV
    use kernel_const_mod, only: DGFV_QMAX
# endif
    use kernel_share_mod, only: coarseMat
    use kernel_share_mod, only: coarseTam
    use share_mod, only: rt => runtime
    use quad_utils_mod, only: quad_get_vertices

    type(quad_t), intent(in)    :: family(MESH_CHILDREN)
    logical                     :: ok

    real(dp) :: loehner
    real(dp), parameter :: cortol = 0.1_dp

    real(dp) :: verts(N_DIMS,MESH_CORNERS)
    real(dp) :: coarsened(N_NODES,N_NODES)

    integer :: q,h,f

    ! ok = .false.
    ! return

    ! do f = 1,MESH_CHILDREN
    !     if (family(f)%aux%hydro%subcell) then
    !         ok = .false.
    !         return
    !     end if
    ! end do

    if (rt%simtime < 0.0025_dp) then
        do f = 1,MESH_CHILDREN
            verts = quad_get_vertices(family(f))
            if (any(abs(SQRT(SUM((verts)**2,dim=1))) < init_amr_min_diameter)) then
                ok = .false.
                return
            end if
        end do
    end if

# if PP_KERNEL_DGFV
    coarsened = 0.0_dp
    do q = DGFV_QMAX,0,-1
        do f = 1,MESH_CHILDREN
            coarsened = coarsened + &
                matmul(coarseMat(:,:,X_DIR,f,q),matmul(family(f)%pld%hydro%state(:,:,DENS_VAR),coarseTam(:,:,Y_DIR,f,q)))
        end do
    end do
# endif

    loehner = calc_loehner(coarsened)

    ok = loehner < cortol

end function

subroutine setup_amr_replace_refine(incoming,outgoing)

    use mesh_types_mod, only: quad_t
    use mesh_const_mod, only: MESH_CHILDREN

    use kernel_share_mod, only: refineMat
    use kernel_share_mod, only: refineTam

# if PP_KERNEL_DGFV
    use kernel_const_mod, only: DGFV_QMAX
# endif
    use equations_mod, only: isvalid

    type(quad_t), intent(in)    :: incoming
    type(quad_t), intent(inout) :: outgoing(MESH_CHILDREN)

    integer :: q,h,f
    logical :: ok

# if PP_KERNEL_DGFV
    !do q = DGFV_QMAX,0,-1
    do q = 0,0,-1
        do f = 1,MESH_CHILDREN
            do h = 1,N_VARS
                outgoing(f)%pld%hydro%state(:,:,h) = &
                    matmul(refineMat(:,:,X_DIR,f,q),matmul(incoming%pld%hydro%state(:,:,h),refineTam(:,:,Y_DIR,f,q)))
            end do
            ok = isvalid(outgoing(f)%pld%hydro%state)
            if (.not.ok) exit
        end do
        if (ok) exit
    end do
# endif

end subroutine

subroutine setup_amr_replace_coarsen(incoming,outgoing)

    use mesh_types_mod, only: quad_t
    use mesh_const_mod, only: MESH_CHILDREN

# if PP_KERNEL_DGFV
    use kernel_const_mod, only: DGFV_QMAX
# endif
    use equations_mod, only: isvalid

    use kernel_share_mod, only: coarseMat
    use kernel_share_mod, only: coarseTam

    type(quad_t), intent(in)    :: incoming(MESH_CHILDREN)
    type(quad_t), intent(inout) :: outgoing

    integer :: q,h,f

# if PP_KERNEL_DGFV
    outgoing%pld%hydro%state = 0.0_dp
    !do q = DGFV_QMAX,0,-1
    do q = 0,0,-1
        do f = 1,MESH_CHILDREN
            do h = 1,N_VARS
                outgoing%pld%hydro%state(:,:,h) = outgoing%pld%hydro%state(:,:,h) + &
                    matmul(coarseMat(:,:,X_DIR,f,q),matmul(incoming(f)%pld%hydro%state(:,:,h),coarseTam(:,:,Y_DIR,f,q)))
            end do
        end do
        if (isvalid(outgoing%pld%hydro%state)) exit
    end do
# endif

end subroutine

end module
