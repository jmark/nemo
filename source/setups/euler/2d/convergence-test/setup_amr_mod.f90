module setup_amr_mod

use constants_mod

contains

subroutine hook_pre_timestep_300

    use share_mod, only: mesh

    use mesh_utils_mod, only: mesh_adapt
    use mesh_utils_mod, only: mesh_partition

    return

    call mesh_adapt(mesh,setup_amr_probe_refine,setup_amr_replace_refine,setup_amr_probe_coarsen,setup_amr_replace_coarsen)
    call mesh_partition(mesh)

end subroutine

function setup_amr_probe_refine(quad) result(ok)

    use setup_config_mod
    use equations_const_mod
    use mesh_types_mod, only: quad_t

    type(quad_t), intent(in)    :: quad
    logical                     :: ok

    ok = .false.

    ! if (any(quad%pld%hydro%state(:,:,:,DENS_VAR) > dens_refine_threshold)) then
    !     ok = .true.
    ! else
    !     ok = .false.
    ! end if

end function

function setup_amr_probe_coarsen(family) result(ok)

    use setup_config_mod
    use equations_const_mod
    use mesh_types_mod, only: quadptr_t
    use mesh_types_mod, only: quad_t
    use mesh_const_mod, only: MESH_CHILDREN

    !type(quadptr_t), intent(in) :: family(MESH_CHILDREN)
    type(quad_t), intent(in) :: family(MESH_CHILDREN)
    logical                     :: ok

    integer :: i

    ok = .false.

    ! ok = .true.

    ! do i = 1,MESH_CHILDREN
    !     ! write (*,'((a12),(99(i7)))') 'probe_coarsen', i, family(i)%p%locid,family(i)%p%mesh%nquads
    !     !if (any(family(i)%p%pld%hydro%state(:,:,:,DENS_VAR) > dens_coarse_threshold)) then
    !     if (any(family(i)%pld%hydro%state(:,:,:,DENS_VAR) > dens_coarse_threshold)) then
    !         ok = .false.
    !         return
    !     end if
    ! end do

end function

subroutine setup_amr_replace_refine(incoming,outgoing)

    use mesh_types_mod, only: quad_t,quadptr_t
    use mesh_const_mod, only: MESH_CHILDREN
    use kernel_share_mod, only: refineMat
    use kernel_const_mod, only: DGFV_QMAX
    use tensor_mod, only: blockmulXYZ
    use equations_mod, only: isvalid

    type(quad_t), intent(in)    :: incoming
    type(quad_t), intent(inout) :: outgoing(MESH_CHILDREN)

    integer, parameter :: N = N_NODES
    integer, parameter :: P = N_NODES+1
    integer, parameter :: M = 2*N_NODES

    real(dp) :: temp(2*N_NODES,2*N_NODES,2*N_NODES,N_VARS)

    integer :: q,h,f
    logical :: ok

# if 0
    do h = 1,N_VARS
        call scale_up(incoming%pld%hydro%state(:,:,:,h),temp(:,:,:,h))
    end do

    forall (h = 1:N_VARS)
        outgoing(1)%pld%hydro%state(:,:,:,h) = temp(1:N,1:N,1:N,h)
        outgoing(2)%pld%hydro%state(:,:,:,h) = temp(P:M,1:N,1:N,h)
        outgoing(3)%pld%hydro%state(:,:,:,h) = temp(1:N,P:M,1:N,h)
        outgoing(4)%pld%hydro%state(:,:,:,h) = temp(P:M,P:M,1:N,h)

        outgoing(5)%pld%hydro%state(:,:,:,h) = temp(1:N,1:N,P:M,h)
        outgoing(6)%pld%hydro%state(:,:,:,h) = temp(P:M,1:N,P:M,h)
        outgoing(7)%pld%hydro%state(:,:,:,h) = temp(1:N,P:M,P:M,h)
        outgoing(8)%pld%hydro%state(:,:,:,h) = temp(P:M,P:M,P:M,h)
    end forall

# else

    ! do q = DGFV_QMAX,0,-1
    !     do f = 1,MESH_CHILDREN
    !         do h = 1,N_VARS
    !             outgoing(f)%pld%hydro%state(:,:,:,h) = &
    !                 blockmulXYZ(refineMat(:,:,f,q,X_DIR),refineMat(:,:,f,q,Y_DIR),refineMat(:,:,f,q,Z_DIR),&
    !                     incoming%pld%hydro%state(:,:,:,h))
    !         end do
    !         ok = isvalid(outgoing(f)%pld%hydro%state)
    !         if (.not.ok) exit
    !     end do
    !     if (ok) exit
    ! end do

# endif

end subroutine

subroutine setup_amr_replace_coarsen(incoming,outgoing)

    use mesh_types_mod, only: quad_t,quadptr_t
    use mesh_const_mod, only: MESH_CHILDREN

    use kernel_share_mod, only: coarseMat
    use kernel_const_mod, only: DGFV_QMAX
    use tensor_mod, only: blockmulXYZ
    use equations_mod, only: isvalid

    type(quad_t), intent(in)    :: incoming(MESH_CHILDREN)
    type(quad_t), intent(inout) :: outgoing

    integer, parameter :: N = N_NODES
    integer, parameter :: P = N_NODES+1
    integer, parameter :: M = 2*N_NODES

    real(dp) :: temp(2*N_NODES,2*N_NODES,2*N_NODES,N_VARS)

    integer :: q,h,f

# if 0
    forall (h = 1:N_VARS)
        temp(1:N,1:N,1:N,h) = incoming(1)%pld%hydro%state(:,:,:,h)
        temp(P:M,1:N,1:N,h) = incoming(2)%pld%hydro%state(:,:,:,h)
        temp(1:N,P:M,1:N,h) = incoming(3)%pld%hydro%state(:,:,:,h)
        temp(P:M,P:M,1:N,h) = incoming(4)%pld%hydro%state(:,:,:,h)

        temp(1:N,1:N,P:M,h) = incoming(5)%pld%hydro%state(:,:,:,h)
        temp(P:M,1:N,P:M,h) = incoming(6)%pld%hydro%state(:,:,:,h)
        temp(1:N,P:M,P:M,h) = incoming(7)%pld%hydro%state(:,:,:,h)
        temp(P:M,P:M,P:M,h) = incoming(8)%pld%hydro%state(:,:,:,h)
    end forall

    do h = 1,N_VARS
        call scale_down(temp(:,:,:,h),outgoing%pld%hydro%state(:,:,:,h))
    end do

# else

    ! outgoing%pld%hydro%state = 0.0_dp
    ! do q = DGFV_QMAX,0,-1
    !     do f = 1,MESH_CHILDREN
    !         do h = 1,N_VARS
    !             outgoing%pld%hydro%state(:,:,:,h) = outgoing%pld%hydro%state(:,:,:,h) + &
    !                 blockmulXYZ(coarseMat(:,:,f,q,X_DIR),coarseMat(:,:,f,q,Y_DIR),coarseMat(:,:,3,f,q,Z_DIR),&
    !                     incoming(f)%pld%hydro%state(:,:,:,h))
    !         end do
    !     end do
    !     if (isvalid(outgoing%pld%hydro%state)) exit
    ! end do
# endif

end subroutine

pure subroutine scale_up(bar,foo)
    
    !! R^(N x N) -> R^(2*N x 2*N)

    integer, parameter :: N =   N_NODES
    integer, parameter :: M = 2*N_NODES

    real(dp), intent(in)    :: bar(N,N,N)
    real(dp), intent(out)   :: foo(M,M,M)

    integer :: i,j,k, ii,jj,kk

    do k = 1,N; do j = 1,N; do i = 1,N
        do kk = 0,1; do jj = 0,1; do ii = 0,1
            foo(2*i-ii,2*j-jj,2*k-kk) = bar(i,j,k)
        end do; end do; end do
    end do; end do; end do

end subroutine

pure subroutine scale_down(foo,bar)

    !! R^(2*N x 2*N) -> R^(N x N)

    integer, parameter :: N =   N_NODES
    integer, parameter :: M = 2*N_NODES

    real(dp), intent(in)    :: foo(M,M,M)
    real(dp), intent(out)   :: bar(N,N,N)

    integer :: i,j,k, ii,jj,kk

    bar   = 0.0

    do k = 1,N; do j = 1,N; do i = 1,N
        do kk = 0,1; do jj = 0,1; do ii = 0,1
            bar(i,j,k) = bar(i,j,k) + 0.125*foo(2*i-ii,2*j-jj,2*k-kk)
        end do; end do; end do
    end do; end do; end do

end subroutine

end module
