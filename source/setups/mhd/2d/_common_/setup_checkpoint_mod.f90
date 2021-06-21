!! # define SINGLE_PRECISION 0

module setup_checkpoint_mod

use constants_mod

contains

subroutine hook_checkpoint_add_meta_300

    use share_mod, only: chkpt
    use checkpoint_utils_mod, only: checkpoint_add_unknown_name

# if PP_CHECKPOINT_FLASH
    call checkpoint_add_unknown_name(chkpt,"dens")

    call checkpoint_add_unknown_name(chkpt,"momx")
    call checkpoint_add_unknown_name(chkpt,"momy")
    call checkpoint_add_unknown_name(chkpt,"momz")

    call checkpoint_add_unknown_name(chkpt,"ener")

    call checkpoint_add_unknown_name(chkpt,"magx")
    call checkpoint_add_unknown_name(chkpt,"magy")
    call checkpoint_add_unknown_name(chkpt,"magz")

    call checkpoint_add_unknown_name(chkpt,"glmp")
    call checkpoint_add_unknown_name(chkpt,"divb")

# if PP_KERNEL_DGFV
    call checkpoint_add_unknown_name(chkpt,"blnd")
# endif

    ! call dset_setattr(chkpt,"/dens",H5T_NATIVE_DOUBLE,"minimum",1.0_dp)
    ! call dset_setattr(chkpt,"/dens",H5T_NATIVE_DOUBLE,"maximum",1.0_dp)
# endif

end subroutine

subroutine hook_checkpoint_create_100

    use hdf5
    use share_mod, only: chkpt
# if PP_KERNEL_DGFV
    use kernel_const_mod, only: DGFV_QMAX
# endif
    use checkpoint_utils_mod, only: dset_create => checkpoint_utils_dataset_create
    use checkpoint_utils_mod, only: dset_setattr => checkpoint_utils_dataset_set_attribute

    call dset_create(chkpt,"/mesh/mpirank",H5T_NATIVE_INTEGER)

    call dset_create(chkpt,"/dens",H5T_NATIVE_DOUBLE,(/N_NODES,N_NODES,1/))

    call dset_create(chkpt,"/momx",H5T_NATIVE_DOUBLE,(/N_NODES,N_NODES,1/))
    call dset_create(chkpt,"/momy",H5T_NATIVE_DOUBLE,(/N_NODES,N_NODES,1/))
    call dset_create(chkpt,"/momz",H5T_NATIVE_DOUBLE,(/N_NODES,N_NODES,1/))

    call dset_create(chkpt,"/ener",H5T_NATIVE_DOUBLE,(/N_NODES,N_NODES,1/))

    call dset_create(chkpt,"/magx",H5T_NATIVE_DOUBLE,(/N_NODES,N_NODES,1/))
    call dset_create(chkpt,"/magy",H5T_NATIVE_DOUBLE,(/N_NODES,N_NODES,1/))
    call dset_create(chkpt,"/magz",H5T_NATIVE_DOUBLE,(/N_NODES,N_NODES,1/))

    call dset_create(chkpt,"/glmp",H5T_NATIVE_DOUBLE,(/N_NODES,N_NODES,1/))

    call dset_create(chkpt,"/divb",H5T_NATIVE_DOUBLE,(/N_NODES,N_NODES,1/))

# if PP_KERNEL_DGFV
    call dset_create(chkpt,"/blnd",H5T_NATIVE_DOUBLE,(/N_NODES,N_NODES,1/))
    call dset_create(chkpt,"/dgfv/blend",H5T_NATIVE_DOUBLE,(/N_NODES,N_NODES,DGFV_QMAX/))
# endif

# if MODULE_SELF_GRAVITY
    ! call dset_create(chkpt,"/gpot",H5T_NATIVE_DOUBLE,(/N_NODES,N_NODES,N_NODES/))
    ! call dset_create(chkpt,"/accx",H5T_NATIVE_DOUBLE,(/N_NODES,N_NODES,N_NODES/))
    ! call dset_create(chkpt,"/accy",H5T_NATIVE_DOUBLE,(/N_NODES,N_NODES,N_NODES/))
    ! call dset_create(chkpt,"/accz",H5T_NATIVE_DOUBLE,(/N_NODES,N_NODES,N_NODES/))
# endif

    ! call dset_create(chkpt,"/minloc",H5T_NATIVE_INTEGER)
    ! call dset_create(chkpt,"/dtmin",H5T_NATIVE_DOUBLE,(/N_NODES,N_NODES,N_NODES/))

end subroutine

subroutine hook_checkpoint_fill_100

    use share_mod, only: chkpt
    use share_mod, only: mesh
    use checkpoint_utils_mod, only: dset_fill => checkpoint_utils_dataset_iterate

    use mesh_utils_mod, only: mesh_exchange
    use mesh_utils_mod, only: mesh_iterate_sides
    use mesh_utils_mod, only: mesh_iterate_quads

    call mesh_exchange(mesh)
# if PP_KERNEL_DGFV
    call mesh_iterate_quads(mesh,checkpoint_fill_blend_cb)
# endif
    call mesh_iterate_sides(mesh,checkpoint_fill_halo_cb)

    ! call dset_fill(chkpt,"/mesh/mpirank",fill_mpirank_cb)

    call dset_fill(chkpt,"/dens",fill_dens_cb)

    call dset_fill(chkpt,"/momx",fill_momx_cb)
    call dset_fill(chkpt,"/momy",fill_momy_cb)
    call dset_fill(chkpt,"/momz",fill_momz_cb)

    call dset_fill(chkpt,"/ener",fill_ener_cb)

    call dset_fill(chkpt,"/magx",fill_magx_cb)
    call dset_fill(chkpt,"/magy",fill_magy_cb)
    call dset_fill(chkpt,"/magz",fill_magz_cb)

    call dset_fill(chkpt,"/glmp",fill_glmp_cb)

# if PP_KERNEL_DGFV
    call dset_fill(chkpt,"/blnd",fill_blnd_cb)
    call dset_fill(chkpt,"/dgfv/blend",fill_blend_cb)
# endif

# if MODULE_SELF_GRAVITY
    ! call dset_fill(chkpt,"/gpot",fill_potential_cb)
    ! call dset_fill(chkpt,"/accx",fill_accx_cb)
    ! call dset_fill(chkpt,"/accy",fill_accy_cb)
    ! call dset_fill(chkpt,"/accz",fill_accz_cb)
# endif

    ! call dset_fill(chkpt,"/minloc",fill_dtminloc_cb)
    ! call dset_fill(chkpt,"/dtmin",fill_dtmin_cb)

    call dset_fill(chkpt,"/divb",fill_divb_cb)

end subroutine

subroutine fill_mpirank_cb(dset,quad)
    
    use share_mod, only: rt => runtime

    use mesh_types_mod, only: quad_t
    use checkpoint_types_mod, only: dset_t => checkpoint_dataset_t
    use checkpoint_utils_mod, only: dset_write => checkpoint_utils_dataset_write

    type(dset_t), intent(inout) :: dset
    type(quad_t), intent(inout) :: quad

    call dset_write(dset,quad,rt%mpi%rank)

end subroutine

subroutine fill_dens_cb(dset,quad)

    use mesh_types_mod, only: quad_t
    use checkpoint_types_mod, only: dset_t => checkpoint_dataset_t
    use checkpoint_utils_mod, only: dset_write => checkpoint_utils_dataset_write

    type(dset_t), intent(inout) :: dset
    type(quad_t), intent(inout) :: quad

    call dset_write(dset,quad,quad%pld%hydro%state(:,:,1))

end subroutine

subroutine fill_momx_cb(dset,quad)

    use mesh_types_mod, only: quad_t
    use checkpoint_types_mod, only: dset_t => checkpoint_dataset_t
    use checkpoint_utils_mod, only: dset_write => checkpoint_utils_dataset_write

    type(dset_t), intent(inout) :: dset
    type(quad_t), intent(inout) :: quad

    call dset_write(dset,quad,quad%pld%hydro%state(:,:,2))

end subroutine

subroutine fill_momy_cb(dset,quad)

    use mesh_types_mod, only: quad_t
    use checkpoint_types_mod, only: dset_t => checkpoint_dataset_t
    use checkpoint_utils_mod, only: dset_write => checkpoint_utils_dataset_write

    type(dset_t), intent(inout) :: dset
    type(quad_t), intent(inout) :: quad

    call dset_write(dset,quad,quad%pld%hydro%state(:,:,3))

end subroutine

subroutine fill_momz_cb(dset,quad)

    use mesh_types_mod, only: quad_t
    use checkpoint_types_mod, only: dset_t => checkpoint_dataset_t
    use checkpoint_utils_mod, only: dset_write => checkpoint_utils_dataset_write

    type(dset_t), intent(inout) :: dset
    type(quad_t), intent(inout) :: quad

    call dset_write(dset,quad,quad%pld%hydro%state(:,:,4))

end subroutine


subroutine fill_ener_cb(dset,quad)

    use mesh_types_mod, only: quad_t
    use checkpoint_types_mod, only: dset_t => checkpoint_dataset_t
    use checkpoint_utils_mod, only: dset_write => checkpoint_utils_dataset_write

    type(dset_t), intent(inout) :: dset
    type(quad_t), intent(inout) :: quad
    logical                     :: ok

    call dset_write(dset,quad,quad%pld%hydro%state(:,:,5))

end subroutine

subroutine fill_magx_cb(dset,quad)

    use mesh_types_mod, only: quad_t
    use checkpoint_types_mod, only: dset_t => checkpoint_dataset_t
    use checkpoint_utils_mod, only: dset_write => checkpoint_utils_dataset_write

    type(dset_t), intent(inout) :: dset
    type(quad_t), intent(inout) :: quad

    call dset_write(dset,quad,quad%pld%hydro%state(:,:,6))

end subroutine

subroutine fill_magy_cb(dset,quad)

    use mesh_types_mod, only: quad_t
    use checkpoint_types_mod, only: dset_t => checkpoint_dataset_t
    use checkpoint_utils_mod, only: dset_write => checkpoint_utils_dataset_write

    type(dset_t), intent(inout) :: dset
    type(quad_t), intent(inout) :: quad

    call dset_write(dset,quad,quad%pld%hydro%state(:,:,7))

end subroutine

subroutine fill_magz_cb(dset,quad)

    use mesh_types_mod, only: quad_t
    use checkpoint_types_mod, only: dset_t => checkpoint_dataset_t
    use checkpoint_utils_mod, only: dset_write => checkpoint_utils_dataset_write

    type(dset_t), intent(inout) :: dset
    type(quad_t), intent(inout) :: quad

    call dset_write(dset,quad,quad%pld%hydro%state(:,:,8))

end subroutine

subroutine fill_glmp_cb(dset,quad)

    use mesh_types_mod, only: quad_t
    use checkpoint_types_mod, only: dset_t => checkpoint_dataset_t
    use checkpoint_utils_mod, only: dset_write => checkpoint_utils_dataset_write

    type(dset_t), intent(inout) :: dset
    type(quad_t), intent(inout) :: quad
    logical                     :: ok

    call dset_write(dset,quad,quad%pld%hydro%state(:,:,9))

end subroutine

! subroutine fill_dtminloc_cb(dset,quad)
! 
!     use share_mod, only: rt => runtime
!     use mesh_types_mod, only: quad_t
!     use checkpoint_types_mod, only: dset_t => checkpoint_dataset_t
!     use checkpoint_utils_mod, only: dset_write => checkpoint_utils_dataset_write
! 
!     type(dset_t), intent(inout) :: dset
!     type(quad_t), intent(inout) :: quad
!     logical                     :: ok
! 
!     integer :: dtminloc
! 
!     dtminloc = INT(quad%aux%hydro%dtmin_count)
!     !dtminloc = quad%level
!     !dtminloc = rt%mpi%rank
! 
!     call dset_write(dset,quad,dtminloc)
! 
!     quad%aux%hydro%dtmin_count = 0
! 
! end subroutine

! subroutine fill_dtmin_cb(dset,quad)
! 
!     use share_mod, only: rt => runtime
!     use mesh_types_mod, only: quad_t
!     use checkpoint_types_mod, only: dset_t => checkpoint_dataset_t
!     use checkpoint_utils_mod, only: dset_write => checkpoint_utils_dataset_write
! 
!     type(dset_t), intent(inout) :: dset
!     type(quad_t), intent(inout) :: quad
!     logical                     :: ok
! 
!     call dset_write(dset,quad,quad%aux%hydro%dtmin)
! 
! end subroutine

# if PP_KERNEL_DGFV
subroutine fill_blnd_cb(dset,quad)

    use mesh_types_mod, only: quad_t
    use checkpoint_types_mod, only: dset_t => checkpoint_dataset_t
    use checkpoint_utils_mod, only: dset_write => checkpoint_utils_dataset_write

    use kernel_const_mod, only: DGFV_QMAX
    use kernel_utils_mod, only: inflate

    type(dset_t), intent(inout) :: dset
    type(quad_t), intent(inout) :: quad
    logical                     :: ok

    real(dp) :: blend(N_NODES,N_NODES,DGFV_QMAX)
    real(dp) :: weighted(N_NODES,N_NODES)

    integer :: q

    blend = inflate(quad%aux%hydro%blend,1,DGFV_QMAX)
    ! blend = inflate(quad%aux%hydro%cummulative_blend/quad%aux%hydro%blend_count,1,DGFV_QMAX)

    weighted = 1.0_dp
    do q = 1,DGFV_QMAX
        weighted = weighted*(1.0_dp - blend(:,:,q)) + REAL(2**q,dp) * blend(:,:,q)
    end do

    ! weighted = blend(:,:,DGFV_QMAX)

    call dset_write(dset,quad,weighted)

    !quad%aux%hydro%cummulative_blend = 0.0_dp
    !quad%aux%hydro%blend_count = 0.0_dp

end subroutine

subroutine fill_blend_cb(dset,quad)

    use mesh_types_mod, only: quad_t
    use checkpoint_types_mod, only: dset_t => checkpoint_dataset_t
    use checkpoint_utils_mod, only: dset_write => checkpoint_utils_dataset_write

    use kernel_const_mod, only: DGFV_QMAX
    use kernel_utils_mod, only: inflate

    type(dset_t), intent(inout) :: dset
    type(quad_t), intent(inout) :: quad
    logical                     :: ok

    real(dp) :: blend(N_NODES,N_NODES,DGFV_QMAX)

    blend = inflate(quad%aux%hydro%blend,1,DGFV_QMAX)
    ! blend = inflate(quad%aux%hydro%cummulative_blend/quad%aux%hydro%blend_count,1,DGFV_QMAX)

    call dset_write(dset,quad,blend)

    quad%aux%hydro%cummulative_blend = 0.0_dp
    quad%aux%hydro%blend_count = 0.0_dp

end subroutine
# endif

# if PP_KERNEL_DGFV
subroutine fill_divb_cb(dset,quad)

    use equations_const_mod
    use mesh_types_mod, only: quad_t
    use checkpoint_types_mod, only: dset_t => checkpoint_dataset_t
    use checkpoint_utils_mod, only: dset_write => checkpoint_utils_dataset_write

    use kernel_const_mod, only: DGFV_QMAX
    use kernel_share_mod, only: diffMat,diffTam
    use kernel_share_mod, only: dg2fvMat,dg2fvTam
    use kernel_share_mod, only: fv2dgMat,fv2dgTam

    use kernel_utils_mod, only: inflate

    type(dset_t), intent(inout) :: dset
    type(quad_t), intent(inout) :: quad

    real(dp) :: magx(0:N_NODES+1,0:N_NODES+1)
    real(dp) :: magy(0:N_NODES+1,0:N_NODES+1)
    real(dp) :: magz(0:N_NODES+1,0:N_NODES+1)

    real(dp) :: divB(N_NODES,N_NODES)
    real(dp) :: dgdv(N_NODES,N_NODES)

    real(dp) :: dgstate(1:N_NODES,1:N_NODES,N_VARS,DGFV_QMAX)

    integer :: i,j

    real(dp) :: sdx,sdy

    integer, parameter :: N = N_NODES
    integer, parameter :: P = N_NODES+1

    integer :: h,q,qmin,qmax

    real(dp) :: blend(N_NODES,N_NODES,DGFV_QMAX)
    real(dp) :: dnebl(N_NODES,N_NODES,DGFV_QMAX)

    blend = inflate(quad%aux%hydro%blend,1,DGFV_QMAX)
    dnebl = 1.0_dp - blend
 
    qmin = quad%aux%hydro%qmin
    qmax = quad%aux%hydro%qmax

    magx(1:N,1:N) = quad%pld%hydro%state(:,:,MAGX_VAR)
    magy(1:N,1:N) = quad%pld%hydro%state(:,:,MAGY_VAR)
    ! magz(1:N,1:N) = quad%pld%hydro%state(:,:,MAGZ_VAR)

    magx(0,1:N) = quad%aux%analyze%halo(:,MAGX_VAR,NOR)
    magx(P,1:N) = quad%aux%analyze%halo(:,MAGX_VAR,SOU)
    magx(1:N,0) = quad%aux%analyze%halo(:,MAGX_VAR,WES)
    magx(1:N,P) = quad%aux%analyze%halo(:,MAGX_VAR,EAS)

    magy(0,1:N) = quad%aux%analyze%halo(:,MAGY_VAR,NOR)
    magy(P,1:N) = quad%aux%analyze%halo(:,MAGY_VAR,SOU)
    magy(1:N,0) = quad%aux%analyze%halo(:,MAGY_VAR,WES)
    magy(1:N,P) = quad%aux%analyze%halo(:,MAGY_VAR,EAS)

    ! magz(0,1:N) = quad%aux%analyze%halo(:,MAGZ_VAR,NOR)
    ! magz(P,1:N) = quad%aux%analyze%halo(:,MAGZ_VAR,SOU)
    ! magz(1:N,0) = quad%aux%analyze%halo(:,MAGZ_VAR,WES)
    ! magz(1:N,P) = quad%aux%analyze%halo(:,MAGZ_VAR,EAS)

    sdx = N_NODES/quad%delta(1)
    sdy = N_NODES/quad%delta(2)

    if (quad%aux%hydro%subcell) then
        do j = 1,N_NODES; do i = 1,N_NODES;
            divB(i,j) = 0.5_dp*sdx*(magx(i+1,j) - magx(i-1,j)) + 0.5_dp*sdy*(magy(i,j+1) - magy(i,j-1))
        end do; end do;

        do q = qmin,qmax
            do h = MAGX_VAR,MAGY_VAR
                dgstate(:,:,h,q) = matmul(fv2dgMat(:,:,q),matmul(quad%pld%hydro%state(:,:,h),fv2dgTam(:,:,q)))
            end do
            dgdv = sdx*matmul(diffMat(:,:,q),dgstate(:,:,MAGX_VAR,q)) + sdy*matmul(dgstate(:,:,MAGY_VAR,q),diffTam(:,:,q))
            divB = dnebl(:,:,q)*divB + blend(:,:,q)*matmul(dg2fvMat(:,:,q),matmul(dgdv,dg2fvTam(:,:,q)))
        end do
    else
        q = qmax
        do h = MAGX_VAR,MAGY_VAR
            dgstate(:,:,h,q) = matmul(fv2dgMat(:,:,q),matmul(quad%pld%hydro%state(:,:,h),fv2dgTam(:,:,q)))
        end do
        dgdv = sdx*matmul(diffMat(:,:,q),dgstate(:,:,MAGX_VAR,q)) + sdy*matmul(dgstate(:,:,MAGY_VAR,q),diffTam(:,:,q))
        divB = matmul(dg2fvMat(:,:,q),matmul(dgdv,dg2fvTam(:,:,q)))
    end if

    call dset_write(dset,quad,divB)

end subroutine

# else

subroutine fill_divb_cb(dset,quad)

    use equations_const_mod
    use mesh_types_mod, only: quad_t
    use checkpoint_types_mod, only: dset_t => checkpoint_dataset_t
    use checkpoint_utils_mod, only: dset_write => checkpoint_utils_dataset_write

    type(dset_t), intent(inout) :: dset
    type(quad_t), intent(inout) :: quad

    real(dp) :: magx(0:N_NODES+1,0:N_NODES+1)
    real(dp) :: magy(0:N_NODES+1,0:N_NODES+1)
    real(dp) :: magz(0:N_NODES+1,0:N_NODES+1)

    real(dp) :: divB(N_NODES,N_NODES)

    integer :: i,j

    real(dp) :: sdx,sdy

    integer, parameter :: N = N_NODES
    integer, parameter :: P = N_NODES+1

    magx(1:N,1:N) = quad%pld%hydro%state(:,:,MAGX_VAR)
    magy(1:N,1:N) = quad%pld%hydro%state(:,:,MAGY_VAR)
    magz(1:N,1:N) = quad%pld%hydro%state(:,:,MAGZ_VAR)

    magx(0,1:N) = quad%aux%analyze%halo(:,MAGX_VAR,NOR)
    magx(P,1:N) = quad%aux%analyze%halo(:,MAGX_VAR,SOU)
    magx(1:N,0) = quad%aux%analyze%halo(:,MAGX_VAR,WES)
    magx(1:N,P) = quad%aux%analyze%halo(:,MAGX_VAR,EAS)

    magy(0,1:N) = quad%aux%analyze%halo(:,MAGY_VAR,NOR)
    magy(P,1:N) = quad%aux%analyze%halo(:,MAGY_VAR,SOU)
    magy(1:N,0) = quad%aux%analyze%halo(:,MAGY_VAR,WES)
    magy(1:N,P) = quad%aux%analyze%halo(:,MAGY_VAR,EAS)

    magz(0,1:N) = quad%aux%analyze%halo(:,MAGZ_VAR,NOR)
    magz(P,1:N) = quad%aux%analyze%halo(:,MAGZ_VAR,SOU)
    magz(1:N,0) = quad%aux%analyze%halo(:,MAGZ_VAR,WES)
    magz(1:N,P) = quad%aux%analyze%halo(:,MAGZ_VAR,EAS)

    sdx = 0.5_dp*N_NODES/quad%delta(1)
    sdy = 0.5_dp*N_NODES/quad%delta(2)

    do j = 1,N_NODES; do i = 1,N_NODES;
        divB(i,j) = sdx*(magx(i+1,j) - magx(i-1,j)) + sdy*(magy(i,j+1) - magy(i,j-1))
    end do; end do;

    call dset_write(dset,quad,divB)

end subroutine
# endif

subroutine checkpoint_fill_halo_cb(side)

    use mesh_types_mod, only: quad_t, side_t
    use mesh_const_mod, only: MESH_HALF, MESH_FACES

    use setup_boundary_mod, only: setup_side_boundary

    integer, parameter :: N = N_NODES
    integer, parameter :: V = N_VARS

    type(side_t), intent(inout) :: side

    integer :: iside, iquad

    real(dp), parameter :: dummy(N_NODES,N_VARS) = 0.0_dp

    if (side%boundary > 0) then

        call setup_side_boundary(side,dummy,side%quads(1,side%inner)%p%aux%analyze%halo(:,:,side%faceid(side%inner)))

    else

        !! Only support for uniform grid now!
        do iside = 1,2
            select case(side%faceid(3-iside)) !! opposite side
                case (NOR); do iquad = 1,1; side%quads(iquad,iside)%p%aux%analyze%halo(:,:,side%faceid(iside)) = side%quads(iquad,3-iside)%p%pld%hydro%state(1,:,:); end do
                case (SOU); do iquad = 1,1; side%quads(iquad,iside)%p%aux%analyze%halo(:,:,side%faceid(iside)) = side%quads(iquad,3-iside)%p%pld%hydro%state(N,:,:); end do
                case (WES); do iquad = 1,1; side%quads(iquad,iside)%p%aux%analyze%halo(:,:,side%faceid(iside)) = side%quads(iquad,3-iside)%p%pld%hydro%state(:,1,:); end do
                case (EAS); do iquad = 1,1; side%quads(iquad,iside)%p%aux%analyze%halo(:,:,side%faceid(iside)) = side%quads(iquad,3-iside)%p%pld%hydro%state(:,N,:); end do
            end select
        end do
    end if

end subroutine

# if PP_KERNEL_DGFV
subroutine checkpoint_fill_blend_cb(quad)

    use kernel_const_mod, only: DGFV_QMAX
    use mesh_const_mod, only: MESH_FACES

    use mesh_types_mod, only: quad_t

    use kernel_const_mod, only: L_L,L_R
    #:for q in range(1,PP_KERNEL_DGFV_QMAX+1)
    use kernel_const_mod, only: L${q}$L,L${q}$R
    use kernel_limiter_mod, only: kernel_blend_limit_${q}$
    #:endfor

    type(quad_t), intent(inout) :: quad

    integer, parameter :: PROBE_QMIN = 1
    integer, parameter :: PROBE_QMAX = DGFV_QMAX

    integer :: q,qmin,qmax

    real(dp) :: blend(L_L:L_R)
    real(dp) :: limit(L_L:L_R)

    logical :: subcell, allfine, allkaput

    qmin = PROBE_QMAX
    qmax = PROBE_QMAX

    limit = 0.0_dp
    blend = 0.0_dp

    do q = PROBE_QMAX,PROBE_QMIN,-1
        qmin = q

        select case (q)
            #:for q in range(PP_KERNEL_DGFV_QMAX,0,-1)
            case (${q}$)
                call kernel_blend_limit_${q}$(quad%pld%hydro%state,blend(L${q}$L:L${q}$R),limit(L${q}$L:L${q}$R),allkaput,allfine)
            #:endfor
        end select

        if (allfine) exit
        ! if (allkaput) qmax = q-1
    end do

    subcell = (qmin < DGFV_QMAX) .or. (qmin == 1)

    quad%aux%hydro%subcell = subcell

    quad%aux%hydro%qmin = qmin
    quad%aux%hydro%qmax = qmax

    quad%aux%hydro%limit = limit
    quad%aux%hydro%blend = blend

end subroutine
# endif
 
end module
