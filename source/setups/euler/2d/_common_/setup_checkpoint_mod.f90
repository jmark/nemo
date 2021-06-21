!! # define SINGLE_PRECISION 0

module setup_checkpoint_mod

use constants_mod
use kernel_visualize_mod, only: visualize => kernel_visualize

contains

subroutine hook_checkpoint_add_meta_300

    use share_mod, only: chkpt
    use checkpoint_utils_mod, only: checkpoint_add_unknown_name

# if PP_CHECKPOINT_FLASH
    call checkpoint_add_unknown_name(chkpt,"dens")
    call checkpoint_add_unknown_name(chkpt,"momx")
    call checkpoint_add_unknown_name(chkpt,"momy")
    call checkpoint_add_unknown_name(chkpt,"ener")

# if PP_KERNEL_DGFV
    call checkpoint_add_unknown_name(chkpt,"blnd")
# endif

    ! call dset_setattr(chkpt,"/dens",H5T_NATIVE_DOUBLE,"minimum",1.0_dp)
    ! call dset_setattr(chkpt,"/dens",H5T_NATIVE_DOUBLE,"maximum",1.0_dp)

    ! call checkpoint_add_unknown_name(chkpt,"entr")
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
    call dset_create(chkpt,"/ener",H5T_NATIVE_DOUBLE,(/N_NODES,N_NODES,1/))

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

    ! call dset_create(chkpt,"/entr",H5T_NATIVE_DOUBLE,(/N_NODES,N_NODES,1/))

end subroutine

subroutine hook_checkpoint_fill_100

    use share_mod, only: chkpt
    use checkpoint_utils_mod, only: dset_fill => checkpoint_utils_dataset_iterate

    ! call dset_fill(chkpt,"/mesh/mpirank",fill_mpirank_cb)

    call dset_fill(chkpt,"/dens",fill_dens_cb)

    call dset_fill(chkpt,"/momx",fill_momx_cb)
    call dset_fill(chkpt,"/momy",fill_momy_cb)
    call dset_fill(chkpt,"/ener",fill_ener_cb)

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

    ! call dset_fill(chkpt,"/entr",fill_entr_cb)

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

    call dset_write(dset,quad,visualize(quad%pld%hydro%state(:,:,1)))

end subroutine

subroutine fill_momx_cb(dset,quad)

    use mesh_types_mod, only: quad_t
    use checkpoint_types_mod, only: dset_t => checkpoint_dataset_t
    use checkpoint_utils_mod, only: dset_write => checkpoint_utils_dataset_write

    type(dset_t), intent(inout) :: dset
    type(quad_t), intent(inout) :: quad

    call dset_write(dset,quad,visualize(quad%pld%hydro%state(:,:,2)))

end subroutine

subroutine fill_momy_cb(dset,quad)

    use mesh_types_mod, only: quad_t
    use checkpoint_types_mod, only: dset_t => checkpoint_dataset_t
    use checkpoint_utils_mod, only: dset_write => checkpoint_utils_dataset_write

    type(dset_t), intent(inout) :: dset
    type(quad_t), intent(inout) :: quad

    call dset_write(dset,quad,visualize(quad%pld%hydro%state(:,:,3)))

end subroutine

subroutine fill_ener_cb(dset,quad)

    use mesh_types_mod, only: quad_t
    use checkpoint_types_mod, only: dset_t => checkpoint_dataset_t
    use checkpoint_utils_mod, only: dset_write => checkpoint_utils_dataset_write

    type(dset_t), intent(inout) :: dset
    type(quad_t), intent(inout) :: quad
    logical                     :: ok

    call dset_write(dset,quad,visualize(quad%pld%hydro%state(:,:,4)))

end subroutine

! subroutine fill_entr_cb(dset,quad)
! 
!     use mesh_types_mod, only: quad_t
!     use checkpoint_types_mod, only: dset_t => checkpoint_dataset_t
!     use checkpoint_utils_mod, only: dset_write => checkpoint_utils_dataset_write
! 
!     type(dset_t), intent(inout) :: dset
!     type(quad_t), intent(inout) :: quad
!     logical                     :: ok
! 
!     call dset_write(dset,quad,quad%aux%hydro%entropy)
! 
! end subroutine

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

end module
