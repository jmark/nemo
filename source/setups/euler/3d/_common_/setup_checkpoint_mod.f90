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

    call dset_create(chkpt,"/dens",H5T_NATIVE_DOUBLE,(/N_NODES,N_NODES,N_NODES/))
    call dset_create(chkpt,"/momx",H5T_NATIVE_DOUBLE,(/N_NODES,N_NODES,N_NODES/))
    call dset_create(chkpt,"/momy",H5T_NATIVE_DOUBLE,(/N_NODES,N_NODES,N_NODES/))
    call dset_create(chkpt,"/momz",H5T_NATIVE_DOUBLE,(/N_NODES,N_NODES,N_NODES/))
    call dset_create(chkpt,"/ener",H5T_NATIVE_DOUBLE,(/N_NODES,N_NODES,N_NODES/))

# if PP_KERNEL_DGFV
    call dset_create(chkpt,"/blnd",H5T_NATIVE_DOUBLE,(/N_NODES,N_NODES,N_NODES/))
    call dset_create(chkpt,"/dgfv/blend",H5T_NATIVE_DOUBLE,(/N_NODES,N_NODES,N_NODES,DGFV_QMAX/))
# endif

# if MODULE_SELF_GRAVITY
    ! call dset_create(chkpt,"/gpot",H5T_NATIVE_DOUBLE,(/N_NODES,N_NODES,N_NODES/))
    ! call dset_create(chkpt,"/accx",H5T_NATIVE_DOUBLE,(/N_NODES,N_NODES,N_NODES/))
    ! call dset_create(chkpt,"/accy",H5T_NATIVE_DOUBLE,(/N_NODES,N_NODES,N_NODES/))
    ! call dset_create(chkpt,"/accz",H5T_NATIVE_DOUBLE,(/N_NODES,N_NODES,N_NODES/))
# endif

    call dset_create(chkpt,"/minloc",H5T_NATIVE_INTEGER)
    call dset_create(chkpt,"/dtmin",H5T_NATIVE_DOUBLE)

end subroutine

subroutine hook_checkpoint_fill_100

    use share_mod, only: chkpt
    use checkpoint_utils_mod, only: dset_fill => checkpoint_utils_dataset_iterate

    ! call dset_fill(chkpt,"/mesh/mpirank",fill_mpirank_cb)

    call dset_fill(chkpt,"/dens",fill_dens_cb)

    call dset_fill(chkpt,"/momx",fill_momx_cb)
    call dset_fill(chkpt,"/momy",fill_momy_cb)
    call dset_fill(chkpt,"/momz",fill_momz_cb)
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

    call dset_fill(chkpt,"/minloc",fill_dtminloc_cb)
    call dset_fill(chkpt,"/dtmin",fill_dtmin_cb)

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

    call dset_write(dset,quad,quad%pld%hydro%state(:,:,:,1))

end subroutine

subroutine fill_momx_cb(dset,quad)

    use mesh_types_mod, only: quad_t
    use checkpoint_types_mod, only: dset_t => checkpoint_dataset_t
    use checkpoint_utils_mod, only: dset_write => checkpoint_utils_dataset_write

    type(dset_t), intent(inout) :: dset
    type(quad_t), intent(inout) :: quad

    call dset_write(dset,quad,quad%pld%hydro%state(:,:,:,2))

end subroutine

subroutine fill_momy_cb(dset,quad)

    use mesh_types_mod, only: quad_t
    use checkpoint_types_mod, only: dset_t => checkpoint_dataset_t
    use checkpoint_utils_mod, only: dset_write => checkpoint_utils_dataset_write

    type(dset_t), intent(inout) :: dset
    type(quad_t), intent(inout) :: quad

    call dset_write(dset,quad,quad%pld%hydro%state(:,:,:,3))

end subroutine

subroutine fill_momz_cb(dset,quad)

    use mesh_types_mod, only: quad_t
    use checkpoint_types_mod, only: dset_t => checkpoint_dataset_t
    use checkpoint_utils_mod, only: dset_write => checkpoint_utils_dataset_write

    type(dset_t), intent(inout) :: dset
    type(quad_t), intent(inout) :: quad

    call dset_write(dset,quad,quad%pld%hydro%state(:,:,:,4))

end subroutine

subroutine fill_ener_cb(dset,quad)

    use mesh_types_mod, only: quad_t
    use checkpoint_types_mod, only: dset_t => checkpoint_dataset_t
    use checkpoint_utils_mod, only: dset_write => checkpoint_utils_dataset_write

    type(dset_t), intent(inout) :: dset
    type(quad_t), intent(inout) :: quad
    logical                     :: ok

    call dset_write(dset,quad,quad%pld%hydro%state(:,:,:,5))

end subroutine

# if MODULE_SELF_GRAVITY
subroutine fill_potential_cb(dset,quad)

    use mesh_types_mod, only: quad_t
    use checkpoint_types_mod, only: dset_t => checkpoint_dataset_t
    use checkpoint_utils_mod, only: dset_write => checkpoint_utils_dataset_write

    type(dset_t), intent(inout) :: dset
    type(quad_t), intent(inout) :: quad
    logical                     :: ok

    call dset_write(dset,quad,quad%aux%gravity%pot)

end subroutine

subroutine fill_accx_cb(dset,quad)

    use mesh_types_mod, only: quad_t
    use checkpoint_types_mod, only: dset_t => checkpoint_dataset_t
    use checkpoint_utils_mod, only: dset_write => checkpoint_utils_dataset_write

    type(dset_t), intent(inout) :: dset
    type(quad_t), intent(inout) :: quad
    logical                     :: ok

    call dset_write(dset,quad,quad%aux%gravity%acc(:,:,:,1))

end subroutine

subroutine fill_accy_cb(dset,quad)

    use mesh_types_mod, only: quad_t
    use checkpoint_types_mod, only: dset_t => checkpoint_dataset_t
    use checkpoint_utils_mod, only: dset_write => checkpoint_utils_dataset_write

    type(dset_t), intent(inout) :: dset
    type(quad_t), intent(inout) :: quad
    logical                     :: ok

    call dset_write(dset,quad,quad%aux%gravity%acc(:,:,:,2))

end subroutine

subroutine fill_accz_cb(dset,quad)

    use mesh_types_mod, only: quad_t
    use checkpoint_types_mod, only: dset_t => checkpoint_dataset_t
    use checkpoint_utils_mod, only: dset_write => checkpoint_utils_dataset_write

    type(dset_t), intent(inout) :: dset
    type(quad_t), intent(inout) :: quad
    logical                     :: ok

    call dset_write(dset,quad,quad%aux%gravity%acc(:,:,:,3))

end subroutine
# endif

# if 1
subroutine fill_dtminloc_cb(dset,quad)

    use share_mod, only: rt => runtime
    use mesh_types_mod, only: quad_t
    use checkpoint_types_mod, only: dset_t => checkpoint_dataset_t
    use checkpoint_utils_mod, only: dset_write => checkpoint_utils_dataset_write

    type(dset_t), intent(inout) :: dset
    type(quad_t), intent(inout) :: quad
    logical                     :: ok

    integer :: dtminloc

    dtminloc = INT(quad%aux%hydro%dtmin_count)
    !dtminloc = quad%level
    !dtminloc = rt%mpi%rank

    call dset_write(dset,quad,dtminloc)

    quad%aux%hydro%dtmin_count = 0

end subroutine
# endif

# if 1
subroutine fill_dtmin_cb(dset,quad)

    use share_mod, only: rt => runtime
    use mesh_types_mod, only: quad_t
    use checkpoint_types_mod, only: dset_t => checkpoint_dataset_t
    use checkpoint_utils_mod, only: dset_write => checkpoint_utils_dataset_write

    type(dset_t), intent(inout) :: dset
    type(quad_t), intent(inout) :: quad
    logical                     :: ok

    call dset_write(dset,quad,quad%aux%hydro%dtmin)
    !call dset_write(dset,quad,quad%aux%hydro%dtmin/quad%aux%hydro%dtmincount)
    !quad%aux%hydro%dtmincount = 0

end subroutine
# endif

# if 0
subroutine fill_fluxes_cb(dset,quad)

    use mesh_types_mod, only: quad_t
    use checkpoint_types_mod, only: dset_t => checkpoint_dataset_t
    use checkpoint_utils_mod, only: dset_write => checkpoint_utils_dataset_write

    use setup_types_mod, only: auxiliary_t
    use quad_utils_mod, only: quad_get_auxiliary

    type(dset_t), intent(inout) :: dset
    type(quad_t), intent(in)    :: quad
    logical                     :: ok

    type(auxiliary_t), pointer  :: aux
    
    call dset_write(dset,quad,quad%aux%hydro%fluxes)

end subroutine
# endif

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

    real(dp) :: blend(N_NODES,N_NODES,N_NODES,DGFV_QMAX)
    real(dp) :: weighted(N_NODES,N_NODES,N_NODES)

    integer :: q

    blend = inflate(quad%aux%hydro%blend,1,DGFV_QMAX)
    ! blend = inflate(quad%aux%hydro%cummulative_blend/quad%aux%hydro%blend_count,1,DGFV_QMAX)

    weighted = 1.0_dp
    do q = 1,DGFV_QMAX
        weighted = weighted*(1.0_dp - blend(:,:,:,q)) + REAL(2**q,dp) * blend(:,:,:,q)
    end do

    ! weighted = blend(:,:,:,DGFV_QMAX)

    call dset_write(dset,quad,weighted)

end subroutine
# endif

# if PP_KERNEL_DGFV
subroutine fill_blend_cb(dset,quad)

    use mesh_types_mod, only: quad_t
    use checkpoint_types_mod, only: dset_t => checkpoint_dataset_t
    use checkpoint_utils_mod, only: dset_write => checkpoint_utils_dataset_write

    use kernel_const_mod, only: DGFV_QMAX
    use kernel_utils_mod, only: inflate

    type(dset_t), intent(inout) :: dset
    type(quad_t), intent(inout) :: quad
    logical                     :: ok

    real(dp) :: blend(N_NODES,N_NODES,N_NODES,DGFV_QMAX)

    !blend = inflate(quad%aux%hydro%blend,1,DGFV_QMAX)
    blend = inflate(quad%aux%hydro%cummulative_blend/quad%aux%hydro%blend_count,1,DGFV_QMAX)

    call dset_write(dset,quad,blend)

    quad%aux%hydro%cummulative_blend = 0.0_dp
    quad%aux%hydro%blend_count = 0.0_dp

end subroutine
# endif

# if 0
subroutine fill_sensor_cb(dset,quad)

    use mesh_types_mod, only: quad_t
    use checkpoint_types_mod, only: dset_t => checkpoint_dataset_t
    use checkpoint_utils_mod, only: dset_write => checkpoint_utils_dataset_write_float64

    use setup_types_mod, only: auxiliary_t
    use quad_utils_mod, only: quad_get_auxiliary

    type(dset_t), intent(inout) :: dset
    type(quad_t), intent(in)    :: quad
    logical                     :: ok

    type(auxiliary_t), pointer  :: aux
    
    real(dp) :: sensor(N_NODES,N_NODES)

    aux => quad_get_auxiliary(quad)

    sensor = aux%amr%sensor

    call dset_write(dset,quad,sensor)

end subroutine

# if 0
subroutine fill_modal_cb(dset,quad)

    use mesh_types_mod, only: quad_t
    use checkpoint_types_mod, only: dset_t => checkpoint_dataset_t
    use checkpoint_utils_mod, only: dset_write => checkpoint_utils_dataset_write_float64

    use equations_const_mod
    use hydro_const_mod, only: DGFV_QMAX

    use hydro_utils_mod, only: inflate

    use setup_types_mod, only: payload_t
    use setup_types_mod, only: auxiliary_t

    use quad_utils_mod, only: quad_get_payload
    use quad_utils_mod, only: quad_get_auxiliary

    use hydro_share_mod, only: toModal, otModal
    use hydro_share_mod, only: fv2dgMat, fv2dgTam

    #:for q in range(1,PP_KERNEL_DGFV_QMAX+1)
    use hydro_utils_mod, only: glue${q}$
    use hydro_utils_mod, only: split${q}$
    use hydro_utils_mod, only: calc_edof${q}$
    #:endfor

    use equations_mod, only: pressure

    type(dset_t), intent(inout) :: dset
    type(quad_t), intent(in)    :: quad
    logical                     :: ok

    type(payload_t), pointer    :: pld
    type(auxiliary_t), pointer  :: aux
    
    real(dp) :: nodal(N_NODES,N_NODES)
    real(dp) :: modal(N_NODES,N_NODES,DGFV_QMAX)
    real(dp) :: edofv(N_NODES,N_NODES,DGFV_QMAX)

    real(dp) :: edof1(2,2,4,4)
    real(dp) :: edof2(4,4,2,2)
    real(dp) :: edof3(8,8)

    integer :: q, i, j

    pld => quad_get_payload(quad)

    nodal = pressure(pld%hydro%state)

    do q = 1,DGFV_QMAX
        edofv(:,:,q) = matmul(fv2dgMat(:,:,q),matmul(nodal,fv2dgTam(:,:,q)))
    end do

    do q = 1,DGFV_QMAX
        modal(:,:,q) = matmul(toModal(:,:,q),matmul(edofv(:,:,q),otModal(:,:,q)))
    end do

    edof1 = split1(modal(:,:,1))
    edof2 = split2(modal(:,:,2))
    edof3 = modal(:,:,3)

    do i = 1,4; do j = 1,4
        edof1(:,:,i,j) = calc_edof1(edof1(:,:,i,j))
    end do; end do

    do i = 1,2; do j = 1,2
        edof2(:,:,i,j) = calc_edof2(edof2(:,:,i,j))
    end do; end do

    edof3 = calc_edof3(edof3)

    modal(:,:,1) = glue1(edof1)
    modal(:,:,2) = glue2(edof2)
    modal(:,:,3) = edof3 

    ok = dset_write(dset,quad,modal)

end subroutine
# endif
# endif

end module
