module setup_checkpoint_mod

use constants_mod

contains

subroutine hook_checkpoint_create_100

    use hdf5
    use share_mod, only: chkpt
    use kernel_const_mod, only: DGFV_QMAX
    use checkpoint_utils_mod, only: dset_create => checkpoint_utils_dataset_create

    call dset_create(chkpt,"/mesh/levels",H5T_STD_I8LE)
    call dset_create(chkpt,"/mesh/morton",H5T_STD_I32LE,(/N_DIMS/))

    call dset_create(chkpt,"/mesh/coordinates",H5T_NATIVE_DOUBLE,(/N_DIMS/))
    call dset_create(chkpt,"/mesh/block size",H5T_NATIVE_DOUBLE,(/N_DIMS/))
    !call dset_create(chkpt,"/mesh/bounding box",H5T_NATIVE_DOUBLE,(/2,N_DIMS/))

    call dset_create(chkpt,"/data/states/dens",H5T_NATIVE_DOUBLE,(/N_NODES,N_NODES,N_NODES/))
    call dset_create(chkpt,"/data/states/momx",H5T_NATIVE_DOUBLE,(/N_NODES,N_NODES,N_NODES/))
    call dset_create(chkpt,"/data/states/momy",H5T_NATIVE_DOUBLE,(/N_NODES,N_NODES,N_NODES/))
    call dset_create(chkpt,"/data/states/momz",H5T_NATIVE_DOUBLE,(/N_NODES,N_NODES,N_NODES/))
    call dset_create(chkpt,"/data/states/ener",H5T_NATIVE_DOUBLE,(/N_NODES,N_NODES,N_NODES/))

    call dset_create(chkpt,"/data/kernel/blend",H5T_NATIVE_DOUBLE,(/N_NODES,N_NODES,N_NODES,DGFV_QMAX/))

    call dset_create(chkpt,"/data/gravity/analytic",H5T_NATIVE_DOUBLE,(/N_NODES,N_NODES,N_NODES/))
    call dset_create(chkpt,"/data/gravity/analacc",H5T_NATIVE_DOUBLE,(/N_NODES,N_NODES,N_NODES,N_DIMS/))

    call dset_create(chkpt,"/data/gravity/potential",H5T_NATIVE_DOUBLE,(/N_NODES,N_NODES,N_NODES/))
    call dset_create(chkpt,"/data/gravity/acceleration",H5T_NATIVE_DOUBLE,(/N_NODES,N_NODES,N_NODES,N_DIMS/))

    call dset_create(chkpt,"/data/timestep/minloc",H5T_NATIVE_DOUBLE,(/2,2,2/))

end subroutine

subroutine hook_checkpoint_fill_100

    use share_mod, only: chkpt
    use checkpoint_utils_mod, only: dset_fill => checkpoint_utils_dataset_iterate

    call dset_fill(chkpt,"/mesh/levels",fill_levels_cb)
    call dset_fill(chkpt,"/mesh/morton",fill_morton_cb)

    call dset_fill(chkpt,"/mesh/coordinates",fill_coordinates_cb)
    call dset_fill(chkpt,"/mesh/block size",fill_block_size_cb)

    call dset_fill(chkpt,"/data/states/dens",fill_dens_cb)
    call dset_fill(chkpt,"/data/states/momx",fill_momx_cb)
    call dset_fill(chkpt,"/data/states/momy",fill_momy_cb)
    call dset_fill(chkpt,"/data/states/momz",fill_momz_cb)
    call dset_fill(chkpt,"/data/states/ener",fill_ener_cb)

    call dset_fill(chkpt,"/data/kernel/blend",fill_blend_cb)

    call dset_fill(chkpt,"/data/gravity/analytic",fill_analytic_cb)
    call dset_fill(chkpt,"/data/gravity/analacc",fill_analacc_cb)

    call dset_fill(chkpt,"/data/gravity/potential",fill_potential_cb)
    call dset_fill(chkpt,"/data/gravity/acceleration",fill_acceleration_cb)
    call dset_fill(chkpt,"/data/timestep/minloc",fill_dtminloc_cb)

end subroutine

subroutine fill_levels_cb(dset,quad)

    use mesh_types_mod, only: quad_t
    use checkpoint_types_mod, only: dset_t => checkpoint_dataset_t
    use checkpoint_utils_mod, only: dset_write => checkpoint_utils_dataset_write

    type(dset_t), intent(inout) :: dset
    type(quad_t), intent(inout) :: quad

    call dset_write(dset,quad,quad%level)

end subroutine

subroutine fill_morton_cb(dset,quad)

    use mesh_types_mod, only: quad_t
    use checkpoint_types_mod, only: dset_t => checkpoint_dataset_t
    use checkpoint_utils_mod, only: dset_write => checkpoint_utils_dataset_write

    type(dset_t), intent(inout) :: dset
    type(quad_t), intent(inout) :: quad

    call dset_write(dset,quad,quad%morton)

end subroutine

subroutine fill_coordinates_cb(dset,quad)

    use mesh_types_mod, only: quad_t
    use checkpoint_types_mod, only: dset_t => checkpoint_dataset_t
    use checkpoint_utils_mod, only: dset_write => checkpoint_utils_dataset_write

    type(dset_t), intent(inout) :: dset
    type(quad_t), intent(inout) :: quad

    call dset_write(dset,quad,quad%center)

end subroutine

subroutine fill_block_size_cb(dset,quad)

    use mesh_types_mod, only: quad_t
    use checkpoint_types_mod, only: dset_t => checkpoint_dataset_t
    use checkpoint_utils_mod, only: dset_write => checkpoint_utils_dataset_write

    type(dset_t), intent(inout) :: dset
    type(quad_t), intent(inout) :: quad

    call dset_write(dset,quad,quad%delta)

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

# if 1
subroutine fill_analytic_cb(dset,quad)

    use mesh_types_mod, only: quad_t
    use checkpoint_types_mod, only: dset_t => checkpoint_dataset_t
    use checkpoint_utils_mod, only: dset_write => checkpoint_utils_dataset_write

    type(dset_t), intent(inout) :: dset
    type(quad_t), intent(inout) :: quad
    logical                     :: ok

    call dset_write(dset,quad,quad%aux%gravity%rpot)

end subroutine

subroutine fill_analacc_cb(dset,quad)

    use mesh_types_mod, only: quad_t
    use checkpoint_types_mod, only: dset_t => checkpoint_dataset_t
    use checkpoint_utils_mod, only: dset_write => checkpoint_utils_dataset_write

    type(dset_t), intent(inout) :: dset
    type(quad_t), intent(inout) :: quad
    logical                     :: ok

    call dset_write(dset,quad,quad%aux%gravity%racc)

end subroutine

subroutine fill_potential_cb(dset,quad)

    use mesh_types_mod, only: quad_t
    use checkpoint_types_mod, only: dset_t => checkpoint_dataset_t
    use checkpoint_utils_mod, only: dset_write => checkpoint_utils_dataset_write

    type(dset_t), intent(inout) :: dset
    type(quad_t), intent(inout) :: quad
    logical                     :: ok

    call dset_write(dset,quad,quad%aux%gravity%pot)

end subroutine

subroutine fill_acceleration_cb(dset,quad)

    use mesh_types_mod, only: quad_t
    use checkpoint_types_mod, only: dset_t => checkpoint_dataset_t
    use checkpoint_utils_mod, only: dset_write => checkpoint_utils_dataset_write

    type(dset_t), intent(inout) :: dset
    type(quad_t), intent(inout) :: quad
    logical                     :: ok

    call dset_write(dset,quad,quad%aux%gravity%acc)

end subroutine
# endif

subroutine fill_dtminloc_cb(dset,quad)

    use share_mod, only: rt => runtime
    use mesh_types_mod, only: quad_t
    use checkpoint_types_mod, only: dset_t => checkpoint_dataset_t
    use checkpoint_utils_mod, only: dset_write => checkpoint_utils_dataset_write

    type(dset_t), intent(inout) :: dset
    type(quad_t), intent(inout) :: quad
    logical                     :: ok

    real(dp) :: dtminloc(2,2,2)

    !dtminloc = quad%aux%hydro%dtmin_count
    dtminloc = quad%level
    !dtminloc = rt%mpi%rank

    call dset_write(dset,quad,dtminloc)

    quad%aux%hydro%dtmin_count = 0

end subroutine

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

    blend = inflate(quad%aux%hydro%blend,1,DGFV_QMAX)
    !blend = inflate(quad%aux%hydro%cummulative_blend/quad%aux%hydro%blend_count,1,DGFV_QMAX)

    call dset_write(dset,quad,blend)

    quad%aux%hydro%cummulative_blend = 0.0_dp
    quad%aux%hydro%blend_count = 0.0_dp

end subroutine

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
