!# define PP_CHECKPOINT 1
!# define PP_CHECKPOINT_FLASH 1

module checkpoint_mod

use hdf5
use constants_mod

contains

subroutine hook_checkpoint_100

    use share_mod, only: chkpt
    use posix_mod, only: posix_mkdir
    use share_mod, only: rt => runtime, mesh
    use checkpoint_utils_mod, only: checkpoint_reset
    use mesh_utils_mod, only: mesh_alloc_tree

    character(len=10)   :: chkptstr

    integer :: istat

    call mesh_alloc_tree(mesh)

    call checkpoint_reset(chkpt)

    chkpt%iocount = rt%chkptcount

    write(chkptstr,'(I0.6)') rt%chkptcount
    chkpt%filename = trim(chkpt%directory) // '/' // trim(chkpt%prefix) // trim(chkptstr) // '.' // trim(chkpt%extension)

    if (rt%mpi%root) then
        istat = posix_mkdir(chkpt%directory)
    end if

end subroutine

subroutine hook_checkpoint_200

    use setup_config_mod
    use share_mod, only: rt => runtime,chkpt,mesh

# if PP_KERNEL_DGFV
    use kernel_const_mod, only: DGFV_QMIN
    use kernel_const_mod, only: DGFV_QMAX
# endif

    use checkpoint_utils_mod, only: checkpoint_add_real_runtime_parameter
    use checkpoint_utils_mod, only: checkpoint_add_string_runtime_parameter
    use checkpoint_utils_mod, only: checkpoint_add_logical_runtime_parameter
    use checkpoint_utils_mod, only: checkpoint_add_integer_runtime_parameter

    use checkpoint_utils_mod, only: checkpoint_add_real_scalar
    use checkpoint_utils_mod, only: checkpoint_add_string_scalar
    use checkpoint_utils_mod, only: checkpoint_add_logical_scalar
    use checkpoint_utils_mod, only: checkpoint_add_integer_scalar

    use hook_checkpoint_add_meta_mod, only: hook_checkpoint_add_meta

# if PP_N_DIMS == 2
    real(dp), parameter :: zlim(2) = (/0.0_dp,0.0_dp/)
# endif

    chkpt%setup_call = "none" // char(0)
    chkpt%file_creation_time = "none" // char(0)
    chkpt%flash_version = "FLASH 4.3_release" // char(0)
    chkpt%build_date = "none" // char(0)
    chkpt%build_dir = "none" // char(0)
    chkpt%build_machine = "none" // char(0)
    chkpt%cflags = "none" // char(0)
    chkpt%fflags = "none" // char(0)
    chkpt%setup_time_stamp = "none" // char(0)
    chkpt%build_time_stamp = "none" // char(0)

    call checkpoint_add_integer_runtime_parameter(chkpt,"lrefine_min",mesh%minlevel+1)
    call checkpoint_add_integer_runtime_parameter(chkpt,"lrefine_max",mesh%maxlevel+1)

    call checkpoint_add_integer_runtime_parameter(chkpt,"checkpointfilenumber",chkpt%iocount)
    call checkpoint_add_integer_runtime_parameter(chkpt,"checkpointfileintervalstep",0)
    call checkpoint_add_integer_runtime_parameter(chkpt,"fileformatversion",chkpt%file_format_version)
    call checkpoint_add_integer_runtime_parameter(chkpt,"forcedplotfilenumber",0)

    call checkpoint_add_real_runtime_parameter(chkpt,"xmin",xlim(1))
    call checkpoint_add_real_runtime_parameter(chkpt,"xmax",xlim(2))

    call checkpoint_add_real_runtime_parameter(chkpt,"ymin",ylim(1))
    call checkpoint_add_real_runtime_parameter(chkpt,"ymax",ylim(2))

    call checkpoint_add_real_runtime_parameter(chkpt,"zmin",zlim(1))
    call checkpoint_add_real_runtime_parameter(chkpt,"zmax",zlim(2))

    call checkpoint_add_real_runtime_parameter(chkpt,"cfl",cfl)
    call checkpoint_add_real_runtime_parameter(chkpt,"checkpointfileintervaltime",dt_checkpoint)

    call checkpoint_add_real_runtime_parameter(chkpt,"tinitial",init_time)
    call checkpoint_add_real_runtime_parameter(chkpt,"tmax",stop_time)

    call checkpoint_add_real_runtime_parameter(chkpt,"gamma",kappa)

# if MODULE_SELF_GRAVITY
    call checkpoint_add_real_runtime_parameter(chkpt,"grav_constant",grav)
# endif

    call checkpoint_add_logical_runtime_parameter(chkpt,"restart",.false.)

    call checkpoint_add_string_runtime_parameter(chkpt,"unitsystem","none")
    call checkpoint_add_string_runtime_parameter(chkpt,"geometry","cartesian")
    call checkpoint_add_string_runtime_parameter(chkpt,"basen",chkpt%prefix)

# if PP_MESH_PERIODIC
    call checkpoint_add_string_runtime_parameter(chkpt,"xl_boundary_type","periodic")
    call checkpoint_add_string_runtime_parameter(chkpt,"xr_boundary_type","periodic")

    call checkpoint_add_string_runtime_parameter(chkpt,"yl_boundary_type","periodic")
    call checkpoint_add_string_runtime_parameter(chkpt,"yr_boundary_type","periodic")

# if PP_N_DIMS > 2
    call checkpoint_add_string_runtime_parameter(chkpt,"zl_boundary_type","periodic")
    call checkpoint_add_string_runtime_parameter(chkpt,"zr_boundary_type","periodic")
# endif
# else
    call checkpoint_add_string_runtime_parameter(chkpt,"xl_boundary_type","outflow")
    call checkpoint_add_string_runtime_parameter(chkpt,"xr_boundary_type","outflow")

    call checkpoint_add_string_runtime_parameter(chkpt,"yl_boundary_type","outflow")
    call checkpoint_add_string_runtime_parameter(chkpt,"yr_boundary_type","outflow")

# if PP_N_DIMS > 2
    call checkpoint_add_string_runtime_parameter(chkpt,"zl_boundary_type","outflow")
    call checkpoint_add_string_runtime_parameter(chkpt,"zr_boundary_type","outflow")
# endif
# endif

    call checkpoint_add_integer_scalar(chkpt,"nstep",rt%timesteps)

    call checkpoint_add_integer_scalar(chkpt,"nxb",N_NODES)
    call checkpoint_add_integer_scalar(chkpt,"nyb",N_NODES)
# if PP_N_DIMS > 2
    call checkpoint_add_integer_scalar(chkpt,"nzb",N_NODES)
# else
    call checkpoint_add_integer_scalar(chkpt,"nzb",1)
# endif

    call checkpoint_add_integer_scalar(chkpt,"dimensionality",N_DIMS)
    call checkpoint_add_integer_scalar(chkpt,"globalnumblocks",INT(mesh%tree%numblocks,kind=4))
    call checkpoint_add_integer_scalar(chkpt,"checkpointfilenumber",chkpt%iocount)

    call checkpoint_add_real_scalar(chkpt,"time",rt%simtime)
    call checkpoint_add_real_scalar(chkpt,"dt",rt%timestep)
    call checkpoint_add_real_scalar(chkpt,"dtold",rt%timestep)
    call checkpoint_add_real_scalar(chkpt,"dtnew",rt%timestep)

    call checkpoint_add_real_scalar(chkpt,"nextcheckpointtime",rt%nextdump)

    call checkpoint_add_logical_scalar(chkpt,"doubleprecision",.true.)

    call checkpoint_add_string_scalar(chkpt,"geometry","cartesian")
    call checkpoint_add_string_scalar(chkpt,"nemo","test")

# if PP_KERNEL_DGFV
    call checkpoint_add_integer_scalar(chkpt,"DGFV_QMIN",DGFV_QMIN)
    call checkpoint_add_integer_scalar(chkpt,"DGFV_QMAX",DGFV_QMAX)
# endif

    call hook_checkpoint_add_meta

end subroutine

subroutine hook_checkpoint_300

    use share_mod, only: mesh, chkpt, rt => runtime
    use h5simple_mod, only: h5simple_make_dset
    use h5simple_mod, only: h5simple_create_path
    use hook_checkpoint_create_mod, only: hook_checkpoint_create
    use checkpoint_utils_mod, only: checkpoint_write_sim_info

    use checkpoint_utils_mod, only: checkpoint_write_real_runtime_parameters
    use checkpoint_utils_mod, only: checkpoint_write_string_runtime_parameters
    use checkpoint_utils_mod, only: checkpoint_write_logical_runtime_parameters
    use checkpoint_utils_mod, only: checkpoint_write_integer_runtime_parameters

    use checkpoint_utils_mod, only: checkpoint_write_real_scalars
    use checkpoint_utils_mod, only: checkpoint_write_string_scalars
    use checkpoint_utils_mod, only: checkpoint_write_logical_scalars
    use checkpoint_utils_mod, only: checkpoint_write_integer_scalars

    use checkpoint_utils_mod, only: checkpoint_write_unknown_names

    use checkpoint_utils_mod, only: dset_create => checkpoint_utils_dataset_create

    use mesh_const_mod, only: MESH_FACES, MESH_CHILDREN
    use setup_config_mod

    use morton_mod, only: morton_to_vertex
    use morton_mod, only: morton_get_child_id

    use h5simple_mod, only: h5simple_make_dset_double
    use h5simple_mod, only: h5simple_make_dset_integr

    integer(HID_T)      :: fid,pid
    integer             :: ierr

    integer(kind=8) :: iquad

    real(dp), allocatable :: block_size(:,:)
    real(dp), allocatable :: bounding_box(:,:,:)
    real(dp), allocatable :: coordinates(:,:)

    integer, allocatable :: which_child(:)
    integer, allocatable :: bflags(:)
    integer, allocatable :: nodetype(:)

    logical :: ok

# if PP_N_DIMS == 2
    real(dp), parameter :: zlim(2) = (/0.0_dp,0.0_dp/)
# endif
    real(dp), parameter :: domsize(3) = (/xlim(2)-xlim(1),ylim(2)-ylim(1),zlim(2)-zlim(1)/)

    if (.not.rt%mpi%root) return

    call h5open_f(ierr)

    !! ------------------------------------------------------------------ !!
    !! Open hdf5 file and define some common dataset properties.
    !! In case, truncate already existing file.

    call h5fcreate_f(trim(chkpt%filename) // trim(chkpt%partial_suffix), H5F_ACC_TRUNC_F, fid, ierr)

    !! Sets the timing for storage space allocation. Here, immediately
    !! allocate the needed memory space, which is important for parallel write.
    call h5pcreate_f(H5P_DATASET_CREATE_F, pid, ierr)
    call h5pset_alloc_time_f(pid, H5D_ALLOC_TIME_EARLY_F, ierr)

    !! found in compound dtype example
    CALL h5pcreate_f(H5P_DATASET_XFER_F, pid, ierr)
    CALL h5pset_preserve_f(pid, .true., ierr)

    chkpt%fileid = fid
    chkpt%propertyid = pid

    !! ------------------------------------------------------------------ !!

    call checkpoint_write_sim_info(chkpt)

    call checkpoint_write_real_runtime_parameters(chkpt)
    call checkpoint_write_string_runtime_parameters(chkpt)
    call checkpoint_write_logical_runtime_parameters(chkpt)
    call checkpoint_write_integer_runtime_parameters(chkpt)

    call checkpoint_write_real_scalars(chkpt)
    call checkpoint_write_string_scalars(chkpt)
    call checkpoint_write_logical_scalars(chkpt)
    call checkpoint_write_integer_scalars(chkpt)

    call checkpoint_write_unknown_names(chkpt)

    !! ------------------------------------------------------------------ !!

    allocate(block_size(3,mesh%tree%numblocks))
    allocate(coordinates(3,mesh%tree%numblocks))
    allocate(bounding_box(2,3,mesh%tree%numblocks))

    allocate(which_child(mesh%tree%numblocks))
    allocate(bflags(mesh%tree%numblocks))
    allocate(nodetype(mesh%tree%numblocks))

    do iquad = 1,mesh%tree%numblocks
        block_size(:,iquad) =  domsize / REAL(2**mesh%tree%levels(iquad),dp)
        coordinates(:,iquad) = (/xlim(1),ylim(1),zlim(1)/) &
# if PP_N_DIMS == 2
            + domsize * (/morton_to_vertex(mesh%tree%morton(:,iquad)),0.0_dp/) &
# else
            + domsize * morton_to_vertex(mesh%tree%morton(:,iquad)) &
# endif
            + 0.5_dp * block_size(:,iquad)

        bounding_box(1,:,iquad) = coordinates(:,iquad) - 0.5_dp*block_size(:,iquad)
        bounding_box(2,:,iquad) = coordinates(:,iquad) + 0.5_dp*block_size(:,iquad)

        which_child(iquad) = morton_get_child_id(INT(mesh%tree%levels(iquad),kind=4),mesh%tree%morton(:,iquad))
    end do

    bflags = -1
    which_child(1) = -1
    nodetype = 2
    nodetype(mesh%tree%noffset(mesh%tree%maxlevel)+1:) = 1

    ok = h5simple_make_dset_double(chkpt%fileid,"/coordinates",coordinates,shape(coordinates))
    ok = h5simple_make_dset_double(chkpt%fileid,"/block size",block_size,shape(block_size))
    ok = h5simple_make_dset_double(chkpt%fileid,"/bounding box",bounding_box,shape(bounding_box))

    ok = h5simple_make_dset_integr(chkpt%fileid,"/gid",mesh%tree%gid,shape(mesh%tree%gid))
    ok = h5simple_make_dset_integr(chkpt%fileid,"/node type",nodetype,shape(nodetype))
    ok = h5simple_make_dset_integr(chkpt%fileid,"/refine level",INT(mesh%tree%levels+1,kind=4),shape(mesh%tree%levels))
    ok = h5simple_make_dset_integr(chkpt%fileid,"/which child",which_child,shape(which_child))
    ok = h5simple_make_dset_integr(chkpt%fileid,"/bflags",bflags,shape(bflags))

    deallocate(block_size)
    deallocate(coordinates)
    deallocate(bounding_box)

    deallocate(bflags)
    deallocate(nodetype)
    deallocate(which_child)

    !! ------------------------------------------------------------------ !!

    call hook_checkpoint_create

    !! ------------------------------------------------------------------ !!

    call h5pclose_f(pid, ierr)
    call h5fclose_f(fid, ierr)
    call h5close_f(ierr)

end subroutine

subroutine hook_checkpoint_500

    use share_mod, only: chkpt
    use parallel_mod, only: barrier
    use share_mod, only: rt => runtime
    use posix_mod, only: posix_rename
    use hook_checkpoint_fill_mod, only: hook_checkpoint_fill

    use checkpoint_utils_mod, only: dset_fill => checkpoint_utils_dataset_iterate

    integer(HID_T) :: fid, pid, dpid
    integer :: istat
    integer :: ierr

    call h5open_f(ierr)

    !! Parallel access of the file in read-write mode.
    CALL h5pcreate_f(H5P_FILE_ACCESS_F, pid, ierr)
    CALL h5pset_fapl_mpio_f(pid, rt%mpi%comm%mpi_val, rt%mpi%info%mpi_val, ierr)

    call barrier(istat)

    CALL h5fopen_f(trim(chkpt%filename) // trim(chkpt%partial_suffix), H5F_ACC_RDWR_F, fid, ierr, pid)

    !! Set parallel write mode to "independent".
    CALL h5pcreate_f(H5P_DATASET_XFER_F, pid, ierr)
    CALL h5pset_dxpl_mpio_f(pid, H5FD_MPIO_INDEPENDENT_F, ierr)

    chkpt%fileid = fid
    chkpt%propertyid = pid

    call hook_checkpoint_fill

    call h5pclose_f(pid, ierr)
    call h5fclose_f(fid, ierr)
    call h5close_f(ierr)

    call barrier(istat)

    if (rt%mpi%root) then
        istat = posix_rename(trim(chkpt%filename) // trim(chkpt%partial_suffix), trim(chkpt%filename))
    end if

    rt%chkptcount = rt%chkptcount + 1

end subroutine

!! subroutine fill_levels_cb(dset,quad)
!! 
!!     use mesh_types_mod, only: quad_t
!!     use checkpoint_types_mod, only: dset_t => checkpoint_dataset_t
!!     use checkpoint_utils_mod, only: dset_write => checkpoint_utils_dataset_write
!! 
!!     type(dset_t), intent(inout) :: dset
!!     type(quad_t), intent(inout) :: quad
!! 
!!     call dset_write(dset,quad,quad%level)
!! 
!! end subroutine
!! 
!! subroutine fill_node_type_cb(dset,quad)
!! 
!!     use mesh_types_mod, only: quad_t
!!     use checkpoint_types_mod, only: dset_t => checkpoint_dataset_t
!!     use checkpoint_utils_mod, only: dset_write => checkpoint_utils_dataset_write
!! 
!!     type(dset_t), intent(inout) :: dset
!!     type(quad_t), intent(inout) :: quad
!! 
!!     call dset_write(dset,quad,INT(1,kind=1))
!! 
!! end subroutine
!! 
!! subroutine fill_which_child_cb(dset,quad)
!! 
!!     use mesh_types_mod, only: quad_t
!!     use checkpoint_types_mod, only: dset_t => checkpoint_dataset_t
!!     use checkpoint_utils_mod, only: dset_write => checkpoint_utils_dataset_write
!!     use p8wrap_interfaces_mod, only: p8wrap_quadrant_child_id
!! 
!!     type(dset_t), intent(inout) :: dset
!!     type(quad_t), intent(inout) :: quad
!! 
!!     integer(kind=1) :: child
!! 
!!     child = INT(1 + p8wrap_quadrant_child_id(quad%level,quad%morton(1),quad%morton(2),quad%morton(3)),kind=1)
!!     call dset_write(dset,quad,child)
!! 
!! end subroutine
!! 
!! subroutine fill_gid_cb(dset,quad)
!! 
!!     use mesh_types_mod, only: quad_t
!!     use checkpoint_types_mod, only: dset_t => checkpoint_dataset_t
!!     use checkpoint_utils_mod, only: dset_write => checkpoint_utils_dataset_write
!! 
!!     use mesh_const_mod, only: MESH_FACES, MESH_CHILDREN
!! 
!!     type(dset_t), intent(inout) :: dset
!!     type(quad_t), intent(inout) :: quad
!! 
!!     integer(kind=8) :: gid(MESH_FACES+1+MESH_CHILDREN)
!! 
!!     gid = quad%gloid
!! 
!!     call dset_write(dset,quad,gid)
!! 
!! end subroutine

!! subroutine fill_morton_cb(dset,quad)
!! 
!!     use mesh_types_mod, only: quad_t
!!     use checkpoint_types_mod, only: dset_t => checkpoint_dataset_t
!!     use checkpoint_utils_mod, only: dset_write => checkpoint_utils_dataset_write
!! 
!!     type(dset_t), intent(inout) :: dset
!!     type(quad_t), intent(inout) :: quad
!! 
!!     call dset_write(dset,quad,quad%morton)
!! 
!! end subroutine

!! subroutine fill_coordinates_cb(dset,quad)
!! 
!!     use mesh_types_mod, only: quad_t
!!     use checkpoint_types_mod, only: dset_t => checkpoint_dataset_t
!!     use checkpoint_utils_mod, only: dset_write => checkpoint_utils_dataset_write
!! 
!!     type(dset_t), intent(inout) :: dset
!!     type(quad_t), intent(inout) :: quad
!! 
!!     call dset_write(dset,quad,quad%center)
!! 
!! end subroutine
!! 
!! subroutine fill_block_size_cb(dset,quad)
!! 
!!     use mesh_types_mod, only: quad_t
!!     use checkpoint_types_mod, only: dset_t => checkpoint_dataset_t
!!     use checkpoint_utils_mod, only: dset_write => checkpoint_utils_dataset_write
!! 
!!     type(dset_t), intent(inout) :: dset
!!     type(quad_t), intent(inout) :: quad
!! 
!!     call dset_write(dset,quad,quad%delta)
!! 
!! end subroutine
!! 
!! subroutine fill_bounding_box_cb(dset,quad)
!! 
!!     use mesh_types_mod, only: quad_t
!!     use checkpoint_types_mod, only: dset_t => checkpoint_dataset_t
!!     use checkpoint_utils_mod, only: dset_write => checkpoint_utils_dataset_write
!! 
!!     type(dset_t), intent(inout) :: dset
!!     type(quad_t), intent(inout) :: quad
!! 
!!     real(dp) :: bounding_box(2,3)
!! 
!!     bouding_box = 0.0_dp
!!     bounding_box(1,1:N_DIMS) = quad%center - 0.5_dp*quad%delta
!!     bounding_box(2,1:N_DIMS) = quad%center + 0.5_dp*quad%delta
!! 
!!     call dset_write(dset,quad,bounding_box)
!! 
!! end subroutine

end module
