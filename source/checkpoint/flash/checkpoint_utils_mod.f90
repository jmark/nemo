#:def expand(n,fmtstr)
${','.join([fmtstr.format(i) for i in range(1,n+1)])}$
#:enddef expand

#:set typenames = 'int8 int32 int64 real64'.split()
#:set typeids   = 'integer(c_int8_t) integer(c_int32_t) integer(c_int64_t) real(dp)'.split()
#:set ranks     = (0,5+1)

module checkpoint_utils_mod

use constants_mod

interface checkpoint_utils_dataset_write
#:for typename,typeid in zip(typenames,typeids) 
#:for rank in range(ranks[0],ranks[1])
    procedure checkpoint_utils_dataset_write_${typename}$_${rank}$d
#:endfor
#:endfor
end interface

contains

subroutine checkpoint_reset(chkpt)

    use checkpoint_types_mod, only: checkpoint_t

    type(checkpoint_t), intent(inout) :: chkpt

    chkpt%index_real_runtime_parameters = 0
    chkpt%index_string_runtime_parameters = 0
    chkpt%index_integer_runtime_parameters = 0
    chkpt%index_logical_runtime_parameters = 0

    chkpt%index_real_scalars = 0
    chkpt%index_string_scalars = 0
    chkpt%index_integer_scalars = 0
    chkpt%index_logical_scalars = 0

    chkpt%index_unknown_names = 0

end subroutine

subroutine checkpoint_utils_dataset_create(chkpt,dpath,dtype,dshape)

    use hdf5
    use share_mod, only: mesh
    use h5simple_mod, only: h5simple_create_path
    use checkpoint_types_mod, only: chkpt_t => checkpoint_t

    type(chkpt_t), intent(inout)    :: chkpt
    character(len=*), intent(in)    :: dpath
    integer(HID_T), intent(in)      :: dtype
    integer, intent(in), optional   :: dshape(:)
    logical                         :: ok

    integer(HID_T)  :: pid
    integer(kind=8) :: buffer_shape(16)
    integer(kind=4) :: buffer_rank

    integer(HID_T) :: dsetid, spaceid

    integer :: ierr

    !! SZIP compression (not working so far)
    !INTEGER :: szip_options_mask
    !INTEGER :: szip_pixels_per_block

    if (present(dshape)) then
        buffer_rank = size(dshape) + 1
        buffer_shape = 0
        buffer_shape(1:size(dshape)) = dshape
        ! buffer_shape(size(dshape)+1) = mesh%nglobal
        buffer_shape(size(dshape)+1) = mesh%tree%numblocks
    else
        buffer_rank = 1
        buffer_shape = 0
        ! buffer_shape(1) = mesh%nglobal
        buffer_shape(1) = mesh%tree%numblocks
    end if

    !! Sets the timing for storage space allocation. Here, immediately
    !! allocate the needed memory space, which is important for parallel write.
    call h5pcreate_f(H5P_DATASET_CREATE_F, pid, ierr)
    call h5pset_alloc_time_f(pid, H5D_ALLOC_TIME_EARLY_F, ierr)

    ok = h5simple_create_path(chkpt%fileid,trim(dpath))

    call h5screate_simple_f(buffer_rank,INT(buffer_shape,kind=HSIZE_T),spaceid,ierr)
    call h5dcreate_f(chkpt%fileid,trim(dpath),dtype,spaceid,dsetid,ierr,dcpl_id=pid)

    call h5sclose_f(spaceid,ierr)
    call h5dclose_f(dsetid,ierr)

    call h5pclose_f(pid, ierr)

end subroutine

subroutine checkpoint_utils_dataset_set_attribute(chkpt,dpath,dtype,attrkey,attrval)

    use hdf5
    use share_mod, only: mesh
    use h5simple_mod, only: h5simple_create_path
    use h5simple_mod, only: h5simple_create_dataset
    use checkpoint_types_mod, only: chkpt_t => checkpoint_t
    use iso_c_binding, only: c_loc

    type(chkpt_t), intent(inout)    :: chkpt
    character(len=*), intent(in)    :: dpath
    integer(HID_T), intent(in)      :: dtype
    character(len=*), intent(in)    :: attrkey
    real(dp), intent(in), target    :: attrval

    integer(HID_T)  :: pid
    integer(kind=8) :: buffer_shape(16)
    integer(kind=4) :: buffer_rank

    integer(HID_T) :: dsetid, spaceid, attrid

    integer :: ierr

    call h5dopen_f(chkpt%fileid, trim(dpath), dsetid, ierr)
    CALL h5screate_simple_f(1, (/INT(1,kind=HSIZE_T)/), spaceid, ierr)
    CALL H5Acreate_f(dsetid,attrkey,dtype,spaceid,attrid, ierr)
    CALL H5awrite_f(attrid,H5T_NATIVE_DOUBLE,c_loc(attrval),ierr)

    call h5aclose_f(attrid,ierr)
    call h5sclose_f(spaceid,ierr)
    call h5dclose_f(dsetid,ierr)

end subroutine

!! ========================================================================== !!

subroutine checkpoint_utils_dataset_iterate(chkpt,dpath,callback)

    use mpi_f08
    use hdf5
    use share_mod, only: mesh
    use share_mod, only: rt => runtime
    use checkpoint_types_mod

    type(checkpoint_t), intent(inout), target   :: chkpt
    character(len=*), intent(in)                :: dpath
    procedure(checkpoint_iterate_cb_t)          :: callback

    type(checkpoint_dataset_t) :: dset

    integer          :: drank,ierr
    integer(HSIZE_T) :: offset(16), tsize
    integer(HSIZE_T) :: dims(16),maxdims(16)
    integer(HSIZE_T) :: dshape(16),bshape(16)
    integer(HID_T)   :: dsetid, dspaceid, dtypeid

    integer(kind=8)  :: i

    dset%chkpt => chkpt
    dset%path = trim(dpath)

    call h5dopen_f(chkpt%fileid, trim(dpath), dsetid, ierr)
    call H5Dget_type_f(dsetid, dtypeid, ierr)

    call h5tget_size_f(dtypeid,tsize,ierr) 
    call h5dget_space_f(dsetid, dspaceid, ierr)
    call h5sget_simple_extent_ndims_f(dspaceid,drank,ierr)
    call h5sget_simple_extent_dims_f(dspaceid,dims,maxdims,ierr)

    dshape          = 1
    dshape(1:drank) = dims(1:drank)

    bshape           = 1
    bshape(1:drank)  = dshape
    bshape(  drank)  = mesh%nquads

    !! datatype size
    dset%tsize = tsize
    !! buffer size
    dset%bsize = tsize * product(bshape(1:drank))
    !! total dataset size
    dset%dsize = tsize * merge(INT(1,kind=8),product(bshape(1:drank-1)),drank < 2)

    allocate(dset%buffer(dset%bsize))

    !!$OMP PARALLEL DO private(i) shared(dset,chkpt)
    do i = 1,mesh%nquads
        call callback(dset,mesh%quads(i))
    end do
    !!$OMP END PARALLEL DO

    offset        = 0
    ! offset(drank) = mesh%iglobal
    offset(drank) = mesh%iglobal + mesh%tree%numblocks - mesh%nglobal

    ! write (*,'(a10,99(i7))') 'offset', rt%mpi%rank, mesh%tree%numblocks, mesh%nglobal, mesh%iglobal, offset(1:9)
    call checkpoint_utils_write_dataset(chkpt%fileid, dtypeid, drank, &
        offset(1:drank), bshape(1:drank), dset%buffer, trim(dset%path), chkpt%propertyid)

    deallocate(dset%buffer)

    call H5Sclose_f(dspaceid, ierr)
    call H5Tclose_f(dtypeid, ierr)
    call H5Dclose_f(dsetid, ierr)

end subroutine

!! ========================================================================== !!

#:for typename,typeid in zip(typenames,typeids) 
#:for rank in range(ranks[0],ranks[1])
subroutine checkpoint_utils_dataset_write_${typename}$_${rank}$d(dset,quad,cdata)

    use iso_c_binding
    use mesh_types_mod, only: quad_t
    use checkpoint_types_mod, only: dset_t => checkpoint_dataset_t

    type(dset_t), intent(inout)     :: dset
    type(quad_t), intent(in)        :: quad
#:if rank < 1
    ${typeid}$, intent(in)          :: cdata
#:else
    ${typeid}$, intent(in)          :: cdata(${expand(rank,":")}$)
#:endif

    integer(kind=8)                 :: a
    character(len=1), parameter     :: dummy = 'a'

    a = dset%dsize*(quad%locid-1)
    dset%buffer(a + 1:a + dset%dsize) = transfer(cdata,mold=dummy,size=dset%dsize)

end subroutine
#:endfor
#:endfor

!! ========================================================================== !!

subroutine checkpoint_utils_write_dataset(fileid, typeid, drank, offset, dshape, dbuffer, path, plist)

    use hdf5

    integer(HID_T), intent(in)      :: fileid
    integer(HID_T), intent(in)      :: typeid
    integer, intent(in)             :: drank
    integer(HSIZE_T), intent(in)    :: offset(:)
    integer(HSIZE_T), intent(in)    :: dshape(:)
    character(len=1), intent(in)    :: dbuffer(*)

    character(len=*), intent(in)            :: path
    integer(HID_T), intent(in), optional    :: plist

    integer             :: ierr

    integer(HID_T)      :: dsetid
    integer(HID_T)      :: dspaceid
    integer(HID_T)      :: memspaceid

    call h5dopen_f(fileid, trim(path), dsetid, ierr)
    call h5dget_space_f(dsetid, dspaceid, ierr)

    call h5sselect_hyperslab_f(dspaceid, H5S_SELECT_SET_F, offset, dshape, ierr) 
    call h5screate_simple_f(drank, dshape, memspaceid, ierr)

    if (present(plist)) then
        call h5dwrite_f(dsetid, typeid, dbuffer, dshape, ierr, &
            mem_space_id=memspaceid, file_space_id = dspaceid, xfer_prp = plist)
    else
        call h5dwrite_f(dsetid, typeid, dbuffer, dshape, ierr, &
            mem_space_id=memspaceid, file_space_id = dspaceid)
    end if

    call h5sclose_f(memspaceid, ierr)
    call h5sclose_f(dspaceid, ierr)
    call h5dclose_f(dsetid, ierr)

end subroutine

subroutine checkpoint_add_integer_runtime_parameter(chkpt,key,val)

    use checkpoint_types_mod, only: checkpoint_t

    type(checkpoint_t), intent(inout)   :: chkpt
    character(len=*), intent(in)        :: key
    integer, intent(in)                 :: val

    chkpt%index_integer_runtime_parameters = chkpt%index_integer_runtime_parameters + 1
    chkpt%keys_integer_runtime_parameters(chkpt%index_integer_runtime_parameters) = trim(key) // char(0)
    chkpt%vals_integer_runtime_parameters(chkpt%index_integer_runtime_parameters) = val

end subroutine

subroutine checkpoint_add_real_runtime_parameter(chkpt,key,val)

    use checkpoint_types_mod, only: checkpoint_t

    type(checkpoint_t), intent(inout)   :: chkpt
    character(len=*), intent(in)        :: key
    real(dp), intent(in)                :: val

    chkpt%index_real_runtime_parameters = chkpt%index_real_runtime_parameters + 1
    chkpt%keys_real_runtime_parameters(chkpt%index_real_runtime_parameters) = trim(key) // char(0)
    chkpt%vals_real_runtime_parameters(chkpt%index_real_runtime_parameters) = val

end subroutine

subroutine checkpoint_add_string_runtime_parameter(chkpt,key,val)

    use checkpoint_types_mod, only: checkpoint_t

    type(checkpoint_t), intent(inout)   :: chkpt
    character(len=*), intent(in)        :: key
    character(len=*), intent(in)        :: val

    chkpt%index_string_runtime_parameters = chkpt%index_string_runtime_parameters + 1
    chkpt%keys_string_runtime_parameters(chkpt%index_string_runtime_parameters) = trim(key) // char(0)
    chkpt%vals_string_runtime_parameters(chkpt%index_string_runtime_parameters) = trim(val) // char(0)

end subroutine

subroutine checkpoint_add_logical_runtime_parameter(chkpt,key,val)

    use checkpoint_types_mod, only: checkpoint_t

    type(checkpoint_t), intent(inout)   :: chkpt
    character(len=*), intent(in)        :: key
    logical, intent(in)                 :: val

    chkpt%index_logical_runtime_parameters = chkpt%index_logical_runtime_parameters + 1
    chkpt%keys_logical_runtime_parameters(chkpt%index_logical_runtime_parameters) = trim(key) // char(0)
    chkpt%vals_logical_runtime_parameters(chkpt%index_logical_runtime_parameters) = val

end subroutine

subroutine checkpoint_add_integer_scalar(chkpt,key,val)

    use checkpoint_types_mod, only: checkpoint_t

    type(checkpoint_t), intent(inout)   :: chkpt
    character(len=*), intent(in)        :: key
    integer, intent(in)                 :: val

    chkpt%index_integer_scalars = chkpt%index_integer_scalars + 1
    chkpt%keys_integer_scalars(chkpt%index_integer_scalars) = trim(key) // char(0)
    chkpt%vals_integer_scalars(chkpt%index_integer_scalars) = val

end subroutine

subroutine checkpoint_add_real_scalar(chkpt,key,val)

    use checkpoint_types_mod, only: checkpoint_t

    type(checkpoint_t), intent(inout)   :: chkpt
    character(len=*), intent(in)        :: key
    real(dp), intent(in)                :: val

    chkpt%index_real_scalars = chkpt%index_real_scalars + 1
    chkpt%keys_real_scalars(chkpt%index_real_scalars) = trim(key) // char(0)
    chkpt%vals_real_scalars(chkpt%index_real_scalars) = val

end subroutine

subroutine checkpoint_add_string_scalar(chkpt,key,val)

    use checkpoint_types_mod, only: checkpoint_t

    type(checkpoint_t), intent(inout)   :: chkpt
    character(len=*), intent(in)        :: key
    character(len=*), intent(in)        :: val

    chkpt%index_string_scalars = chkpt%index_string_scalars + 1
    chkpt%keys_string_scalars(chkpt%index_string_scalars) = trim(key) // char(0)
    chkpt%vals_string_scalars(chkpt%index_string_scalars) = trim(val) // char(0)

end subroutine

subroutine checkpoint_add_logical_scalar(chkpt,key,val)

    use checkpoint_types_mod, only: checkpoint_t

    type(checkpoint_t), intent(inout)   :: chkpt
    character(len=*), intent(in)        :: key
    logical, intent(in)                 :: val

    chkpt%index_logical_scalars = chkpt%index_logical_scalars + 1
    chkpt%keys_logical_scalars(chkpt%index_logical_scalars) = trim(key) // char(0)
    chkpt%vals_logical_scalars(chkpt%index_logical_scalars) = val

end subroutine

subroutine checkpoint_add_unknown_name(chkpt,unkname)

    use checkpoint_types_mod, only: checkpoint_t

    type(checkpoint_t), intent(inout)   :: chkpt
    character(len=*), intent(in)        :: unkname

    chkpt%index_unknown_names = chkpt%index_unknown_names + 1
    chkpt%unknown_names(chkpt%index_unknown_names) = trim(unkname)

end subroutine

subroutine checkpoint_write_sim_info(chkpt)

    use hdf5
    use checkpoint_types_mod, only: checkpoint_t

    type(checkpoint_t), intent(inout) :: chkpt

    INTEGER(HID_T)  :: dset_id   ! Dset identifier
    INTEGER(HID_T)  :: dspace_id ! Dataspace identifier

    INTEGER(HID_T)  :: siminfo_dtype
    INTEGER(HID_T)  :: string_dtype80
    INTEGER(HID_T)  :: string_dtype400

    INTEGER(HID_T)  :: memory_dtypes(11)   ! Memory datatype identifier 

    INTEGER(SIZE_T) :: type_sizeC ! Size of compound datatype

    INTEGER(SIZE_T) :: type_sizeS80 ! Size of name datatype 
    INTEGER(SIZE_T) :: type_sizeS400 ! Size of name datatype 

    INTEGER(SIZE_T) :: type_sizeI ! Size of value datatype
    INTEGER(SIZE_T) :: offset     ! Member's offset

    INTEGER         :: error      ! Error flag

    integer, parameter :: strlens(10) = (/ 400,80,80,80,80,80,400,400,80,80 /)

    integer :: i

    !! Calculate total size by calculating sizes of each member
    call h5tcopy_f(H5T_C_S1,string_dtype80,error)
    call h5tcopy_f(H5T_C_S1,string_dtype400,error)

    call h5tset_size_f(string_dtype80,INT(80,kind=SIZE_T),error)
    call h5tset_size_f(string_dtype400,INT(400,kind=SIZE_T),error)

    call h5tget_size_f(string_dtype80,type_sizeS80,error)
    call h5tget_size_f(string_dtype400,type_sizeS400,error)
    call h5tget_size_f(H5T_NATIVE_INTEGER,type_sizeI,error)

    type_sizeC = type_sizeI + 3*type_sizeS400 + 7*type_sizeS80

    CALL h5tcreate_f(H5T_COMPOUND_F,type_sizeC,siminfo_dtype,error)

    offset = 0
    call h5tinsert_f(siminfo_dtype,"file format version",offset,H5T_NATIVE_INTEGER,error)
     
    offset = offset + type_sizeI
    call h5tinsert_f(siminfo_dtype,"setup call",offset,string_dtype400,error)

    offset = offset + type_sizeS400
    call h5tinsert_f(siminfo_dtype,"file creation time",offset,string_dtype80,error)

    offset = offset + type_sizeS80
    call h5tinsert_f(siminfo_dtype,"flash version",offset,string_dtype80,error)

    offset = offset + type_sizeS80
    call h5tinsert_f(siminfo_dtype,"build date",offset,string_dtype80,error)

    offset = offset + type_sizeS80
    call h5tinsert_f(siminfo_dtype,"build dir",offset,string_dtype80,error)

    offset = offset + type_sizeS80
    call h5tinsert_f(siminfo_dtype,"build machine",offset,string_dtype80,error)

    offset = offset + type_sizeS80
    call h5tinsert_f(siminfo_dtype,"cflags",offset,string_dtype400,error)

    offset = offset + type_sizeS400
    call h5tinsert_f(siminfo_dtype,"fflags",offset,string_dtype400,error)

    offset = offset + type_sizeS400
    call h5tinsert_f(siminfo_dtype,"setup time stamp",offset,string_dtype80,error)

    offset = offset + type_sizeS80
    call h5tinsert_f(siminfo_dtype,"build time stamp",offset,string_dtype80,error)

    !! data space
    call h5screate_simple_f(1,(/INT(1,kind=SIZE_T)/),dspace_id,error)

    !! dataset
    call h5dcreate_f(chkpt%fileid,"sim info",siminfo_dtype,dspace_id,dset_id,error)

    !! Create memory types.
# if 1
    offset = 0
    memory_dtypes = 0

    i = 1
    call h5tcreate_f(H5T_COMPOUND_F,type_sizeI,memory_dtypes(i),error)
    CALL h5tinsert_f(memory_dtypes(i),"file format version",offset,H5T_NATIVE_INTEGER,error)

    i = i + 1
    call h5tcreate_f(H5T_COMPOUND_F,type_sizeS400,memory_dtypes(i),error)
    call h5tinsert_f(memory_dtypes(i),"setup call",offset,string_dtype400,error)

    i = i + 1
    call h5tcreate_f(H5T_COMPOUND_F,type_sizeS80,memory_dtypes(i),error)
    call h5tinsert_f(memory_dtypes(i),"file creation time",offset,string_dtype80,error)

    i = i + 1
    call h5tcreate_f(H5T_COMPOUND_F,type_sizeS80,memory_dtypes(i),error)
    call h5tinsert_f(memory_dtypes(i),"flash version",offset,string_dtype80,error)

    i = i + 1
    call h5tcreate_f(H5T_COMPOUND_F,type_sizeS80,memory_dtypes(i),error)
    call h5tinsert_f(memory_dtypes(i),"build date",offset,string_dtype80,error)

    i = i + 1
    call h5tcreate_f(H5T_COMPOUND_F,type_sizeS80,memory_dtypes(i),error)
    call h5tinsert_f(memory_dtypes(i),"build dir",offset,string_dtype80,error)

    i = i + 1
    call h5tcreate_f(H5T_COMPOUND_F,type_sizeS80,memory_dtypes(i),error)
    call h5tinsert_f(memory_dtypes(i),"build machine",offset,string_dtype80,error)

    i = i + 1
    call h5tcreate_f(H5T_COMPOUND_F,type_sizeS400,memory_dtypes(i),error)
    call h5tinsert_f(memory_dtypes(i),"cflags",offset,string_dtype400,error)

    i = i + 1
    call h5tcreate_f(H5T_COMPOUND_F,type_sizeS400,memory_dtypes(i),error)
    call h5tinsert_f(memory_dtypes(i),"fflags",offset,string_dtype400,error)

    i = i + 1
    call h5tcreate_f(H5T_COMPOUND_F,type_sizeS80,memory_dtypes(i),error)
    call h5tinsert_f(memory_dtypes(i),"setup time stamp",offset,string_dtype80,error)

    i = i + 1
    call h5tcreate_f(H5T_COMPOUND_F,type_sizeS80,memory_dtypes(i),error)
    call h5tinsert_f(memory_dtypes(i),"build time stamp",offset,string_dtype80,error)
# endif

# if 1
    i = 1
    CALL h5dwrite_f(dset_id,memory_dtypes(i),chkpt%file_format_version,(/INT(1,kind=SIZE_T)/),error,xfer_prp = chkpt%propertyid)

    i = i + 1
    CALL h5dwrite_f(dset_id,memory_dtypes(i),chkpt%setup_call,(/INT(1,kind=SIZE_T)/),error,xfer_prp = chkpt%propertyid)

    i = i + 1
    CALL h5dwrite_f(dset_id,memory_dtypes(i),chkpt%file_creation_time,(/INT(1,kind=SIZE_T)/),error,xfer_prp = chkpt%propertyid)

    i = i + 1
    CALL h5dwrite_f(dset_id,memory_dtypes(i),chkpt%flash_version,(/INT(1,kind=SIZE_T)/),error,xfer_prp = chkpt%propertyid)

    i = i + 1
    CALL h5dwrite_f(dset_id,memory_dtypes(i),chkpt%build_date,(/INT(1,kind=SIZE_T)/),error,xfer_prp = chkpt%propertyid)

    i = i + 1
    CALL h5dwrite_f(dset_id,memory_dtypes(i),chkpt%build_dir,(/INT(1,kind=SIZE_T)/),error,xfer_prp = chkpt%propertyid)

    i = i + 1
    CALL h5dwrite_f(dset_id,memory_dtypes(i),chkpt%build_machine,(/INT(1,kind=SIZE_T)/),error,xfer_prp = chkpt%propertyid)

    i = i + 1
    CALL h5dwrite_f(dset_id,memory_dtypes(i),chkpt%cflags,(/INT(1,kind=SIZE_T)/),error,xfer_prp = chkpt%propertyid)

    i = i + 1
    CALL h5dwrite_f(dset_id,memory_dtypes(i),chkpt%fflags,(/INT(1,kind=SIZE_T)/),error,xfer_prp = chkpt%propertyid)

    i = i + 1
    CALL h5dwrite_f(dset_id,memory_dtypes(i),chkpt%setup_time_stamp,(/INT(1,kind=SIZE_T)/),error,xfer_prp = chkpt%propertyid)

    i = i + 1
    CALL h5dwrite_f(dset_id,memory_dtypes(i),chkpt%build_time_stamp,(/INT(1,kind=SIZE_T)/),error,xfer_prp = chkpt%propertyid)
# endif

    !! End access to the dataset and release resources used by it.
    call h5dclose_f(dset_id,error)

    !! Terminate access to the data space.
    call h5sclose_f(dspace_id,error)
     
    !! Terminate access to the datatype
    CALL h5tclose_f(siminfo_dtype, error)

    CALL h5tclose_f(string_dtype80, error)
    CALL h5tclose_f(string_dtype400, error)

    do i = 1,size(memory_dtypes)
        if (memory_dtypes(i) > 0) CALL h5tclose_f(memory_dtypes(i), error)
    end do

end subroutine

subroutine checkpoint_write_integer_runtime_parameters(chkpt)

    use hdf5
    use checkpoint_types_mod, only: checkpoint_t

    type(checkpoint_t), intent(inout) :: chkpt

    INTEGER(HID_T) :: dset_id   ! Dset identifier
    INTEGER(HID_T) :: dspace_id ! Dataspace identifier
    INTEGER(HID_T) :: dtype_id  ! Compound datatype identifier

    INTEGER(HID_T) :: dt1_id      ! Memory datatype identifier (for character field)
    INTEGER(HID_T) :: dt2_id      ! Memory datatype identifier (for integer field)
    INTEGER(HID_T) :: dt5_id      ! Memory datatype identifier 

    INTEGER             ::   error      ! Error flag
    INTEGER(SIZE_T)     ::   type_sizec ! Size of compound datatype
    INTEGER(SIZE_T)     ::   type_sizek ! Size of name datatype 
    INTEGER(SIZE_T)     ::   type_sizev ! Size of value datatype
    INTEGER(SIZE_T)     ::   offset     ! Member's offset

    !! Calculate total size by calculating sizes of each member
    call h5tcopy_f(H5T_C_S1,dt5_id,error)

    call h5tset_size_f(dt5_id,INT(80,kind=SIZE_T),error)
    call h5tget_size_f(dt5_id,type_sizek,error)
    call h5tget_size_f(H5T_NATIVE_INTEGER,type_sizev,error)

    type_sizec = type_sizek + type_sizev

    CALL h5tcreate_f(H5T_COMPOUND_F,type_sizec,dtype_id,error)

    !! string memeber
    offset = 0
    call h5tinsert_f(dtype_id,"name",offset,dt5_id,error)
     
    !! integer member
    offset = offset + type_sizek
    call h5tinsert_f(dtype_id,"value",offset,H5T_NATIVE_INTEGER,error)

    !! data space
    call h5screate_simple_f(1,(/INT(chkpt%index_integer_runtime_parameters,kind=SIZE_T)/),dspace_id,error)

    !! dataset
    call h5dcreate_f(chkpt%fileid,"integer runtime parameters",dtype_id,dspace_id,dset_id,error)

    !! Create memory types.
    call h5tcreate_f(H5T_COMPOUND_F,type_sizeK,dt1_id,error)
    offset = 0
    CALL h5tinsert_f(dt1_id,"name",offset,dt5_id,error)

    call h5tcreate_f(H5T_COMPOUND_F,type_sizeV,dt2_id,error)
    offset = 0
    call h5tinsert_f(dt2_id,"value",offset,H5T_NATIVE_INTEGER,error)

    !! Write data by fields in the datatype. Fields order is not important.
    CALL h5dwrite_f(dset_id,dt1_id,chkpt%keys_integer_runtime_parameters,&
        (/INT(chkpt%index_integer_runtime_parameters,kind=SIZE_T)/),error,xfer_prp = chkpt%propertyid)
    CALL h5dwrite_f(dset_id,dt2_id,chkpt%vals_integer_runtime_parameters,&
        (/INT(chkpt%index_integer_runtime_parameters,kind=SIZE_T)/),error,xfer_prp = chkpt%propertyid)

    !! End access to the dataset and release resources used by it.
    call h5dclose_f(dset_id,error)

    !! Terminate access to the data space.
    call h5sclose_f(dspace_id,error)
     
    !! Terminate access to the datatype
    CALL h5tclose_f(dtype_id, error)
    CALL h5tclose_f(dt1_id, error)
    CALL h5tclose_f(dt2_id, error)
    CALL h5tclose_f(dt5_id, error)

end subroutine

subroutine checkpoint_write_real_runtime_parameters(chkpt)

    use hdf5
    use checkpoint_types_mod, only: checkpoint_t

    type(checkpoint_t), intent(inout) :: chkpt

    INTEGER(HID_T)  :: dset_id   ! Dset identifier
    INTEGER(HID_T)  :: dspace_id ! Dataspace identifier
    INTEGER(HID_T)  :: dtype_id  ! Compound datatype identifier

    INTEGER(HID_T)  :: dt1_id      ! Memory datatype identifier (for character field)
    INTEGER(HID_T)  :: dt2_id      ! Memory datatype identifier (for integer field)
    INTEGER(HID_T)  :: dt5_id      ! Memory datatype identifier 

    INTEGER(SIZE_T) :: type_sizec ! Size of compound datatype
    INTEGER(SIZE_T) :: type_sizek ! Size of name datatype 
    INTEGER(SIZE_T) :: type_sizev ! Size of value datatype
    INTEGER(SIZE_T) :: offset     ! Member's offset

    INTEGER         :: error      ! Error flag

    !! Calculate total size by calculating sizes of each member
    call h5tcopy_f(H5T_C_S1,dt5_id,error)

    call h5tset_size_f(dt5_id,INT(80,kind=SIZE_T),error)
    call h5tget_size_f(dt5_id,type_sizek,error)
    call h5tget_size_f(H5T_NATIVE_DOUBLE,type_sizev,error)

    type_sizec = type_sizek + type_sizev

    CALL h5tcreate_f(H5T_COMPOUND_F,type_sizec,dtype_id,error)

    !! string memeber
    offset = 0
    call h5tinsert_f(dtype_id,"name",offset,dt5_id,error)
     
    !! integer member
    offset = offset + type_sizek
    call h5tinsert_f(dtype_id,"value",offset,H5T_NATIVE_DOUBLE,error)

    !! data space
    call h5screate_simple_f(1,(/INT(chkpt%index_real_runtime_parameters,kind=SIZE_T)/),dspace_id,error)

    !! dataset
    call h5dcreate_f(chkpt%fileid,"real runtime parameters",dtype_id,dspace_id,dset_id,error)

    !! Create memory types.
    call h5tcreate_f(H5T_COMPOUND_F,type_sizek,dt1_id,error)
    offset = 0
    CALL h5tinsert_f(dt1_id,"name",offset,dt5_id,error)

    call h5tcreate_f(H5T_COMPOUND_F,type_sizev,dt2_id,error)
    offset = 0
    call h5tinsert_f(dt2_id,"value",offset,H5T_NATIVE_DOUBLE,error)

    !! Write data by fields in the datatype. Fields order is not important.
    CALL h5dwrite_f(dset_id,dt1_id,chkpt%keys_real_runtime_parameters,&
        (/INT(chkpt%index_real_runtime_parameters,kind=SIZE_T)/),error,xfer_prp = chkpt%propertyid)
    CALL h5dwrite_f(dset_id,dt2_id,chkpt%vals_real_runtime_parameters,&
        (/INT(chkpt%index_real_runtime_parameters,kind=SIZE_T)/),error,xfer_prp = chkpt%propertyid)

    !! End access to the dataset and release resources used by it.
    call h5dclose_f(dset_id,error)

    !! Terminate access to the data space.
    call h5sclose_f(dspace_id,error)
     
    !! Terminate access to the datatype
    CALL h5tclose_f(dtype_id, error)
    CALL h5tclose_f(dt1_id, error)
    CALL h5tclose_f(dt2_id, error)
    CALL h5tclose_f(dt5_id, error)

end subroutine

subroutine checkpoint_write_logical_runtime_parameters(chkpt)

    use hdf5
    use checkpoint_types_mod, only: checkpoint_t

    type(checkpoint_t), intent(inout) :: chkpt

    INTEGER(HID_T)  :: dset_id   ! Dset identifier
    INTEGER(HID_T)  :: dspace_id ! Dataspace identifier
    INTEGER(HID_T)  :: dtype_id  ! Compound datatype identifier

    INTEGER(HID_T)  :: dt1_id      ! Memory datatype identifier (for character field)
    INTEGER(HID_T)  :: dt2_id      ! Memory datatype identifier (for integer field)
    INTEGER(HID_T)  :: dt5_id      ! Memory datatype identifier 

    INTEGER(SIZE_T) :: type_sizec ! Size of compound datatype
    INTEGER(SIZE_T) :: type_sizek ! Size of name datatype 
    INTEGER(SIZE_T) :: type_sizev ! Size of value datatype
    INTEGER(SIZE_T) :: offset     ! Member's offset

    INTEGER         :: error      ! Error flag

    !! Calculate total size by calculating sizes of each member
    call h5tcopy_f(H5T_C_S1,dt5_id,error)

    call h5tset_size_f(dt5_id,INT(80,kind=SIZE_T),error)
    call h5tget_size_f(dt5_id,type_sizek,error)
    call h5tget_size_f(H5T_NATIVE_INTEGER,type_sizev,error)

    type_sizec = type_sizek + type_sizev

    CALL h5tcreate_f(H5T_COMPOUND_F,type_sizec,dtype_id,error)

    !! string memeber
    offset = 0
    call h5tinsert_f(dtype_id,"name",offset,dt5_id,error)
     
    !! integer member
    offset = offset + type_sizek
    call h5tinsert_f(dtype_id,"value",offset,H5T_NATIVE_INTEGER,error)

    !! data space
    call h5screate_simple_f(1,(/INT(chkpt%index_logical_runtime_parameters,kind=SIZE_T)/),dspace_id,error)

    !! dataset
    call h5dcreate_f(chkpt%fileid,"logical runtime parameters",dtype_id,dspace_id,dset_id,error)

    !! Create memory types.
    call h5tcreate_f(H5T_COMPOUND_F,type_sizek,dt1_id,error)
    offset = 0
    CALL h5tinsert_f(dt1_id,"name",offset,dt5_id,error)

    call h5tcreate_f(H5T_COMPOUND_F,type_sizev,dt2_id,error)
    offset = 0
    call h5tinsert_f(dt2_id,"value",offset,H5T_NATIVE_INTEGER,error)

    !! Write data by fields in the datatype. Fields order is not important.
    CALL h5dwrite_f(dset_id,dt1_id,chkpt%keys_logical_runtime_parameters,&
        (/INT(chkpt%index_logical_runtime_parameters,kind=SIZE_T)/),error,xfer_prp = chkpt%propertyid)
    CALL h5dwrite_f(dset_id,dt2_id,merge(1,0,chkpt%vals_logical_runtime_parameters),&
        (/INT(chkpt%index_logical_runtime_parameters,kind=SIZE_T)/),error,xfer_prp = chkpt%propertyid)

    !! End access to the dataset and release resources used by it.
    call h5dclose_f(dset_id,error)

    !! Terminate access to the data space.
    call h5sclose_f(dspace_id,error)
     
    !! Terminate access to the datatype
    CALL h5tclose_f(dtype_id, error)
    CALL h5tclose_f(dt1_id, error)
    CALL h5tclose_f(dt2_id, error)
    CALL h5tclose_f(dt5_id, error)

end subroutine

subroutine checkpoint_write_string_runtime_parameters(chkpt)

    use hdf5
    use checkpoint_types_mod, only: checkpoint_t

    type(checkpoint_t), intent(inout) :: chkpt

    INTEGER(HID_T)  :: dset_id   ! Dset identifier
    INTEGER(HID_T)  :: dspace_id ! Dataspace identifier
    INTEGER(HID_T)  :: dtype_id  ! Compound datatype identifier

    INTEGER(HID_T)  :: dt1_id      ! Memory datatype identifier (for character field)
    INTEGER(HID_T)  :: dt2_id      ! Memory datatype identifier (for integer field)
    INTEGER(HID_T)  :: dt5_id      ! Memory datatype identifier 

    INTEGER(SIZE_T) :: type_sizec ! Size of compound datatype
    INTEGER(SIZE_T) :: type_sizek ! Size of name datatype 
    INTEGER(SIZE_T) :: type_sizev ! Size of value datatype
    INTEGER(SIZE_T) :: offset     ! Member's offset

    INTEGER         :: error      ! Error flag

    !! Calculate total size by calculating sizes of each member
    call h5tcopy_f(H5T_C_S1,dt5_id,error)

    call h5tset_size_f(dt5_id,INT(80,kind=SIZE_T),error)
    call h5tget_size_f(dt5_id,type_sizek,error)
    call h5tget_size_f(dt5_id,type_sizev,error)

    type_sizec = type_sizek + type_sizev

    CALL h5tcreate_f(H5T_COMPOUND_F,type_sizec,dtype_id,error)

    !! string memeber
    offset = 0
    call h5tinsert_f(dtype_id,"name",offset,dt5_id,error)
     
    !! integer member
    offset = offset + type_sizek
    call h5tinsert_f(dtype_id,"value",offset,dt5_id,error)

    !! data space
    call h5screate_simple_f(1,(/INT(chkpt%index_string_runtime_parameters,kind=SIZE_T)/),dspace_id,error)

    !! dataset
    call h5dcreate_f(chkpt%fileid,"string runtime parameters",dtype_id,dspace_id,dset_id,error)

    !! Create memory types.
    call h5tcreate_f(H5T_COMPOUND_F,type_sizek,dt1_id,error)
    offset = 0
    CALL h5tinsert_f(dt1_id,"name",offset,dt5_id,error)

    call h5tcreate_f(H5T_COMPOUND_F,type_sizev,dt2_id,error)
    offset = 0
    call h5tinsert_f(dt2_id,"value",offset,dt5_id,error)

    !! Write data by fields in the datatype. Fields order is not important.
    CALL h5dwrite_f(dset_id,dt1_id,chkpt%keys_string_runtime_parameters,&
        (/INT(chkpt%index_string_runtime_parameters,kind=SIZE_T)/),error,xfer_prp = chkpt%propertyid)
    CALL h5dwrite_f(dset_id,dt2_id,chkpt%vals_string_runtime_parameters,&
        (/INT(chkpt%index_string_runtime_parameters,kind=SIZE_T)/),error,xfer_prp = chkpt%propertyid)

    !! End access to the dataset and release resources used by it.
    call h5dclose_f(dset_id,error)

    !! Terminate access to the data space.
    call h5sclose_f(dspace_id,error)
     
    !! Terminate access to the datatype
    CALL h5tclose_f(dtype_id, error)
    CALL h5tclose_f(dt1_id, error)
    CALL h5tclose_f(dt2_id, error)
    CALL h5tclose_f(dt5_id, error)

end subroutine

subroutine checkpoint_write_integer_scalars(chkpt)

    use hdf5
    use checkpoint_types_mod, only: checkpoint_t

    type(checkpoint_t), intent(inout) :: chkpt

    INTEGER(HID_T) :: dset_id   ! Dset identifier
    INTEGER(HID_T) :: dspace_id ! Dataspace identifier
    INTEGER(HID_T) :: dtype_id  ! Compound datatype identifier

    INTEGER(HID_T) :: dt1_id      ! Memory datatype identifier (for character field)
    INTEGER(HID_T) :: dt2_id      ! Memory datatype identifier (for integer field)
    INTEGER(HID_T) :: dt5_id      ! Memory datatype identifier 

    INTEGER             ::   error      ! Error flag
    INTEGER(SIZE_T)     ::   type_sizec ! Size of compound datatype
    INTEGER(SIZE_T)     ::   type_sizek ! Size of name datatype 
    INTEGER(SIZE_T)     ::   type_sizev ! Size of value datatype
    INTEGER(SIZE_T)     ::   offset     ! Member's offset

    !! Calculate total size by calculating sizes of each member
    call h5tcopy_f(H5T_C_S1,dt5_id,error)

    call h5tset_size_f(dt5_id,INT(80,kind=SIZE_T),error)
    call h5tget_size_f(dt5_id,type_sizek,error)
    call h5tget_size_f(H5T_NATIVE_INTEGER,type_sizev,error)

    type_sizec = type_sizek + type_sizev

    CALL h5tcreate_f(H5T_COMPOUND_F,type_sizec,dtype_id,error)

    !! string memeber
    offset = 0
    call h5tinsert_f(dtype_id,"name",offset,dt5_id,error)
     
    !! integer member
    offset = offset + type_sizek
    call h5tinsert_f(dtype_id,"value",offset,H5T_NATIVE_INTEGER,error)

    !! data space
    call h5screate_simple_f(1,(/INT(chkpt%index_integer_scalars,kind=SIZE_T)/),dspace_id,error)

    !! dataset
    call h5dcreate_f(chkpt%fileid,"integer scalars",dtype_id,dspace_id,dset_id,error)

    !! Create memory types.
    call h5tcreate_f(H5T_COMPOUND_F,type_sizeK,dt1_id,error)
    offset = 0
    CALL h5tinsert_f(dt1_id,"name",offset,dt5_id,error)

    call h5tcreate_f(H5T_COMPOUND_F,type_sizeV,dt2_id,error)
    offset = 0
    call h5tinsert_f(dt2_id,"value",offset,H5T_NATIVE_INTEGER,error)

    !! Write data by fields in the datatype. Fields order is not important.
    CALL h5dwrite_f(dset_id,dt1_id,chkpt%keys_integer_scalars,&
        (/INT(chkpt%index_integer_scalars,kind=SIZE_T)/),error,xfer_prp = chkpt%propertyid)
    CALL h5dwrite_f(dset_id,dt2_id,chkpt%vals_integer_scalars,&
        (/INT(chkpt%index_integer_scalars,kind=SIZE_T)/),error,xfer_prp = chkpt%propertyid)

    !! End access to the dataset and release resources used by it.
    call h5dclose_f(dset_id,error)

    !! Terminate access to the data space.
    call h5sclose_f(dspace_id,error)
     
    !! Terminate access to the datatype
    CALL h5tclose_f(dtype_id, error)
    CALL h5tclose_f(dt1_id, error)
    CALL h5tclose_f(dt2_id, error)
    CALL h5tclose_f(dt5_id, error)

end subroutine

subroutine checkpoint_write_real_scalars(chkpt)

    use hdf5
    use checkpoint_types_mod, only: checkpoint_t

    type(checkpoint_t), intent(inout) :: chkpt

    INTEGER(HID_T)  :: dset_id   ! Dset identifier
    INTEGER(HID_T)  :: dspace_id ! Dataspace identifier
    INTEGER(HID_T)  :: dtype_id  ! Compound datatype identifier

    INTEGER(HID_T)  :: dt1_id      ! Memory datatype identifier (for character field)
    INTEGER(HID_T)  :: dt2_id      ! Memory datatype identifier (for integer field)
    INTEGER(HID_T)  :: dt5_id      ! Memory datatype identifier 

    INTEGER(SIZE_T) :: type_sizec ! Size of compound datatype
    INTEGER(SIZE_T) :: type_sizek ! Size of name datatype 
    INTEGER(SIZE_T) :: type_sizev ! Size of value datatype
    INTEGER(SIZE_T) :: offset     ! Member's offset

    INTEGER         :: error      ! Error flag

    !! Calculate total size by calculating sizes of each member
    call h5tcopy_f(H5T_C_S1,dt5_id,error)

    !call h5tset_size_f(dt5_id,INT(kind=SIZE_T),error)
    call h5tset_size_f(dt5_id,INT(80,kind=SIZE_T),error)
    call h5tget_size_f(dt5_id,type_sizek,error)
    call h5tget_size_f(H5T_NATIVE_DOUBLE,type_sizev,error)

    type_sizec = type_sizek + type_sizev

    CALL h5tcreate_f(H5T_COMPOUND_F,type_sizec,dtype_id,error)

    !! string memeber
    offset = 0
    call h5tinsert_f(dtype_id,"name",offset,dt5_id,error)
     
    !! integer member
    offset = offset + type_sizek
    call h5tinsert_f(dtype_id,"value",offset,H5T_NATIVE_DOUBLE,error)

    !! data space
    call h5screate_simple_f(1,(/INT(chkpt%index_real_scalars,kind=SIZE_T)/),dspace_id,error)

    !! dataset
    call h5dcreate_f(chkpt%fileid,"real scalars",dtype_id,dspace_id,dset_id,error)

    !! Create memory types.
    call h5tcreate_f(H5T_COMPOUND_F,type_sizek,dt1_id,error)
    offset = 0
    CALL h5tinsert_f(dt1_id,"name",offset,dt5_id,error)

    call h5tcreate_f(H5T_COMPOUND_F,type_sizev,dt2_id,error)
    offset = 0
    call h5tinsert_f(dt2_id,"value",offset,H5T_NATIVE_DOUBLE,error)

    !! Write data by fields in the datatype. Fields order is not important.
    CALL h5dwrite_f(dset_id,dt1_id,chkpt%keys_real_scalars,&
        (/INT(chkpt%index_real_scalars,kind=SIZE_T)/),error,xfer_prp = chkpt%propertyid)
    CALL h5dwrite_f(dset_id,dt2_id,chkpt%vals_real_scalars,&
        (/INT(chkpt%index_real_scalars,kind=SIZE_T)/),error,xfer_prp = chkpt%propertyid)

    !! End access to the dataset and release resources used by it.
    call h5dclose_f(dset_id,error)

    !! Terminate access to the data space.
    call h5sclose_f(dspace_id,error)
     
    !! Terminate access to the datatype
    CALL h5tclose_f(dtype_id, error)
    CALL h5tclose_f(dt1_id, error)
    CALL h5tclose_f(dt2_id, error)
    CALL h5tclose_f(dt5_id, error)

end subroutine

subroutine checkpoint_write_logical_scalars(chkpt)

    use hdf5
    use checkpoint_types_mod, only: checkpoint_t

    type(checkpoint_t), intent(inout) :: chkpt

    INTEGER(HID_T)  :: dset_id   ! Dset identifier
    INTEGER(HID_T)  :: dspace_id ! Dataspace identifier
    INTEGER(HID_T)  :: dtype_id  ! Compound datatype identifier

    INTEGER(HID_T)  :: dt1_id      ! Memory datatype identifier (for character field)
    INTEGER(HID_T)  :: dt2_id      ! Memory datatype identifier (for integer field)
    INTEGER(HID_T)  :: dt5_id      ! Memory datatype identifier 

    INTEGER(SIZE_T) :: type_sizec ! Size of compound datatype
    INTEGER(SIZE_T) :: type_sizek ! Size of name datatype 
    INTEGER(SIZE_T) :: type_sizev ! Size of value datatype
    INTEGER(SIZE_T) :: offset     ! Member's offset

    INTEGER         :: error      ! Error flag

    !! Calculate total size by calculating sizes of each member
    call h5tcopy_f(H5T_C_S1,dt5_id,error)

    !call h5tset_size_f(dt5_id,INT(kind=SIZE_T),error)
    call h5tset_size_f(dt5_id,INT(80,kind=SIZE_T),error)
    call h5tget_size_f(dt5_id,type_sizek,error)
    call h5tget_size_f(H5T_NATIVE_INTEGER,type_sizev,error)

    type_sizec = type_sizek + type_sizev

    CALL h5tcreate_f(H5T_COMPOUND_F,type_sizec,dtype_id,error)

    !! string memeber
    offset = 0
    call h5tinsert_f(dtype_id,"name",offset,dt5_id,error)
     
    !! integer member
    offset = offset + type_sizek
    call h5tinsert_f(dtype_id,"value",offset,H5T_NATIVE_INTEGER,error)

    !! data space
    call h5screate_simple_f(1,(/INT(chkpt%index_logical_scalars,kind=SIZE_T)/),dspace_id,error)

    !! dataset
    call h5dcreate_f(chkpt%fileid,"logical scalars",dtype_id,dspace_id,dset_id,error)

    !! Create memory types.
    call h5tcreate_f(H5T_COMPOUND_F,type_sizek,dt1_id,error)
    offset = 0
    CALL h5tinsert_f(dt1_id,"name",offset,dt5_id,error)

    call h5tcreate_f(H5T_COMPOUND_F,type_sizev,dt2_id,error)
    offset = 0
    call h5tinsert_f(dt2_id,"value",offset,H5T_NATIVE_INTEGER,error)

    !! Write data by fields in the datatype. Fields order is not important.
    CALL h5dwrite_f(dset_id,dt1_id,chkpt%keys_logical_scalars,&
        (/INT(chkpt%index_logical_scalars,kind=SIZE_T)/),error,xfer_prp = chkpt%propertyid)
    CALL h5dwrite_f(dset_id,dt2_id,merge(1,0,chkpt%vals_logical_scalars),&
        (/INT(chkpt%index_logical_scalars,kind=SIZE_T)/),error,xfer_prp = chkpt%propertyid)

    !! End access to the dataset and release resources used by it.
    call h5dclose_f(dset_id,error)

    !! Terminate access to the data space.
    call h5sclose_f(dspace_id,error)
     
    !! Terminate access to the datatype
    CALL h5tclose_f(dtype_id, error)
    CALL h5tclose_f(dt1_id, error)
    CALL h5tclose_f(dt2_id, error)
    CALL h5tclose_f(dt5_id, error)

end subroutine

subroutine checkpoint_write_string_scalars(chkpt)

    use hdf5
    use checkpoint_types_mod, only: checkpoint_t

    type(checkpoint_t), intent(inout) :: chkpt

    INTEGER(HID_T)  :: dset_id   ! Dset identifier
    INTEGER(HID_T)  :: dspace_id ! Dataspace identifier
    INTEGER(HID_T)  :: dtype_id  ! Compound datatype identifier

    INTEGER(HID_T)  :: dt1_id      ! Memory datatype identifier (for character field)
    INTEGER(HID_T)  :: dt2_id      ! Memory datatype identifier (for integer field)
    INTEGER(HID_T)  :: dt5_id      ! Memory datatype identifier 

    INTEGER(SIZE_T) :: type_sizec ! Size of compound datatype
    INTEGER(SIZE_T) :: type_sizek ! Size of name datatype 
    INTEGER(SIZE_T) :: type_sizev ! Size of value datatype
    INTEGER(SIZE_T) :: offset     ! Member's offset

    INTEGER         :: error      ! Error flag

    !! Calculate total size by calculating sizes of each member
    call h5tcopy_f(H5T_C_S1,dt5_id,error)

    call h5tset_size_f(dt5_id,INT(80,kind=SIZE_T),error)
    call h5tget_size_f(dt5_id,type_sizek,error)
    call h5tget_size_f(dt5_id,type_sizev,error)

    type_sizec = type_sizek + type_sizev

    CALL h5tcreate_f(H5T_COMPOUND_F,type_sizec,dtype_id,error)

    !! string memeber
    offset = 0
    call h5tinsert_f(dtype_id,"name",offset,dt5_id,error)
     
    !! integer member
    offset = offset + type_sizek
    call h5tinsert_f(dtype_id,"value",offset,dt5_id,error)

    !! data space
    call h5screate_simple_f(1,(/INT(chkpt%index_string_scalars,kind=SIZE_T)/),dspace_id,error)

    !! dataset
    call h5dcreate_f(chkpt%fileid,"string scalars",dtype_id,dspace_id,dset_id,error)

    !! Create memory types.
    call h5tcreate_f(H5T_COMPOUND_F,type_sizek,dt1_id,error)
    offset = 0
    CALL h5tinsert_f(dt1_id,"name",offset,dt5_id,error)

    call h5tcreate_f(H5T_COMPOUND_F,type_sizev,dt2_id,error)
    offset = 0
    call h5tinsert_f(dt2_id,"value",offset,dt5_id,error)

    !! Write data by fields in the datatype. Fields order is not important.
    CALL h5dwrite_f(dset_id,dt1_id,chkpt%keys_string_scalars,&
        (/INT(chkpt%index_string_scalars,kind=SIZE_T)/),error,xfer_prp = chkpt%propertyid)
    CALL h5dwrite_f(dset_id,dt2_id,chkpt%vals_string_scalars,&
        (/INT(chkpt%index_string_scalars,kind=SIZE_T)/),error,xfer_prp = chkpt%propertyid)

    !! End access to the dataset and release resources used by it.
    call h5dclose_f(dset_id,error)

    !! Terminate access to the data space.
    call h5sclose_f(dspace_id,error)
     
    !! Terminate access to the datatype
    CALL h5tclose_f(dtype_id, error)
    CALL h5tclose_f(dt1_id, error)
    CALL h5tclose_f(dt2_id, error)
    CALL h5tclose_f(dt5_id, error)

end subroutine

subroutine checkpoint_write_unknown_names(chkpt)

    use hdf5
    use checkpoint_types_mod, only: checkpoint_t

    type(checkpoint_t), intent(inout) :: chkpt

    INTEGER(HID_T)  :: dset_id   ! Dset identifier
    INTEGER(HID_T)  :: dspace_id ! Dataspace identifier
    INTEGER(HID_T)  :: dtype_id  ! Compound datatype identifier

    INTEGER(HID_T)  :: dt1_id      ! Memory datatype identifier (for character field)
    INTEGER(HID_T)  :: dt2_id      ! Memory datatype identifier (for integer field)
    INTEGER(HID_T)  :: dt5_id      ! Memory datatype identifier 

    INTEGER(SIZE_T) :: type_sizec ! Size of compound datatype
    INTEGER(SIZE_T) :: type_sizek ! Size of name datatype 
    INTEGER(SIZE_T) :: type_sizev ! Size of value datatype
    INTEGER(SIZE_T) :: offset     ! Member's offset

    INTEGER         :: error      ! Error flag

    !! Calculate total size by calculating sizes of each member
    call h5tcopy_f(H5T_C_S1,dt5_id,error)

    call h5tset_size_f(dt5_id,INT(4,kind=SIZE_T),error)
    call h5tget_size_f(dt5_id,type_sizek,error)

    !! data space
    call h5screate_simple_f(2,(/INT(1,kind=SIZE_T),INT(chkpt%index_unknown_names,kind=SIZE_T)/),dspace_id,error)

    !! dataset
    call h5dcreate_f(chkpt%fileid,"unknown names",dt5_id,dspace_id,dset_id,error)

    !! Write data by fields in the datatype. Fields order is not important.
    CALL h5dwrite_f(dset_id,dt5_id,chkpt%unknown_names,&
        (/INT(chkpt%index_unknown_names,kind=SIZE_T)/),error,xfer_prp = chkpt%propertyid)

    !! End access to the dataset and release resources used by it.
    call h5dclose_f(dset_id,error)

    !! Terminate access to the data space.
    call h5sclose_f(dspace_id,error)
     
    !! Terminate access to the datatype
    CALL h5tclose_f(dt5_id, error)

end subroutine

end module
