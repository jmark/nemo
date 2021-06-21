# define tint integer
# define bool logical
# define char character
# define thid integer(HID_T)
# define thiz integer(HSIZE_T)

module h5simple_mod

use hdf5 !! pmake: ignore
use constants_mod

integer, parameter :: h5simple_dtype_int8 = 1
integer, parameter :: h5simple_dtype_int32 = 2
integer, parameter :: h5simple_dtype_float32 = 3
integer, parameter :: h5simple_dtype_float64 = 4

interface h5simple_create_path
    procedure h5simple_create_path_string
    !procedure h5simple_create_path_optionals
end interface

interface h5simple_make_dset
    procedure h5simple_make_dset_integr
    procedure h5simple_make_dset_double
    procedure h5simple_make_dset_string
end interface

interface h5simple_make_attr
    procedure h5simple_make_attr_integr
    procedure h5simple_make_attr_string
end interface

interface h5simple_read_attr
    procedure h5simple_read_attr_string
    procedure h5simple_read_attr_integr
    procedure h5simple_read_attr_double
end interface

interface h5simple_fill_dset
    procedure h5simple_fill_dset_char
    procedure h5simple_fill_dset_double
    procedure h5simple_fill_dset_integr
end interface

interface h5simple_read_dset
    procedure h5simple_read_dset_string
    procedure h5simple_read_dset_integr
    procedure h5simple_read_dset_double
end interface

contains

function h5simple_path_split(path, components, isdir, errstring) result(ncomps)

    character(len=*), intent(in)    :: path
    character(len=*), intent(out)   :: components(:)
    logical, intent(out)            :: isdir
    character(len=*), intent(out)   :: errstring

    integer :: ncomps

    integer :: i,a,b,pathlen

    isdir = .false.
    ncomps = 0
    components = ""

    if (path(1:1) /= '/') then
        errstring = "does not start with '/'"
        ncomps = -1
        return
    end if

    pathlen = len_trim(path)

    if (path(pathlen:pathlen) == '/') then
        isdir = .true.
        pathlen = pathlen - 1
    end if

    a = 1
    b = 1
    ncomps = 0

    do i = 2,pathlen
        if (path(i:i) == '/') then
            b = i
            if (a+1 > b-1) then
                errstring = "invalid path"
                ncomps = -1
                return
            end if
            ncomps = ncomps + 1
            if (ncomps >= size(components)) then
                errstring = "components array to small"
                ncomps = -1
                return
            end if
            components(ncomps) = trim(path(a+1:b-1))
            a = b
        end if
    end do
            
    ncomps = ncomps + 1
    components(ncomps) = trim(path(a+1:pathlen))

end function

function h5simple_path_exists(fileid, path) result(retval)

    thid, intent(in)        :: fileid
    char(len=*), intent(in) :: path
    tint                    :: stat
    logical                 :: retval

    retval = .false.
    call h5lexists_f(fileid, path, retval, stat)

end function

function h5simple_create_path_string(locid, path) result(ok)

    thid, intent(in)        :: locid 
    char(len=*), intent(in) :: path
    logical                 :: ok

    integer(HID_T)          :: gids(0:16)

    character(len=64) :: components(16), errstring
    integer :: ncomponents, i, stat
    logical :: isdir

    stat = 0

    ncomponents = h5simple_path_split(path,components,isdir,errstring)    

    if (ncomponents < 1) then
        write (*,*) errstring, path
        ok = .false.
        return
    end if

    !write (*,*) 'isdir', isdir
    !do i = 1,ncomponents
    !    write (*,*) i, components(i)
    !end do

    if (.not.isdir) ncomponents = ncomponents - 1

    gids(0) = locid

    do i = 1,ncomponents
        if (.not.h5simple_path_exists(gids(i-1), trim(components(i)))) then
            call h5gcreate_f(gids(i-1), trim(components(i)), gids(i), stat)
        else
            call h5gopen_f(gids(i-1), trim(components(i)), gids(i), stat)
        end if
    end do

    do i = 1,ncomponents
        call h5gclose_f(gids(i), stat)
    end do

    ok = stat == 0
  
end function

# if 0
function h5simple_create_path_optionals(locid, path1, path2, path3, path4, path5, path6, path7) result(stat)

    thid, intent(in)                    :: locid 
    char(len=*), intent(in), optional   :: path1
    char(len=*), intent(in), optional   :: path2
    char(len=*), intent(in), optional   :: path3
    char(len=*), intent(in), optional   :: path4
    char(len=*), intent(in), optional   :: path5
    char(len=*), intent(in), optional   :: path6
    char(len=*), intent(in), optional   :: path7

    tint                    :: stat
    integer(HID_T)          :: gids(0:10)
    tint                    :: i

    stat = 0
    gids = -1
    gids(0) = locid

# define PP_CREATE(PNAME, IDX) \
    if (present(PNAME)) then; \
        if (.not.h5simple_path_exists(gids(IDX-1), PNAME)) then; \
            call h5gcreate_f(gids(IDX-1), PNAME, gids(IDX), stat); \
        else; \
            call h5gopen_f(gids(IDX-1), PNAME, gids(IDX), stat); \
        end if; \
    end if

    PP_CREATE(path1, 1)
    PP_CREATE(path2, 2)
    PP_CREATE(path3, 3)
    PP_CREATE(path4, 4)
    PP_CREATE(path5, 5)
    PP_CREATE(path6, 6)
    PP_CREATE(path7, 7)

# undef PP_APPEND

    do i = 1, size(gids)-1
        if (gids(i) > 0) call h5gclose_f(gids(i), stat)
    end do
  
end function
# endif

function h5simple_create_dataset(fileid, dname, typeid, dshape, propertyid) result(ok)

    char(len=*), intent(in)     :: dname
    thid, intent(in)            :: fileid
    thid, intent(in)            :: typeid
    integer(kind=8), intent(in) :: dshape(:)
    thid, intent(in), optional  :: propertyid

    thiz, allocatable       :: shapeinfo(:)
    thid                    :: dsetid, spaceid
    tint                    :: ierr
    logical                 :: ok

    allocate(shapeinfo(SIZE(dshape)))
    shapeinfo = dshape

    call h5screate_simple_f(size(shapeinfo), shapeinfo, spaceid, ierr)

    if (PRESENT(propertyid)) then
        call h5dcreate_f(fileid, trim(dname), typeid, spaceid, dsetid, ierr, dcpl_id=propertyid)
    else
        call h5dcreate_f(fileid, trim(dname), typeid, spaceid, dsetid, ierr)
    end if

    call h5sclose_f(spaceid, ierr)
    call h5dclose_f(dsetid, ierr)

    ok = (ierr == 0)
    
end function

function h5simple_prep_dset(fileid, dname, typeid, dshape, propertyid, doclose) result(dsetid)

    char(len=*), intent(in) :: dname
    thid, intent(in)        :: fileid
    thid, intent(in)        :: typeid
    tint, intent(in)        :: dshape(:)

    thid, intent(in), optional  :: propertyid
    bool, intent(in), optional  :: doclose

    thiz, allocatable       :: shapeinfo(:)
    thid                    :: dsetid, spaceid
    tint                    :: ierr

    allocate(shapeinfo(SIZE(dshape)))
    shapeinfo = dshape

    call h5screate_simple_f(size(shapeinfo), shapeinfo, spaceid, ierr)

    if (PRESENT(propertyid)) then
        call h5dcreate_f(fileid, trim(dname), typeid, spaceid, dsetid, ierr, dcpl_id=propertyid)
    else
        call h5dcreate_f(fileid, trim(dname), typeid, spaceid, dsetid, ierr)
    end if

    call h5sclose_f(spaceid, ierr)
    if (PRESENT(doclose)) then
        if (doclose) then
            call h5dclose_f(dsetid, ierr)
        end if
    end if

end function

function h5simple_make_dset_string(fileid, dname, dvalue, dshape, strlen, propertyid) result(ok)

    character(len=*), intent(in)    :: dname
    integer(HID_T), intent(in)      :: fileid
    integer, intent(in)             :: strlen
    character(len=strlen)           :: dvalue(*)
    integer, intent(in)             :: dshape(:)

    integer(HID_T), intent(in), optional    :: propertyid

    integer(HSIZE_T), allocatable   :: shapeinfo(:)
    integer(HID_T)                  :: dsetid, strid

    integer(SIZE_T)                 :: tstrlen
    integer                         :: ierr

    logical                         :: ok

    tstrlen = strlen

    call h5tcopy_f(H5T_NATIVE_CHARACTER, strid, ierr)
    call h5tset_size_f(strid, tstrlen, ierr)

    allocate(shapeinfo(SIZE(dshape)))
    shapeinfo = dshape

    if (PRESENT(propertyid)) then
        dsetid = h5simple_prep_dset(fileid, dname, strid, dshape, propertyid=propertyid)
    else                                    
        dsetid = h5simple_prep_dset(fileid, dname, strid, dshape)
    end if

    call h5dwrite_f(dsetid, strid, dvalue, shapeinfo, ierr)
    call h5dclose_f(dsetid, ierr)

    ok = (ierr == 0)
    
end function

! function h5simple_make_dset_int8(fileid, dname, dvalue, dshape, propertyid) result(ok)
! 
!     integer(HID_T), intent(in)      :: fileid
!     character(len=*), intent(in)    :: dname
!     integer, intent(in)             :: dvalue(*)
!     integer, intent(in)             :: dshape(:)
! 
!     integer(HID_T), intent(in), optional :: propertyid
! 
!     integer(HSIZE_T), allocatable   :: shapeinfo(:)
!     integer(HID_T)                  :: dsetid
!     integer                         :: ierr
! 
!     logical                         :: ok
! 
!     allocate(shapeinfo(SIZE(dshape)))
!     shapeinfo = dshape
! 
!     if (PRESENT(propertyid)) then
!         dsetid = h5simple_prep_dset(fileid, dname, H5T_NATIVE_INTEGER, dshape, propertyid=propertyid)
!     else                                    
!         dsetid = h5simple_prep_dset(fileid, dname, H5T_NATIVE_INTEGER, dshape)
!     end if
! 
!     call h5dwrite_f(dsetid, H5T_NATIVE_INTEGER, dvalue, shapeinfo, ierr)
!     call h5dclose_f(dsetid, ierr)
!     
!     ok = (ierr == 0)
! 
! end function

function h5simple_make_dset_integr(fileid, dname, dvalue, dshape, propertyid) result(ok)

    integer(HID_T), intent(in)      :: fileid
    character(len=*), intent(in)    :: dname
    integer, intent(in)             :: dvalue(*)
    integer, intent(in)             :: dshape(:)

    integer(HID_T), intent(in), optional :: propertyid

    integer(HSIZE_T), allocatable   :: shapeinfo(:)
    integer(HID_T)                  :: dsetid
    integer                         :: ierr

    logical                         :: ok

    allocate(shapeinfo(SIZE(dshape)))
    shapeinfo = dshape

    if (PRESENT(propertyid)) then
        dsetid = h5simple_prep_dset(fileid, dname, H5T_NATIVE_INTEGER, dshape, propertyid=propertyid)
    else                                    
        dsetid = h5simple_prep_dset(fileid, dname, H5T_NATIVE_INTEGER, dshape)
    end if

    call h5dwrite_f(dsetid, H5T_NATIVE_INTEGER, dvalue, shapeinfo, ierr)
    call h5dclose_f(dsetid, ierr)
    
    ok = (ierr == 0)

end function

function h5simple_make_dset_double(fileid, dname, dvalue, dshape, propertyid, doclose) result(ok)

    integer(HID_T), intent(in)      :: fileid
    character(len=*), intent(in)    :: dname
    real(dp), intent(in)                :: dvalue(*)
    integer, intent(in)             :: dshape(:)

    integer(HID_T), intent(in), optional    :: propertyid
    logical, intent(in), optional           :: doclose

    integer(HSIZE_T), allocatable   :: shapeinfo(:)
    integer(HID_T)                  :: dsetid
    integer                         :: ierr

    logical                         :: ok

    allocate(shapeinfo(SIZE(dshape)))
    shapeinfo = dshape

    if (PRESENT(propertyid)) then
        dsetid = h5simple_prep_dset(fileid, dname, H5T_NATIVE_DOUBLE, dshape, propertyid=propertyid)
    else                                    
        dsetid = h5simple_prep_dset(fileid, dname, H5T_NATIVE_DOUBLE, dshape)
    end if

    call h5dwrite_f(dsetid, H5T_NATIVE_DOUBLE, dvalue, shapeinfo, ierr)
    call h5dclose_f(dsetid, ierr)
    
    ok = (ierr == 0)

end function

function h5simple_make_attr_integr(fid, path, dvalue, dshape) result(ierr)

    integer(HID_T), intent(in)      :: fid
    character(len=*), intent(in)    :: path

    integer                         :: dvalue(*)
    integer, intent(in)             :: dshape(:)

    character(len=32)                :: dname
    character(len=32)                :: aname

    integer(HID_T)                  :: dsetid
    integer(HSIZE_T), allocatable   :: shapeinfo(:)
    integer(HID_T)                  :: attrid, spacid

    integer                         :: ierr

    dname = path(:index(path,"/@")-1)
    aname = path(index(path,"/@")+2:)

    CALL h5dopen_f(fid, dname, dsetid, ierr)

    allocate(shapeinfo(SIZE(dshape)))
    shapeinfo = dshape

    call h5screate_simple_f(size(shapeinfo), shapeinfo, spacid, ierr)

    call h5acreate_f(dsetid, trim(aname), H5T_NATIVE_INTEGER, spacid, attrid, ierr)
    call h5awrite_f(attrid, H5T_NATIVE_INTEGER, dvalue, shapeinfo, ierr)

    call h5aclose_f(attrid, ierr)
    call h5sclose_f(spacid, ierr)
    call h5dclose_f(dsetid, ierr)
    
end function


function h5simple_make_attr_string(fid, path, dvalue, dshape, strlen) result(ierr)

    integer(HID_T), intent(in)      :: fid
    character(len=*), intent(in)    :: path
    integer, intent(in)             :: strlen

    character(len=strlen)           :: dvalue(*)
    integer, intent(in)             :: dshape(:)

    character(len=32)                :: dname
    character(len=32)                :: aname

    integer(HID_T)                  :: dsetid
    integer(HSIZE_T), allocatable   :: shapeinfo(:)
    integer(HID_T)                  :: attrid, spacid, strid
    integer(SIZE_T)                 :: tstrlen

    integer                         :: ierr

    dname = path(:index(path,"/@")-1)
    aname = path(index(path,"/@")+2:)

    CALL h5dopen_f(fid, dname, dsetid, ierr)
    tstrlen = strlen

    call h5tcopy_f(H5T_NATIVE_CHARACTER, strid, ierr)
    call h5tset_size_f(strid, tstrlen, ierr)

    allocate(shapeinfo(SIZE(dshape)))
    shapeinfo = dshape

    call h5screate_simple_f(size(shapeinfo), shapeinfo, spacid, ierr)

    call h5acreate_f(dsetid, trim(aname), strid, spacid, attrid, ierr)
    call h5awrite_f(attrid, strid, dvalue, shapeinfo, ierr)

    call h5aclose_f(attrid, ierr)
    call h5sclose_f(spacid, ierr)
    call h5dclose_f(dsetid, ierr)
    
end function

function h5simple_read_attr_string(fid, path, dvalue, dshape) result(ierr)

    integer(HID_T), intent(in)      :: fid
    character(len=*), intent(in)    :: path

    character(len=*),intent(inout)  :: dvalue(*)
    integer, intent(in)             :: dshape(:)

    character(len=32)               :: dname
    character(len=32)               :: aname

    integer(HID_T)                  :: dsetid
    integer(HSIZE_T), allocatable   :: shapeinfo(:)
    integer(HID_T)                  :: attrid, strid
    integer(SIZE_T)                 :: tstrlen

    integer                         :: ierr

    dname = path(:index(path,"/@")-1)
    aname = path(index(path,"/@")+2:)

    tstrlen = len(dvalue(1))
    call h5tcopy_f(H5T_NATIVE_CHARACTER, strid, ierr)
    call h5tset_size_f(strid, tstrlen, ierr)

    allocate(shapeinfo(SIZE(dshape)))
    shapeinfo = dshape

    call h5dopen_f(fid, dname, dsetid, ierr)

    call h5aopen_f(dsetid, aname, attrid, ierr, H5P_DEFAULT_F)
    call h5aread_f(attrid, strid, dvalue, shapeinfo, ierr)

    call h5aclose_f(attrid, ierr)
    call h5dclose_f(dsetid, ierr)
    
end function

function h5simple_read_attr_integr(fid, path, dvalue, dshape) result(ierr)

    integer(HID_T), intent(in)      :: fid
    character(len=*), intent(in)    :: path

    integer ,intent(inout)          :: dvalue(*)
    integer, intent(in)             :: dshape(:)

    character(len=32)               :: dname
    character(len=32)               :: aname

    integer(HID_T)                  :: dsetid
    integer(HSIZE_T), allocatable   :: shapeinfo(:)
    integer(HID_T)                  :: attrid

    integer                         :: ierr

    dname = path(:index(path,"/@")-1)
    aname = path(index(path,"/@")+2:)

    allocate(shapeinfo(SIZE(dshape)))
    shapeinfo = dshape

    call h5dopen_f(fid, dname, dsetid, ierr)

    call h5aopen_f(dsetid, aname, attrid, ierr, H5P_DEFAULT_F)
    call h5aread_f(attrid, H5T_NATIVE_INTEGER, dvalue, shapeinfo, ierr)

    call h5aclose_f(attrid, ierr)
    call h5dclose_f(dsetid, ierr)
    
end function

function h5simple_read_attr_double(fid, path, dvalue, dshape) result(ierr)

    integer(HID_T), intent(in)      :: fid
    character(len=*), intent(in)    :: path

    real(dp),intent(inout)             :: dvalue(*)
    integer, intent(in)             :: dshape(:)

    character(len=32)               :: dname
    character(len=32)               :: aname

    integer(HID_T)                  :: dsetid
    integer(HSIZE_T), allocatable   :: shapeinfo(:)
    integer(HID_T)                  :: attrid

    integer                         :: ierr

    dname = path(:index(path,"/@")-1)
    aname = path(index(path,"/@")+2:)

    allocate(shapeinfo(SIZE(dshape)))
    shapeinfo = dshape

    call h5dopen_f(fid, dname, dsetid, ierr)

    call h5aopen_f(dsetid, aname, attrid, ierr, H5P_DEFAULT_F)
    call h5aread_f(attrid, H5T_NATIVE_DOUBLE, dvalue, shapeinfo, ierr)

    call h5aclose_f(attrid, ierr)
    call h5dclose_f(dsetid, ierr)
    
end function

function h5simple_fill_dset_byte(fileid, dname, dvalue, dshape) result(ierr)

    integer(HID_T), intent(in)      :: fileid
    character(len=*), intent(in)    :: dname
    character(len=1), intent(in)    :: dvalue(*)
    integer, intent(in)             :: dshape(:)

    integer(HSIZE_T), allocatable   :: shapeinfo(:)
    integer(HID_T)                  :: dsetid
    integer                         :: ierr

    allocate(shapeinfo(SIZE(dshape)))
    shapeinfo = dshape

    CALL h5dopen_f(fileid, dname, dsetid, ierr)
    !CALL h5dwrite_f(dsetid, H5T_NATIVE_B8, dvalue, shapeinfo, ierr)
    CALL h5dwrite_f(dsetid, H5T_NATIVE_CHARACTER, dvalue, shapeinfo, ierr)
    call h5dclose_f(dsetid, ierr)

end function


function h5simple_fill_dset_char(fileid, dname, dvalue, dshape) result(ierr)

    integer(HID_T), intent(in)      :: fileid
    character(len=*), intent(in)    :: dname
    character(len=1), intent(in)    :: dvalue(*)
    integer, intent(in)             :: dshape(:)

    integer(HSIZE_T), allocatable   :: shapeinfo(:)
    integer(HID_T)                  :: dsetid
    integer                         :: ierr

    allocate(shapeinfo(SIZE(dshape)))
    shapeinfo = dshape

    CALL h5dopen_f(fileid, dname, dsetid, ierr)
    CALL h5dwrite_f(dsetid, H5T_NATIVE_CHARACTER, dvalue, shapeinfo, ierr)
    call h5dclose_f(dsetid, ierr)

end function

function h5simple_fill_dset_double(fileid, dname, dvalue, dshape) result(ierr)

    integer(HID_T), intent(in)      :: fileid
    character(len=*), intent(in)    :: dname
    real(dp), intent(in)                :: dvalue(*)
    integer, intent(in)             :: dshape(:)

    integer(HSIZE_T), allocatable   :: shapeinfo(:)
    integer(HID_T)                  :: dsetid
    integer                         :: ierr

    allocate(shapeinfo(SIZE(dshape)))
    shapeinfo = dshape

    CALL h5dopen_f(fileid, dname, dsetid, ierr)
    CALL h5dwrite_f(dsetid, H5T_NATIVE_DOUBLE, dvalue, shapeinfo, ierr)
    call h5dclose_f(dsetid, ierr)

end function

function h5simple_fill_dset_integr(fileid, dname, dvalue, dshape) result(ierr)

    integer(HID_T), intent(in)      :: fileid
    character(len=*), intent(in)    :: dname
    integer, intent(in)             :: dvalue(*)
    integer, intent(in)             :: dshape(:)

    integer(HSIZE_T), allocatable   :: shapeinfo(:)
    integer(HID_T)                  :: dsetid
    integer                         :: ierr

    allocate(shapeinfo(SIZE(dshape)))
    shapeinfo = dshape

    CALL h5dopen_f(fileid, dname, dsetid, ierr)
    CALL h5dwrite_f(dsetid, H5T_NATIVE_INTEGER, dvalue, shapeinfo, ierr)
    call h5dclose_f(dsetid, ierr)

end function

function h5simple_fill_dset_integr8(fileid, dname, dvalue, dshape) result(ierr)

    use iso_c_binding !! pmake: ignore
    integer(HID_T), intent(in)      :: fileid
    character(len=*), intent(in)    :: dname
    integer(c_int8_t), intent(in)   :: dvalue(*)
    integer, intent(in)             :: dshape(:)

    integer(HSIZE_T), allocatable   :: shapeinfo(:)
    integer(HID_T)                  :: dsetid
    integer                         :: ierr

    allocate(shapeinfo(SIZE(dshape)))
    shapeinfo = dshape

    CALL h5dopen_f(fileid, dname, dsetid, ierr)
    CALL h5dwrite_f(dsetid, H5T_STD_I8LE, dvalue, shapeinfo, ierr)
    call h5dclose_f(dsetid, ierr)

end function

function h5simple_fill_dset_integr32(fileid, dname, dvalue, dshape) result(ierr)

    use iso_c_binding !! pmake: ignore
    integer(HID_T), intent(in)      :: fileid
    character(len=*), intent(in)    :: dname
    integer(c_int32_t), intent(in)  :: dvalue(*)
    integer, intent(in)             :: dshape(:)

    integer(HSIZE_T), allocatable   :: shapeinfo(:)
    integer(HID_T)                  :: dsetid
    integer                         :: ierr

    allocate(shapeinfo(SIZE(dshape)))
    shapeinfo = dshape

    CALL h5dopen_f(fileid, dname, dsetid, ierr)
    CALL h5dwrite_f(dsetid, H5T_STD_I32LE, dvalue, shapeinfo, ierr)
    call h5dclose_f(dsetid, ierr)
end function

function h5simple_read_dset_string(fileid, dname, dvalue, dshape) result(ierr)

    integer(HID_T), intent(in)      :: fileid
    character(len=*), intent(in)    :: dname
    character(len=*), intent(out)   :: dvalue(*)
    integer, intent(in)             :: dshape(:)

    integer(HSIZE_T), allocatable   :: shapeinfo(:)
    integer(HID_T)                  :: dsetid, strid
    integer                         :: ierr

    integer(SIZE_T)                 :: tstrlen

    tstrlen = len(dvalue(1))

    call h5tcopy_f(H5T_NATIVE_CHARACTER, strid, ierr)
    call h5tset_size_f(strid, tstrlen, ierr)

    allocate(shapeinfo(SIZE(dshape)))
    shapeinfo = dshape

    CALL h5dopen_f(fileid, dname, dsetid, ierr)
    CALL h5dread_f(dsetid, strid, dvalue, shapeinfo, ierr)
    call h5dclose_f(dsetid, ierr)

end function

function h5simple_read_dset_double(fileid, dname, dvalue, dshape) result(ierr)

    integer(HID_T), intent(in)      :: fileid
    character(len=*), intent(in)    :: dname
    real(dp), intent(out)           :: dvalue(*)
    integer, intent(in)             :: dshape(:)

    integer(HSIZE_T), allocatable   :: shapeinfo(:)
    integer(HID_T)                  :: dsetid
    integer                         :: ierr

    allocate(shapeinfo(SIZE(dshape)))
    shapeinfo = dshape

    CALL h5dopen_f(fileid, dname, dsetid, ierr)
    CALL h5dread_f(dsetid, H5T_NATIVE_DOUBLE, dvalue, shapeinfo, ierr)
    call h5dclose_f(dsetid, ierr)

end function

function h5simple_read_dset_integr(fileid, dname, dvalue, dshape) result(ierr)

    integer(HID_T), intent(in)      :: fileid
    character(len=*), intent(in)    :: dname
    integer, intent(out)            :: dvalue(*)
    integer, intent(in)             :: dshape(:)

    integer(HSIZE_T), allocatable   :: shapeinfo(:)
    integer(HID_T)                  :: dsetid
    integer                         :: ierr

    allocate(shapeinfo(SIZE(dshape)))
    shapeinfo = dshape

    CALL h5dopen_f(fileid, dname, dsetid, ierr)
    CALL h5dread_f(dsetid, H5T_NATIVE_INTEGER, dvalue, shapeinfo, ierr)
    call h5dclose_f(dsetid, ierr)

end function

function h5simple_dtype_isequal(a,b) result(isequal)

    integer(HID_T), intent(in) :: a,b

    integer :: errno
    logical :: isequal

    call H5Tequal_f(a,b, isequal, errno)

end function

function h5simple_dset_get_dtype(fileid,dpath,dtype) result(ok)

    integer(HID_T), intent(in)      :: fileid
    character(len=*), intent(in)    :: dpath
    integer(HID_T), intent(out)     :: dtype
    logical                         :: ok

    integer :: errno

    integer(HID_T) :: dsetid
    integer(HID_T) :: dtypeid

    call h5dopen_f(fileid, trim(dpath), dsetid, errno)
    call H5Dget_type_f(dsetid, dtypeid, errno)

    if (h5simple_dtype_isequal(H5T_STD_I8LE, dtypeid)) then
        dtype = h5simple_dtype_int8
    else if (h5simple_dtype_isequal(H5T_STD_I32LE, dtypeid)) then
        dtype = h5simple_dtype_int32
    else if (h5simple_dtype_isequal(H5T_NATIVE_REAL, dtypeid)) then
        dtype = h5simple_dtype_float32
    else if (h5simple_dtype_isequal(H5T_NATIVE_DOUBLE, dtypeid)) then
        dtype = h5simple_dtype_float64
    end if

    call H5Tclose_f(dtypeid, errno)
    call H5Dclose_f(dsetid, errno)

    ok = (errno == 0)

end function

function h5simple_dset_get_shape(fileid,dpath,drank,dshape) result(ok)

    integer(HID_T), intent(in)      :: fileid
    character(len=*), intent(in)    :: dpath
    integer(kind=4), intent(out)    :: drank,dshape(16)
    logical                         :: ok

    integer :: errno

    integer(HID_T) :: dsetid
    integer(HID_T) :: dspaceid

    integer(HSIZE_T) :: dims(16),maxdims(16)

    call h5dopen_f(fileid, trim(dpath), dsetid, errno)
    call h5dget_space_f(dsetid, dspaceid, errno)

    call h5sget_simple_extent_ndims_f(dspaceid, drank, errno)
    call h5sget_simple_extent_dims_f(dspaceid,dims,maxdims,errno)

    call H5Sclose_f(dspaceid, errno)
    call H5Dclose_f(dsetid, errno)

    dshape = -1
    dshape(1:drank) = INT(dims(1:drank),kind=4)

    ok = (errno == 0)

end function

end module
