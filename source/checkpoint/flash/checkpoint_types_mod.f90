module checkpoint_types_mod

use constants_mod
use hdf5, only: HID_T, HSIZE_T

character(len=*), parameter :: checkpoint_directory = 'checkpoints'
character(len=*), parameter :: checkpoint_prefix = 'chkpt_' 

integer, parameter              :: maxentries = 256

type :: checkpoint_t

    integer                         :: maxentries = maxentries

    integer                         :: iocount

    character(len=4096)             :: partial_suffix = '.partial'
    character(len=4096)             :: prefix = checkpoint_prefix
    character(len=4096)             :: directory = checkpoint_directory

    character(len=4096)             :: filename = ''

    character(len=4096)             :: extension = "f5"

    character(len=80)    :: keys_real_runtime_parameters(maxentries)
    character(len=80)    :: keys_string_runtime_parameters(maxentries)
    character(len=80)    :: keys_integer_runtime_parameters(maxentries)
    character(len=80)    :: keys_logical_runtime_parameters(maxentries)

    real(dp)             :: vals_real_runtime_parameters(maxentries)
    character(len=80)    :: vals_string_runtime_parameters(maxentries)
    integer              :: vals_integer_runtime_parameters(maxentries)
    logical              :: vals_logical_runtime_parameters(maxentries)

    integer :: index_real_runtime_parameters = 0
    integer :: index_string_runtime_parameters = 0
    integer :: index_integer_runtime_parameters = 0
    integer :: index_logical_runtime_parameters = 0

    character(len=80)    :: keys_real_scalars(maxentries)
    character(len=80)    :: keys_string_scalars(maxentries)
    character(len=80)    :: keys_integer_scalars(maxentries)
    character(len=80)    :: keys_logical_scalars(maxentries)

    real(dp)             :: vals_real_scalars(maxentries)
    character(len=80)    :: vals_string_scalars(maxentries)
    integer              :: vals_integer_scalars(maxentries)
    logical              :: vals_logical_scalars(maxentries)

    integer :: index_real_scalars = 0
    integer :: index_string_scalars = 0
    integer :: index_integer_scalars = 0
    integer :: index_logical_scalars = 0

    integer             :: index_unknown_names = 0
    character(len=4)    :: unknown_names(maxentries)

    integer(HID_T)      :: fileid
    integer(HID_T)      :: propertyid

    integer :: file_format_version = 9

    character(len=400) :: setup_call
    character(len= 80) :: file_creation_time
    character(len= 80) :: flash_version
    character(len= 80) :: build_date
    character(len= 80) :: build_dir
    character(len= 80) :: build_machine
    character(len=400) :: cflags
    character(len=400) :: fflags
    character(len= 80) :: setup_time_stamp
    character(len= 80) :: build_time_stamp

end type

type :: checkpoint_dataset_t

    type(checkpoint_t), pointer     :: chkpt
    integer(HID_T)                  :: dsetid

    character(len=256)              :: path

    integer(HSIZE_T)                :: tsize !! datatype size in *bytes*
    integer(HSIZE_T)                :: bsize !! buffer size in *bytes*
    integer(HSIZE_T)                :: dsize !! dataset size *per quad* in *bytes*

    character(len=1), allocatable   :: buffer(:)

end type

interface
    subroutine checkpoint_iterate_cb_t(dset,quad)

        use mesh_types_mod, only: quad_t
        import checkpoint_dataset_t

        type(checkpoint_dataset_t), intent(inout)   :: dset
        type(quad_t), intent(inout)                 :: quad

    end subroutine
end interface

end module
