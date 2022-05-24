# if PP_N_DIMS == 3

# define P4_TO_P8 1

# define p4est p8est
# define p4est_quadrant_t p8est_quadrant_t
# define p4est_quadrant_data_t p8est_quadrant_data_t

# define P4EST_HALF P8EST_HALF
# define P4EST_FACES P8EST_FACES
# define P4EST_CORNERS P8EST_CORNERS
# define P4EST_CHILDREN P8EST_CHILDREN

# define p4est_types_mod p8est_types_mod
# define p4est_constants_mod p8est_constants_mod
# define p4est_interfaces_mod p8est_interfaces_mod

# define p4wrap_types_mod p8wrap_types_mod
# define p4wrap_constants_mod p8wrap_constants_mod
# define p4wrap_interfaces_mod p8wrap_interfaces_mod

# define p4wrap_mesh_t p8wrap_mesh_t
# define p4wrap_quad_t p8wrap_quad_t
# define p4wrap_side_t p8wrap_side_t

# define P4WRAP_FLAG_NONE P8WRAP_FLAG_NONE
# define P4WRAP_FLAG_REFINE P8WRAP_FLAG_REFINE
# define P4WRAP_FLAG_COARSEN P8WRAP_FLAG_COARSEN
# define P4WRAP_FLAG_ISGHOST P8WRAP_FLAG_ISGHOST

# define P4WRAP_FLAG_NONE P8WRAP_FLAG_NONE
# define P4WRAP_FLAG_REFINE P8WRAP_FLAG_REFINE
# define P4WRAP_FLAG_COARSEN P8WRAP_FLAG_COARSEN
# define P4WRAP_FLAG_ISGHOST P8WRAP_FLAG_ISGHOST

# define P4WRAP_BOUNDARY_NULL P8WRAP_BOUNDARY_NULL
# define P4WRAP_BOUNDARY_NOR P8WRAP_BOUNDARY_NOR
# define P4WRAP_BOUNDARY_SOU P8WRAP_BOUNDARY_SOU
# define P4WRAP_BOUNDARY_WES P8WRAP_BOUNDARY_WES
# define P4WRAP_BOUNDARY_EAS P8WRAP_BOUNDARY_EAS

# endif

module p4wrap_types_mod

use iso_c_binding, only: c_ptr, c_int, c_int8_t, c_int16_t, c_int32_t, c_size_t

use p4est_constants_mod, only: N_DIMS
use p4est_constants_mod, only: P4EST_HALF
use p4est_constants_mod, only: P4EST_FACES
use p4est_constants_mod, only: P4EST_CORNERS
use p4est_constants_mod, only: P4EST_CHILDREN

type, bind(c) :: p4wrap_quad_t

    !! Do not change order !!
    integer(c_int32_t)  :: morton(N_DIMS)
    integer(c_int8_t)   :: level

    integer(c_size_t)   :: cellid
    integer(c_int32_t)  :: treeid
    integer(c_int32_t)  :: quadid

    integer(c_size_t)   :: isides(P4EST_FACES)
    integer(c_size_t)   :: iverts(P4EST_CORNERS)

    integer(c_int8_t)   :: flags

    integer(c_int32_t)  :: local_num

end type

type, bind(c) :: p4wrap_side_t

    integer(c_int8_t) :: tree_boundary
    integer(c_int8_t) :: faceid(2)
    integer(c_int8_t) :: nquads(2)
    integer(c_size_t) :: iquads(P4EST_HALF,2)

end type

type, bind(c) :: p4wrap_vert_t

    integer(c_int8_t) :: tree_boundary
    integer(c_int8_t) :: nquads = 0
    integer(c_size_t) :: iquads(P4EST_CORNERS) = -1

end type

end module
