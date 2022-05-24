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
# define p4wrap_cell_t p8wrap_cell_t
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

module p4est_constants_mod

use iso_c_binding, only: c_int8_t

integer(c_int8_t), parameter :: P4EST_MAXLEVEL = 30
integer(c_int8_t), parameter :: P8EST_MAXLEVEL = P4EST_MAXLEVEL

integer(c_int8_t), parameter :: P4EST_CONNECT_FACE   = 21
integer(c_int8_t), parameter :: P4EST_CONNECT_CORNER = 22
integer(c_int8_t), parameter :: P4EST_CONNECT_FULL   = P4EST_CONNECT_CORNER

integer(c_int8_t), parameter :: P8EST_CONNECT_FACE   = 31
integer(c_int8_t), parameter :: P8EST_CONNECT_EDGE   = 32
integer(c_int8_t), parameter :: P8EST_CONNECT_CORNER = 33
integer(c_int8_t), parameter :: P8EST_CONNECT_FULL   = P8EST_CONNECT_CORNER

# if defined(P4_TO_P8)
integer, parameter :: N_DIMS = 3
# else
integer, parameter :: N_DIMS = 2
# endif
integer, parameter :: P4EST_HALF = 2**(N_DIMS-1)
integer, parameter :: P4EST_FACES = 2*N_DIMS
integer, parameter :: P4EST_CORNERS = 2**N_DIMS
integer, parameter :: P4EST_CHILDREN = 2**N_DIMS

end module
