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

# define p4wrap_gather_mirroridx p8wrap_gather_mirroridx
# define p4wrap_gather_mirrorptr p8wrap_gather_mirrorptr
# define p4wrap_gather_connection_info p8wrap_gather_connection_info

# define p4wrap_get_num_ghosts p8wrap_get_num_ghosts
# define p4wrap_get_num_mirrors p8wrap_get_num_mirrors
# define p4wrap_get_num_quads p8wrap_get_num_quads
# define p4wrap_get_num_global p8wrap_get_num_global
# define p4wrap_get_idx_global p8wrap_get_idx_global
# define p4wrap_get_num_sides p8wrap_get_num_sides

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

module p4wrap_constants_mod

use iso_c_binding, only: c_int8_t

integer(c_int8_t), parameter :: P4WRAP_FLAG_NONE      = 0
integer(c_int8_t), parameter :: P4WRAP_FLAG_REFINE    = 1
integer(c_int8_t), parameter :: P4WRAP_FLAG_COARSEN   = 2
integer(c_int8_t), parameter :: P4WRAP_FLAG_ISGHOST   = 4

!! integer(c_int8_t), parameter :: P4WRAP_SIDE_NONE   = 0
!! integer(c_int8_t), parameter :: P4WRAP_SIDE_1v1    = 1
!! integer(c_int8_t), parameter :: P4WRAP_SIDE_1vN    = 2
!! integer(c_int8_t), parameter :: P4WRAP_SIDE_Nv1    = 3

integer, parameter :: P4WRAP_BOUNDARY_NULL = 0
integer, parameter :: P4WRAP_BOUNDARY_NOR  = 1
integer, parameter :: P4WRAP_BOUNDARY_SOU  = 2
integer, parameter :: P4WRAP_BOUNDARY_WES  = 4
integer, parameter :: P4WRAP_BOUNDARY_EAS  = 8

# if defined(P4_TO_P8)
integer, parameter :: P8WRAP_BOUNDARY_FRO  = 16
integer, parameter :: P8WRAP_BOUNDARY_BAC  = 32
# endif

end module
