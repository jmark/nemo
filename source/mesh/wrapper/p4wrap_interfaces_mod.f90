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

module p4wrap_interfaces_mod

interface
# if defined(P4_TO_P8)
    subroutine p8wrap_gather_connection_info(p4est,layer,quads,sides,verts) bind(c,name="p8wrap_gather_connection_info")
# else
    subroutine p4wrap_gather_connection_info(p4est,layer,quads,sides,verts) bind(c,name="p4wrap_gather_connection_info")
# endif
        use iso_c_binding, only: c_ptr

        type(c_ptr), value, intent(in) :: p4est,layer
        type(c_ptr), value, intent(in) :: quads,sides,verts
    end subroutine
end interface

interface
# if defined(P4_TO_P8)
    subroutine p8wrap_gather_mirror_proc_mirrors(layer,len,proc_mirrors) bind(c,name="p8wrap_gather_mirror_proc_mirrors")
# else
    subroutine p4wrap_gather_mirror_proc_mirrors(layer,len,proc_mirrors) bind(c,name="p4wrap_gather_mirror_proc_mirrors")
# endif
        use iso_c_binding, only: c_ptr, c_int

        type(c_ptr), value, intent(in) :: layer,proc_mirrors
        integer(c_int), value, intent(in) :: len
    end subroutine
end interface

interface
# if defined(P4_TO_P8)
    subroutine p8wrap_gather_mirror_proc_offsets(layer,len,proc_offsets) bind(c,name="p8wrap_gather_mirror_proc_offsets")
# else
    subroutine p4wrap_gather_mirror_proc_offsets(layer,len,proc_offsets) bind(c,name="p4wrap_gather_mirror_proc_offsets")
# endif
        use iso_c_binding, only: c_ptr, c_int

        type(c_ptr), value, intent(in) :: layer,proc_offsets
        integer(c_int), value, intent(in) :: len
    end subroutine
end interface

interface
# if defined(P4_TO_P8)
    subroutine p8wrap_gather_proc_offsets(layer,len,proc_offsets) bind(c,name="p8wrap_gather_proc_offsets")
# else
    subroutine p4wrap_gather_proc_offsets(layer,len,proc_offsets) bind(c,name="p4wrap_gather_proc_offsets")
# endif
        use iso_c_binding, only: c_ptr, c_int

        type(c_ptr), value, intent(in) :: layer,proc_offsets
        integer(c_int), value, intent(in) :: len
    end subroutine
end interface

interface
# if defined(P4_TO_P8)
    subroutine p8wrap_gather_mirroridx(layer,mirroridx) bind(c,name="p8wrap_gather_mirroridx")
# else
    subroutine p4wrap_gather_mirroridx(layer,mirroridx) bind(c,name="p4wrap_gather_mirroridx")
# endif
        use iso_c_binding, only: c_ptr

        type(c_ptr), value, intent(in) :: layer,mirroridx
    end subroutine
end interface

interface
# if defined(P4_TO_P8)
    subroutine p8wrap_gather_mirrorptr(layer,mirrorptr,esize) bind(c,name="p8wrap_gather_mirrorptr")
# else
    subroutine p4wrap_gather_mirrorptr(layer,mirrorptr,esize) bind(c,name="p4wrap_gather_mirrorptr")
# endif
        use iso_c_binding, only: c_ptr, c_size_t

        type(c_ptr), value, intent(in) :: layer,mirrorptr
        integer(c_size_t), value, intent(in) :: esize
    end subroutine
end interface

interface
# if defined(P4_TO_P8)
    subroutine p8wrap_refine(p4est,flags) bind(c,name="p8wrap_refine")
# else
    subroutine p4wrap_refine(p4est,flags) bind(c,name="p4wrap_refine")
# endif
        use iso_c_binding, only: c_ptr

        type(c_ptr), value, intent(in) :: p4est,flags
    
    end subroutine
end interface

interface
# if defined(P4_TO_P8)
    subroutine p8wrap_coarsen(p4est,flags) bind(c,name="p8wrap_coarsen")
# else
    subroutine p4wrap_coarsen(p4est,flags) bind(c,name="p4wrap_coarsen")
# endif
        use iso_c_binding, only: c_ptr

        type(c_ptr), value, intent(in) :: p4est,flags
    end subroutine
end interface

interface
# if defined(P4_TO_P8)
    subroutine p8wrap_balance(p4est) bind(c,name="p8wrap_balance")
# else
    subroutine p4wrap_balance(p4est) bind(c,name="p4wrap_balance")
# endif
        use iso_c_binding, only: c_ptr

        type(c_ptr), value, intent(in) :: p4est
    end subroutine
end interface

interface
# if defined(P4_TO_P8)
    subroutine p8wrap_partition(p4est) bind(c,name="p8wrap_partition")
# else
    subroutine p4wrap_partition(p4est) bind(c,name="p4wrap_partition")
# endif
        use iso_c_binding, only: c_ptr

        type(c_ptr), value, intent(in) :: p4est
    end subroutine
end interface

interface
# if defined(P4_TO_P8)
    subroutine p8wrap_copy_gfq(p4est,count,gfq) bind(c,name="p8wrap_copy_gfq")
# else
    subroutine p4wrap_copy_gfq(p4est,count,gfq) bind(c,name="p4wrap_copy_gfq")
# endif
        use iso_c_binding, only: c_ptr, c_int

        type(c_ptr), value, intent(in)      :: p4est
        integer(c_int), value, intent(in)   :: count
        type(c_ptr), value, intent(in)      :: gfq

    end subroutine
end interface

interface
# if defined(P4_TO_P8)
    pure function p8wrap_get_num_ghosts(layer) bind(c,name="p8wrap_get_num_ghosts") result(num)
# else
    pure function p4wrap_get_num_ghosts(layer) bind(c,name="p4wrap_get_num_ghosts") result(num)
# endif
        use iso_c_binding, only: c_size_t, c_ptr
        type(c_ptr), value, intent(in)  :: layer
        integer(c_size_t)               :: num
    end function
end interface

interface
# if defined(P4_TO_P8)
    pure function p8wrap_get_num_mirrors(layer) bind(c,name="p8wrap_get_num_mirrors") result(num)
# else
    pure function p4wrap_get_num_mirrors(layer) bind(c,name="p4wrap_get_num_mirrors") result(num)
# endif
        use iso_c_binding, only: c_size_t, c_ptr
        type(c_ptr), value, intent(in)  :: layer
        integer(c_size_t)               :: num
    end function
end interface

interface
# if defined(P4_TO_P8)
    pure function p8wrap_get_num_quads(p4est) bind(c,name="p8wrap_get_num_quads") result(num)
# else
    pure function p4wrap_get_num_quads(p4est) bind(c,name="p4wrap_get_num_quads") result(num)
# endif
        use iso_c_binding, only: c_size_t, c_ptr
        type(c_ptr), value, intent(in)  :: p4est
        integer(c_size_t)               :: num
    end function
end interface

interface
# if defined(P4_TO_P8)
    pure function p8wrap_get_num_global(p4est) bind(c,name="p8wrap_get_num_global") result(num)
# else
    pure function p4wrap_get_num_global(p4est) bind(c,name="p4wrap_get_num_global") result(num)
# endif
        use iso_c_binding, only: c_size_t, c_ptr
        type(c_ptr), value, intent(in)  :: p4est
        integer(c_size_t)               :: num
    end function
end interface

interface
# if defined(P4_TO_P8)
    pure function p8wrap_get_idx_global(p4est) bind(c,name="p8wrap_get_idx_global") result(idx)
# else
    pure function p4wrap_get_idx_global(p4est) bind(c,name="p4wrap_get_idx_global") result(idx)
# endif
        use iso_c_binding, only: c_size_t, c_ptr
        type(c_ptr), value, intent(in)  :: p4est
        integer(c_size_t)               :: idx
    end function
end interface

interface
# if defined(P4_TO_P8)
    pure function p8wrap_get_num_sides(p4est,layer) bind(c,name="p8wrap_get_num_sides") result(num)
# else
    pure function p4wrap_get_num_sides(p4est,layer) bind(c,name="p4wrap_get_num_sides") result(num)
# endif
        use iso_c_binding, only: c_size_t, c_ptr
        type(c_ptr), value, intent(in)  :: p4est,layer
        integer(c_size_t)               :: num
    end function
end interface

interface
# if defined(P4_TO_P8)
    pure function p8wrap_get_num_verts(p4est,layer) bind(c,name="p8wrap_get_num_verts") result(num)
# else
    pure function p4wrap_get_num_verts(p4est,layer) bind(c,name="p4wrap_get_num_verts") result(num)
# endif
        use iso_c_binding, only: c_size_t, c_ptr
        type(c_ptr), value, intent(in)  :: p4est,layer
        integer(c_size_t)               :: num
    end function
end interface

interface
# if defined(P4_TO_P8)
    pure function p8wrap_get_num_verts_hanging(p4est,layer) bind(c,name="p8wrap_get_num_verts_hanging") result(num)
# else
    pure function p4wrap_get_num_verts_hanging(p4est,layer) bind(c,name="p4wrap_get_num_verts_hanging") result(num)
# endif
        use iso_c_binding, only: c_size_t, c_ptr
        type(c_ptr), value, intent(in)  :: p4est,layer
        integer(c_size_t)               :: num
    end function
end interface

interface
# if defined(P4_TO_P8)
    pure function p8wrap_quadrant_child_id(level,x,y,z) bind(c,name="p8wrap_quadrant_child_id") result(id)
# else
    pure function p4wrap_quadrant_child_id(level,x,y) bind(c,name="p4wrap_quadrant_child_id") result(id)
# endif
        use iso_c_binding, only: c_int, c_int8_t, c_int32_t
        integer(c_int8_t), value, intent(in) :: level
        integer(c_int32_t), value, intent(in) :: x,y
# if defined(P4_TO_P8)
        integer(c_int32_t), value, intent(in) :: z
# endif
        integer(c_int)                  :: id
    end function
end interface

end module
