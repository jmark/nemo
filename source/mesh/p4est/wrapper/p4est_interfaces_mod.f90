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

module p4est_interfaces_mod

interface
# if defined(P4_TO_P8)
    function p8est_new_ext(&
        mpicomm, &          !! sc_MPI_Comm mpicomm
        connectivity, &     !! p8est_connectivity_t *connectivity
        min_quads, &        !! p8est_locidx_t min_quadrants
        min_level, &        !! int min_level
        fill_uniform, &     !! int fill_uniform
        data_size, &        !! size_t data_size
        init_fn, &          !! p8est_init_t init_fn
        userptr &           !! void *user_pointer
    ) result(p8est) bind(c,name="p8est_new_ext") !! p8est_t*
# else
    function p4est_new_ext( &
        mpicomm, &          !! sc_MPI_Comm mpicomm
        connectivity, &     !! p4est_connectivity_t *connectivity
        min_quads, &        !! p4est_locidx_t min_quadrants
        min_level, &        !! int min_level
        fill_uniform, &     !! int fill_uniform
        data_size, &        !! size_t data_size
        init_fn, &          !! p4est_init_t init_fn
        userptr &           !! void *user_pointer
    ) result(p4est) bind(c,name="p4est_new_ext") !! p4est_t*
# endif
        use iso_c_binding, only: c_ptr, c_int, c_int32_t, c_size_t

# if defined(__INTEL_COMPILER)
        integer(c_int),     intent(in), value :: mpicomm
# else
        type(c_ptr),        intent(in), value :: mpicomm
# endif
        type(c_ptr),        intent(in), value :: connectivity
        integer(c_int32_t), intent(in), value :: min_quads
        integer(c_int),     intent(in), value :: min_level
        integer(c_int),     intent(in), value :: fill_uniform
        integer(c_size_t),  intent(in), value :: data_size
        type(c_ptr),        intent(in), value :: init_fn
        type(c_ptr),        intent(in), value :: userptr

        type(c_ptr) :: p4est
    end function
end interface

interface
# if defined(P4_TO_P8)
    function p8est_ghost_new(p8est,btype) result(ghost) bind(c,name="p8est_ghost_new")
# else
    function p4est_ghost_new(p4est,btype) result(ghost) bind(c,name="p4est_ghost_new")
# endif
        use iso_c_binding, only: c_ptr, c_int8_t

        type(c_ptr),        intent(in), value   :: p4est
        integer(c_int8_t),  intent(in), value   :: btype
        type(c_ptr)                             :: ghost
    end function
end interface

interface
# if defined(P4_TO_P8)
    subroutine p8est_ghost_destroy(ghost) bind(c,name="p8est_ghost_destroy")
# else
    subroutine p4est_ghost_destroy(ghost) bind(c,name="p4est_ghost_destroy")
# endif
        use iso_c_binding, only: c_ptr

        type(c_ptr), intent(in), value :: ghost
    end subroutine
end interface

interface
# if defined(P4_TO_P8)
    subroutine p8est_destroy(p8est) bind(c,name="p8est_destroy")
# else
    subroutine p4est_destroy(p4est) bind(c,name="p4est_destroy")
# endif
        use iso_c_binding, only: c_ptr

        type(c_ptr), intent(in), value :: p4est
    end subroutine
end interface

interface
# if defined(P4_TO_P8)
    subroutine p8est_connectivity_destroy(connectivity) bind(c,name="p8est_connectivity_destroy")
# else
    subroutine p4est_connectivity_destroy(connectivity) bind(c,name="p4est_connectivity_destroy")
# endif
        use iso_c_binding, only: c_ptr

        type(c_ptr), intent(in), value :: connectivity
    end subroutine
end interface

interface
# if defined(P4_TO_P8)
    subroutine p8est_init(log_handler, log_threshold) bind(c,name="p4est_init")
# else
    subroutine p4est_init(log_handler, log_threshold) bind(c,name="p4est_init")
# endif
        use iso_c_binding, only: c_ptr, c_int

        type(c_ptr), intent(in), value      :: log_handler
        integer(c_int), intent(in), value   :: log_threshold
    end subroutine
end interface

interface
# if defined(P4_TO_P8)
    function p8est_connectivity_new_periodic() result(connectivity) bind(c,name="p8est_connectivity_new_periodic")
# else
    function p4est_connectivity_new_periodic() result(connectivity) bind(c,name="p4est_connectivity_new_periodic")
# endif
        use iso_c_binding, only: c_ptr

        type(c_ptr) :: connectivity
    end function
end interface

interface
# if defined(P4_TO_P8)
    function p8est_connectivity_new_unitcube() result(connectivity) bind(c,name="p8est_connectivity_new_unitcube")
# else
    function p4est_connectivity_new_unitsquare() result(connectivity) bind(c,name="p4est_connectivity_new_unitsquare")
# endif
        use iso_c_binding, only: c_ptr

        type(c_ptr) :: connectivity
    end function
end interface

interface
# if defined(P4_TO_P8)
    subroutine p8est_ghost_exchange_data(p8est, ghost_layer, ghost_data) bind(c,name="p8est_ghost_exchange_data")
# else
    subroutine p4est_ghost_exchange_data(p4est, ghost_layer, ghost_data) bind(c,name="p4est_ghost_exchange_data")
# endif
        use iso_c_binding, only: c_ptr

        type(c_ptr), value, intent(in) :: p4est
        type(c_ptr), value, intent(in) :: ghost_layer
        type(c_ptr), value, intent(in) :: ghost_data
    end subroutine
end interface

interface
# if defined(P4_TO_P8)
    subroutine p8est_ghost_exchange_contiguous(p8est, ghost_layer, data_size, mirror_data, ghost_data) bind(c,name="p8est_ghost_exchange_contiguous")
# else
    subroutine p4est_ghost_exchange_contiguous(p4est, ghost_layer, data_size, mirror_data, ghost_data) bind(c,name="p4est_ghost_exchange_contiguous")
# endif
        use iso_c_binding, only: c_ptr, c_size_t

        type(c_ptr),        value,  intent(in)  :: p4est
        type(c_ptr),        value,  intent(in)  :: ghost_layer
        integer(c_size_t),  value,  intent(in)  :: data_size
        type(c_ptr),        value,  intent(in)  :: mirror_data
        type(c_ptr),        value,  intent(in)  :: ghost_data
    end subroutine
end interface

interface
# if defined(P4_TO_P8)
    subroutine p8est_ghost_exchange_custom(p8est, ghost_layer, data_size, mirror_data_ptrs, ghost_data) bind(c,name="p8est_ghost_exchange_custom")
# else
    subroutine p4est_ghost_exchange_custom(p4est, ghost_layer, data_size, mirror_data_ptrs, ghost_data) bind(c,name="p4est_ghost_exchange_custom")
# endif
        use iso_c_binding, only: c_ptr, c_size_t

        type(c_ptr),        value,  intent(in)  :: p4est
        type(c_ptr),        value,  intent(in)  :: ghost_layer
        integer(c_size_t),  value,  intent(in)  :: data_size
        type(c_ptr),        value,  intent(in)  :: mirror_data_ptrs
        type(c_ptr),        value,  intent(in)  :: ghost_data
    end subroutine
end interface

interface
# if defined(P4_TO_P8)
    subroutine p8est_balance(p8est, btype, init_fn) bind(c,name="p8est_balance")
# else
    subroutine p4est_balance(p4est, btype, init_fn) bind(c,name="p4est_balance")
# endif
        use iso_c_binding, only: c_ptr, c_int8_t

        type(c_ptr), value, intent(in)          :: p4est
        integer(c_int8_t), value, intent(in)    :: btype
        type(c_ptr), value, intent(in)          :: init_fn
    end subroutine
end interface

interface
# if defined(P4_TO_P8)
    subroutine p8est_transfer_fixed(dst_gfq, src_gfq, mpicomm, tag, dst_data, src_data, data_size) bind(c,name="p8est_transfer_fixed")
# else
    subroutine p4est_transfer_fixed(dst_gfq, src_gfq, mpicomm, tag, dst_data, src_data, data_size) bind(c,name="p4est_transfer_fixed")
# endif
        use iso_c_binding, only: c_ptr, c_int, c_size_t

        type(c_ptr), value, intent(in)          :: dst_gfq !! gfq ... global first quadrant
        type(c_ptr), value, intent(in)          :: src_gfq

# if defined(__INTEL_COMPILER)
        integer(c_int), value, intent(in)       :: mpicomm
# else
        type(c_ptr), value, intent(in)          :: mpicomm
# endif
        integer(c_int), value, intent(in)       :: tag

        type(c_ptr), value, intent(in)          :: dst_data
        type(c_ptr), value, intent(in)          :: src_data

        integer(c_size_t), value, intent(in)    :: data_size
    end subroutine
end interface

interface
# if defined(P4_TO_P8)
    pure subroutine p8est_qcoord_to_vertex(connectivity, treeid, x,y,z, coords) bind(c,name="p8est_qcoord_to_vertex")
# else
    pure subroutine p4est_qcoord_to_vertex(connectivity, treeid, x, y, coords) bind(c,name="p4est_qcoord_to_vertex")
# endif
        use iso_c_binding, only: c_ptr, c_int32_t, c_double

        type(c_ptr), value, intent(in)          :: connectivity
        integer(c_int32_t), value, intent(in)   :: treeid
# if defined(P4_TO_P8)
        integer(c_int32_t), value, intent(in)   :: x,y,z
# else
        integer(c_int32_t), value, intent(in)   :: x,y
# endif
        real(c_double), intent(inout)           :: coords(3)
    end subroutine
end interface

interface
# if defined(P4_TO_P8)
    pure function p8est_quadrant_is_family(q0,q1,q2,q3,q4,q5,q6,q7) bind(c,name="p8est_quadrant_is_family") result(yes)
# else
    pure function p4est_quadrant_is_family(q0,q1,q2,q3) bind(c,name="p4est_quadrant_is_family") result(yes)
# endif
        use iso_c_binding, only: c_int, c_ptr
        type(c_ptr), value, intent(in)  :: q0,q1,q2,q3
# if defined(P4_TO_P8)
        type(c_ptr), value, intent(in)  :: q4,q5,q6,q7
# endif
        integer(c_int)                  :: yes
    end function
end interface

interface
# if defined(P4_TO_P8)
    pure function p8est_quadrant_is_familypv(quads) bind(c,name="p8est_quadrant_is_familypv") result(yes)
# else
    pure function p4est_quadrant_is_familypv(quads) bind(c,name="p4est_quadrant_is_familypv") result(yes)
# endif
        use iso_c_binding, only: c_int, c_ptr
        type(c_ptr), value, intent(in)  :: quads
        integer(c_int)                  :: yes
    end function
end interface

interface
# if defined(P4_TO_P8)
    subroutine p8est_save_ext(filename, p8est, save_data, save_partition) bind(c,name="p8est_save_ext")
# else
    subroutine p4est_save_ext(filename, p4est, save_data, save_partition) bind(c,name="p4est_save_ext")
# endif
        use iso_c_binding, only: c_int, c_char, c_ptr
        type(c_ptr), value, intent(in)     :: filename
        type(c_ptr), value, intent(in)     :: p4est
        integer(c_int), value, intent(in)  :: save_data
        integer(c_int), value, intent(in)  :: save_partition 
    end subroutine
end interface

interface
# if defined(P4_TO_P8)
    subroutine p8est_ghost_expand(p8est,layer) bind(c,name="p8est_ghost_expand")
# else
    subroutine p4est_ghost_expand(p4est,layer) bind(c,name="p4est_ghost_expand")
# endif
        use iso_c_binding, only:  c_ptr
        type(c_ptr), value, intent(in)     :: p4est,layer
    end subroutine
end interface

end module
