module mesh_types_mod

use constants_mod
use mesh_const_mod
use morton_mod, only: MORTON_MAXLEVEL

use types_mod, only: payload_t
use types_mod, only: auxiliary_t

use iso_c_binding, only: c_ptr, c_int, c_size_t, c_int64_t

type :: quadptr_t

    type(quad_t), pointer :: p

end type

type :: sideptr_t

    type(side_t), pointer :: p

end type

type :: vertptr_t

    type(vert_t), pointer :: p

end type

type :: quad_t

    type(mesh_t), pointer   :: mesh

    integer(kind=1)         :: level
    integer(kind=4)         :: morton(N_DIMS)

    integer(kind=8)         :: locid,gloid

    type(payload_t)         :: pld
    type(auxiliary_t)       :: aux

    integer                 :: flags

    type(sideptr_t)         :: sides(MESH_FACES)
    type(vertptr_t)         :: verts(MESH_VERTS)

    integer(kind=8)         :: isides(MESH_FACES)
    integer(kind=8)         :: iverts(MESH_VERTS)

    real(dp)                :: delta(N_DIMS)
    real(dp)                :: center(N_DIMS)

    logical                 :: ismirror = .false.
    
end type

type :: side_t

    type(mesh_t), pointer       :: mesh

    integer                     :: direction    !! x,y
    integer                     :: faceid(2)    !! NOR,SOU,WES,EAS
    integer                     :: nquads(2)
    integer                     :: nsides

    integer(kind=8)             :: iquads(MESH_HALF,2)
    type(quadptr_t)             :: quads(MESH_HALF,2)

    integer(kind=8)             :: loc

    integer                     :: boundary
    integer                     :: inner,outer
    integer                     :: plus,minus

    logical                     :: isatzoneborder = .false.

end type

type :: vert_t

    type(mesh_t), pointer       :: mesh

    integer                     :: nquads = 0

    integer(kind=8)             :: iquads(MESH_VERTS) = 0
    type(quadptr_t)             :: quads(MESH_VERTS)

    integer(kind=8)             :: loc = 0

    integer                     :: boundary
    integer                     :: inner,outer
    integer                     :: plus,minus

    integer(kind=1)             :: level
    integer                     :: hanging
    integer                     :: refquad

    logical                     :: amr = .false.

end type

!! Humble beginnings of a proper octree
type :: tree_t
    
    integer(kind=8) :: timestamp = 0

    integer(kind=4) :: maxlevel

    integer(kind=8) :: numblocks

    integer(kind=8) :: nblocks(0:MORTON_MAXLEVEL)
    integer(kind=8) :: noffset(0:MORTON_MAXLEVEL)

    integer(kind=4), allocatable :: gid(:,:)
    integer(kind=4), allocatable :: morton(:,:)
    integer(kind=1), allocatable :: levels(:)

end type

type :: mesh_t

    integer(kind=1) :: minlevel, maxlevel

    integer(kind=8) :: ntotal,nquads,nghosts,nglobal,iglobal
    integer(kind=8) :: nsides,nverts,nverts_with_ghosts
    integer(kind=8) :: nmirrors

    ! integer(kind=8) :: noctree

    type(quad_t), allocatable   :: quads(:)
    type(side_t), allocatable   :: sides(:)
    type(vert_t), allocatable   :: verts(:)

    !! p4est specific fields
    type(c_ptr) :: p4est,ghosts,conn

    character(len=1), allocatable   :: mpi_sendbuf(:)
    character(len=1), allocatable   :: mpi_recvbuf(:)

    integer(c_size_t), allocatable  :: mirroridx(:)
    type(c_ptr), allocatable        :: mirrorptr(:)

    integer(kind=8) :: timestamp = 0

    integer(c_int64_t), allocatable :: gfq(:)

    integer(c_size_t), allocatable :: proc_offsets(:)
    integer(c_size_t), allocatable :: mirror_proc_offsets(:)
    integer(c_size_t), allocatable :: mirror_proc_mirrors(:)

    integer(c_size_t), allocatable :: proc_count(:)
    integer(c_size_t), allocatable :: mirror_proc_count(:)

    type(tree_t) :: tree

end type

abstract interface

    subroutine mesh_iterate_quads_cb_t(quad)

        import quad_t

        type(quad_t), intent(inout) :: quad

    end subroutine

    subroutine mesh_iterate_sides_cb_t(side)

        import side_t

        type(side_t), intent(inout) :: side

    end subroutine

    function mesh_probe_refine_cb_t(quad) result(ok)

        import quad_t

        type(quad_t), intent(in) :: quad
        logical                  :: ok

    end function

    function mesh_probe_coarsen_cb_t(family) result(ok)

        use mesh_const_mod, only: MESH_CHILDREN

        import quad_t

        type(quad_t), intent(in)    :: family(MESH_CHILDREN)
        logical                     :: ok

    end function

    subroutine mesh_replace_refine_cb_t(incoming, outgoing)

        use mesh_const_mod, only: MESH_CHILDREN

        import quad_t

        type(quad_t), intent(in)        :: incoming
        type(quad_t), intent(inout)     :: outgoing(MESH_CHILDREN)

    end subroutine

    subroutine mesh_replace_coarsen_cb_t(incoming, outgoing)

        use mesh_const_mod, only: MESH_CHILDREN

        import quad_t

        type(quad_t), intent(in)    :: incoming(MESH_CHILDREN)
        type(quad_t), intent(inout) :: outgoing

    end subroutine

end interface

end module
