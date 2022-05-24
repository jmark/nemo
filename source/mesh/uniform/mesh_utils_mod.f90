!# define PP_MESH_PERIODIC 1

module mesh_utils_mod

use mpi_f08
use mesh_types_mod

# if DO_PROFILING
    use profiler_mod, only: profiler_init,profiler_start, profiler_stop
# endif

interface mesh_iterate_quads
    procedure :: mesh_iterate_quads_without_flags
    procedure :: mesh_iterate_quads_with_flags
end interface

interface mesh_iterate_sides
    procedure :: mesh_iterate_sides_without_flags
    procedure :: mesh_iterate_sides_with_flags
end interface

interface mesh_refine
    procedure :: mesh_refine_without_replace
    procedure :: mesh_refine_with_replace
end interface

# if DO_PROFILING
    integer, save :: prof_mesh_exchange
    integer, save :: prof_mesh_adapt
    integer, save :: prof_mesh_partition
# endif

contains

subroutine mesh_init(mesh)

    use share_mod, only: rt => runtime
    
    use sc_constants_mod, only: SC_LP_ERROR
    use p4est_interfaces_mod, only: p4est_init
    use p4est_interfaces_mod, only: p4est_new_ext

    use p4est_interfaces_mod, only: p4est_connectivity_new_periodic
    use p4est_interfaces_mod, only: p4est_connectivity_new_unitsquare

    use setup_config_mod, only: mesh_minlevel, mesh_maxlevel
    use iso_c_binding, only: c_ptr, c_null_ptr, c_loc

    type(mesh_t), intent(inout), target :: mesh

# if PP_MESH_SET_LEVELS_AT_RUNTIME
    CHARACTER(len=32) :: nemo_mesh_minlevel
    CHARACTER(len=32) :: nemo_mesh_maxlevel
            
    CALL getenv("NEMO_MESH_MINLEVEL", nemo_mesh_minlevel)
    CALL getenv("NEMO_MESH_MAXLEVEL", nemo_mesh_maxlevel)

    read(nemo_mesh_minlevel,*) mesh%minlevel
    read(nemo_mesh_maxlevel,*) mesh%maxlevel
# else
    mesh%minlevel = mesh_minlevel
    mesh%maxlevel = mesh_maxlevel
# endif

    call p4est_init(c_null_ptr, SC_LP_ERROR)

# if PP_MESH_PERIODIC
    mesh%conn = p4est_connectivity_new_periodic()
# else
    mesh%conn = p4est_connectivity_new_unitsquare()
# endif

    mesh%p4est = p4est_new_ext(&
         mpicomm         = rt%mpi%ccomm, &
         connectivity    = mesh%conn, &
         min_quads       = 0, &
         min_level       = INT(mesh%minlevel,kind=4), &
         fill_uniform    = 1, &
         data_size       = INT(0,kind=8), &
         init_fn         = c_null_ptr, &
         userptr         = c_null_ptr)
    
    call mesh_alloc(mesh)

# if DO_PROFILING
    prof_mesh_exchange  = profiler_init('mesh_exchange')
    prof_mesh_adapt     = profiler_init('mesh_adapt')
    prof_mesh_partition = profiler_init('mesh_partition')
# endif

end subroutine

subroutine mesh_nuke(mesh)

    use p4est_interfaces_mod, only: p4est_destroy
    use p4est_interfaces_mod, only: p4est_connectivity_destroy

    type(mesh_t), intent(inout), target :: mesh

    call mesh_dealloc(mesh)
    call p4est_connectivity_destroy(mesh%conn)
    call p4est_destroy(mesh%p4est)

end subroutine

subroutine mesh_alloc(mesh)

    use share_mod, only: rt => runtime

    use iso_c_binding, only: c_associated
    use iso_c_binding, only: c_loc, c_f_pointer

    use p4est_interfaces_mod, only: p4est_ghost_new

    use p4wrap_interfaces_mod, only: p4wrap_gather_mirroridx
    use p4wrap_interfaces_mod, only: p4wrap_gather_connection_info
    use p4wrap_interfaces_mod, only: p4wrap_gather_proc_offsets
    use p4wrap_interfaces_mod, only: p4wrap_gather_mirror_proc_offsets
    use p4wrap_interfaces_mod, only: p4wrap_gather_mirror_proc_mirrors

    use p4wrap_types_mod, only: p4wrap_quad_t
    use p4wrap_types_mod, only: p4wrap_side_t
    use p4wrap_types_mod, only: p4wrap_vert_t

    use p4est_constants_mod, only: P4EST_CONNECT_FULL

    use p4wrap_interfaces_mod, only: p4wrap_get_num_ghosts
    use p4wrap_interfaces_mod, only: p4wrap_get_num_mirrors
    use p4wrap_interfaces_mod, only: p4wrap_get_num_quads
    use p4wrap_interfaces_mod, only: p4wrap_get_num_global
    use p4wrap_interfaces_mod, only: p4wrap_get_num_sides
    use p4wrap_interfaces_mod, only: p4wrap_get_num_verts
    use p4wrap_interfaces_mod, only: p4wrap_get_num_verts_hanging
    use p4wrap_interfaces_mod, only: p4wrap_get_idx_global
    use p4wrap_interfaces_mod, only: p4wrap_copy_gfq

    use types_mod, only: payload_size

    use mesh_const_mod, only: map_to_opposite_vertid
    use mesh_const_mod, only: map_faceid_to_hanging_vertid_verts
    use mesh_const_mod, only: map_faceid_to_hanging_vertid_quads

    use mesh_const_mod, only: map_faceid_to_direction
    use mesh_const_mod, only: map_faceid_to_opposite

    use quad_utils_mod, only: vert_get_coords

    use quad_utils_mod, only: quad_get_delta
    use quad_utils_mod, only: quad_get_center

    use status_mod, only: raise_error

    use mpi_f08, only: MPI_IN_PLACE
    use mpi_f08, only: MPI_INTEGER1
    use mpi_f08, only: MPI_Allgatherv
    use mpi_f08, only: MPI_DATATYPE_NULL

    use types_mod, only: payload_dummy, auxiliary_dummy

    use setup_config_mod
 
    type(mesh_t), intent(inout), target :: mesh

    type(p4wrap_quad_t), allocatable, target :: quads(:)
    type(p4wrap_side_t), allocatable, target :: sides(:)
    type(p4wrap_vert_t), allocatable, target :: verts(:)

    integer(kind=8) :: i, iproc

    integer :: iquad,iside,ivert

    mesh%ghosts   = p4est_ghost_new(mesh%p4est,P4EST_CONNECT_FULL)
    mesh%nghosts  = p4wrap_get_num_ghosts(mesh%ghosts)
    mesh%nmirrors = p4wrap_get_num_mirrors(mesh%ghosts)

    allocate(mesh%mirroridx(mesh%nmirrors))
    allocate(mesh%mirrorptr(mesh%nmirrors))

    call p4wrap_gather_mirroridx(mesh%ghosts,c_loc(mesh%mirroridx))

    !! C indexing -> Fortran indexing
    mesh%mirroridx = mesh%mirroridx + 1

    do i = 1,mesh%nmirrors
        mesh%mirrorptr(i) = c_loc(mesh%mpi_sendbuf(1+(i-1)*payload_size))
    end do

    mesh%iglobal = p4wrap_get_idx_global(mesh%p4est)
    mesh%nglobal = p4wrap_get_num_global(mesh%p4est)
    mesh%nquads = p4wrap_get_num_quads(mesh%p4est)
    mesh%nsides = p4wrap_get_num_sides(mesh%p4est,mesh%ghosts)
    mesh%ntotal = mesh%nquads + mesh%nghosts
    mesh%nverts = p4wrap_get_num_verts(mesh%p4est,mesh%ghosts)

    allocate(mesh%gfq(rt%mpi%size+1))
    call p4wrap_copy_gfq(mesh%p4est,size(mesh%gfq),c_loc(mesh%gfq))

    allocate(mesh%proc_offsets(rt%mpi%size+1))
    call p4wrap_gather_proc_offsets(mesh%ghosts,size(mesh%proc_offsets),c_loc(mesh%proc_offsets))

    allocate(mesh%mirror_proc_offsets(rt%mpi%size+1))
    call p4wrap_gather_mirror_proc_offsets(mesh%ghosts,size(mesh%mirror_proc_offsets),c_loc(mesh%mirror_proc_offsets))

    allocate(mesh%mirror_proc_mirrors(mesh%mirror_proc_offsets(rt%mpi%size+1)))
    call p4wrap_gather_mirror_proc_mirrors(mesh%ghosts,size(mesh%mirror_proc_mirrors),c_loc(mesh%mirror_proc_mirrors))

    allocate(mesh%proc_count(rt%mpi%size))
    allocate(mesh%mirror_proc_count(rt%mpi%size))

    do i = 1,size(mesh%mirror_proc_count)
        mesh%proc_count(i) = mesh%proc_offsets(i+1) - mesh%proc_offsets(i)
    end do

    do i = 1,size(mesh%mirror_proc_count)
        mesh%mirror_proc_count(i) = mesh%mirror_proc_offsets(i+1) - mesh%mirror_proc_offsets(i)
    end do

    allocate(mesh%mpi_sendbuf(payload_size*mesh%mirror_proc_offsets(rt%mpi%size+1)))
    allocate(mesh%mpi_recvbuf(payload_size*mesh%nghosts))

    allocate(quads(mesh%ntotal))
    allocate(sides(mesh%nsides))
    allocate(verts(mesh%nverts))

    call p4wrap_gather_connection_info(mesh%p4est,mesh%ghosts,c_loc(quads),c_loc(sides),c_loc(verts))

    allocate(mesh%quads(mesh%ntotal))
    allocate(mesh%sides(mesh%nsides))
    allocate(mesh%verts(mesh%nverts + MESH_VERTS*mesh%nghosts))

    do i = 1,mesh%ntotal
        mesh%quads(i)%pld = payload_dummy
        mesh%quads(i)%aux = auxiliary_dummy
    end do

    do i = 1,mesh%nquads
        associate(fquad => mesh%quads(i), cquad => quads(i))
        fquad%mesh => mesh

        fquad%level = cquad%level
        fquad%morton = cquad%morton

        fquad%isides = cquad%isides+1
        fquad%iverts = cquad%iverts+1

        fquad%delta = quad_get_delta(fquad)
        fquad%center = quad_get_center(fquad)

        fquad%locid = i
        fquad%gloid = i + mesh%gfq(rt%mpi%rank+1)

        fquad%ismirror = .false.
        end associate
    end do

    do i = 1,size(mesh%mirroridx)
        mesh%quads(mesh%mirroridx(i))%ismirror = .true.
    end do

    iproc = 1
    do i = 1,mesh%nghosts
        associate(fquad => mesh%quads(i+mesh%nquads),cquad => quads(i+mesh%nquads))
        fquad%mesh => mesh

        fquad%level = cquad%level
        fquad%morton = cquad%morton

        fquad%isides = cquad%isides+1
        fquad%iverts = cquad%iverts+1

        fquad%locid = i+mesh%nquads

        do while (i > mesh%proc_offsets(iproc+1)) 
            iproc = iproc + 1
        end do

        fquad%gloid = cquad%local_num + mesh%gfq(iproc) + 1

        fquad%delta = quad_get_delta(fquad)
        fquad%center = quad_get_center(fquad)

        fquad%ismirror = .false.
        end associate
    end do

    do i = 1,mesh%nsides
        associate(fside => mesh%sides(i), cside => sides(i))
        fside%loc = i
        fside%mesh => mesh

        fside%nquads = cside%nquads
        fside%faceid = cside%faceid+1
        fside%iquads = cside%iquads+1

        fside%isatzoneborder = any(fside%iquads > mesh%nquads)

        do iside = 1,2
            do iquad = 1,fside%nquads(iside)
                fside%quads(iquad,iside)%p => mesh%quads(cside%iquads(iquad,iside)+1)
            end do
        end do

        if (cside%tree_boundary > 0 .and. cside%nquads(2) < 1) then
            fside%faceid(2) = map_faceid_to_opposite(fside%faceid(1))

            fside%boundary = fside%faceid(1)

            fside%inner = 1
            fside%outer = 2
        else
            fside%boundary = 0

            fside%inner = -1
            fside%outer = -1
        end if

        fside%direction = map_faceid_to_direction(minval(fside%faceid))

        fside%plus  = maxloc(fside%faceid,DIM=1)
        fside%minus = minloc(fside%faceid,DIM=1)
        end associate
    end do

    do i = 1,mesh%nverts
        associate(fvert => mesh%verts(i), cvert => verts(i))
        fvert%loc = i
        fvert%mesh => mesh

        fvert%nquads = cvert%nquads
        fvert%iquads = cvert%iquads+1

        fvert%level = 0
        do ivert = 1,MESH_VERTS
            if (fvert%iquads(ivert) > 0) then
                fvert%level = max(fvert%level,mesh%quads(fvert%iquads(ivert))%level)
            end if
        end do

        fvert%amr = .false.
        do ivert = 1,MESH_VERTS
            if (fvert%iquads(ivert) > 0) then
                if (fvert%level /= mesh%quads(fvert%iquads(ivert))%level) then
                    fvert%amr = .true.
                    exit 
                end if
            end if
        end do

        fvert%hanging = 0
        fvert%refquad = maxloc(fvert%iquads,DIM=1)
        end associate
    end do

    block
    integer(kind=8) :: iv_incr
    integer :: iv
    iv_incr = mesh%nverts
    do i = mesh%nquads+1,mesh%nquads+mesh%nghosts
        associate(quad => mesh%quads(i))
        do iv = 1,MESH_VERTS
            if (quad%iverts(iv) > 0) cycle

            iv_incr = iv_incr + 1
            quad%iverts(iv) = iv_incr

            associate(vert => mesh%verts(iv_incr))
            vert%loc = iv_incr
            vert%mesh => mesh
            vert%level = quad%level

            vert%nquads = 1
            vert%iquads(map_to_opposite_vertid(iv)) = i
            vert%refquad = maxloc(vert%iquads,DIM=1)
            end associate
        end do
        end associate
    end do
    mesh%nverts_with_ghosts = iv_incr
    end block

    deallocate(quads)
    deallocate(sides)
    deallocate(verts)

    call SYSTEM_CLOCK(mesh%timestamp)

end subroutine

subroutine mesh_alloc_tree(mesh)

    use share_mod, only: rt => runtime

    use p4est_constants_mod, only: P4EST_CONNECT_FULL

    use mpi_f08, only: MPI_IN_PLACE
    use mpi_f08, only: MPI_INTEGER1
    use mpi_f08, only: MPI_INTEGER
    use mpi_f08, only: MPI_Allgatherv
    use mpi_f08, only: MPI_DATATYPE_NULL

    use morton_mod, only: morton_is_family
    use morton_mod, only: morton_get_parent
    use morton_mod, only: morton_is_ancestor

    type(mesh_t), intent(inout), target :: mesh

    integer(kind=4) :: level,l, morton(N_DIMS)

    integer(kind=4) :: QUAD_LEN
    integer(kind=4) :: i,j,k,d
    integer(kind=4) :: ii,jj,kk
    integer(kind=4) :: zorder

    integer(kind=8) :: iquad,jquad,kquad,incr

    integer(kind=1), allocatable :: alllevels(:)
    integer(kind=4), allocatable :: allmorton(:,:)

    if (mesh%tree%timestamp > mesh%timestamp) return

    allocate(alllevels(       mesh%nglobal))
    allocate(allmorton(N_DIMS,mesh%nglobal))

    do iquad = 1,mesh%nquads
        alllevels(  mesh%gfq(rt%mpi%rank+1)+iquad) = mesh%quads(iquad)%level
        allmorton(:,mesh%gfq(rt%mpi%rank+1)+iquad) = mesh%quads(iquad)%morton
    end do

    if (rt%mpi%size > 1) then
        call MPI_Allgatherv(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,alllevels,&
            INT(mesh%gfq(2:size(mesh%gfq))-mesh%gfq(1:size(mesh%gfq)-1),kind=4),&
            INT(mesh%gfq,kind=4),MPI_INTEGER1,rt%mpi%comm)

        call MPI_Allgatherv(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,allmorton,&
            N_DIMS*INT(mesh%gfq(2:size(mesh%gfq))-mesh%gfq(1:size(mesh%gfq)-1),kind=4),&
            N_DIMS*INT(mesh%gfq,kind=4),MPI_INTEGER,rt%mpi%comm)
    end if

    mesh%tree%maxlevel = maxval(alllevels)

    mesh%tree%nblocks = 0
    mesh%tree%nblocks(0) = 1

    do level = 1,mesh%tree%maxlevel-1
        incr = 0
        do zorder = 1,(2**level)**N_DIMS
            morton = 0
            do l = 0,level
                do d = 1,N_DIMS
                    morton(d) = morton(d) + SHIFTL(SHIFTR(iand(zorder-1,SHIFTL(1,N_DIMS*l+d-1)),N_DIMS*l+d-1),MORTON_MAXLEVEL-level+l)
                end do
            end do
            do iquad = 1,mesh%nglobal
                if (morton_is_ancestor(level,morton,INT(alllevels(iquad),kind=4),allmorton(:,iquad))) then
                    incr = incr + 1
                    exit
                end if
            end do
        end do
        mesh%tree%nblocks(level) = incr
    end do

    mesh%tree%nblocks(mesh%tree%maxlevel) = mesh%nglobal
    mesh%tree%numblocks = sum(mesh%tree%nblocks)

    mesh%tree%noffset = 0
    do level = 1,mesh%tree%maxlevel+1
        do i = 0,level-1
            mesh%tree%noffset(level) = mesh%tree%noffset(level) + mesh%tree%nblocks(i)
        end do
    end do

    if (allocated(mesh%tree%morton)) deallocate(mesh%tree%morton)
    if (allocated(mesh%tree%levels)) deallocate(mesh%tree%levels)

    allocate(mesh%tree%morton(N_DIMS,mesh%tree%numblocks))
    allocate(mesh%tree%levels(       mesh%tree%numblocks))

    mesh%tree%levels(  1) = 0
    mesh%tree%morton(:,1) = 0

    do iquad = 1,mesh%nglobal
        mesh%tree%levels(  mesh%tree%noffset(mesh%tree%maxlevel)+iquad) = alllevels(  iquad)
        mesh%tree%morton(:,mesh%tree%noffset(mesh%tree%maxlevel)+iquad) = allmorton(:,iquad)
    end do

    do level = 1,mesh%tree%maxlevel-1
        incr = 0
        do zorder = 1,(2**level)**N_DIMS
            morton = 0
            do l = 0,level
                do d = 1,N_DIMS
                    morton(d) = morton(d) + SHIFTL(SHIFTR(iand(zorder-1,SHIFTL(1,N_DIMS*l+d-1)),N_DIMS*l+d-1),MORTON_MAXLEVEL-level+l)
                end do
            end do
            do iquad = 1,mesh%nglobal
                if (morton_is_ancestor(level,morton,INT(alllevels(iquad),kind=4),allmorton(:,iquad))) then
                    incr = incr + 1
                    mesh%tree%levels(  mesh%tree%noffset(level)+incr) = INT(level,kind=1)
                    mesh%tree%morton(:,mesh%tree%noffset(level)+incr) = morton
                    exit
                end if
            end do
        end do
    end do

    deallocate(allmorton)
    deallocate(alllevels)

    if (allocated(mesh%tree%gid)) deallocate(mesh%tree%gid)
    allocate(mesh%tree%gid(MESH_FACES+1+MESH_CHILDREN,mesh%tree%numblocks))
    mesh%tree%gid = -1

    do iquad = 1,mesh%tree%noffset(mesh%tree%maxlevel)
        incr = 0
        do jquad = iquad+1,mesh%tree%numblocks
            if (mesh%tree%levels(iquad)+1 == mesh%tree%levels(jquad)) then
                if (morton_is_ancestor(&
                        INT(mesh%tree%levels(iquad),kind=4),mesh%tree%morton(:,iquad),&
                        INT(mesh%tree%levels(jquad),kind=4),mesh%tree%morton(:,jquad))) then

                    incr = incr + 1
                    mesh%tree%gid(MESH_FACES+1+incr,iquad) = INT(jquad,kind=4)
                    mesh%tree%gid(MESH_FACES+1     ,jquad) = INT(iquad,kind=4)
                    if (incr >= MESH_CHILDREN) exit
                end if
            end if
        end do
    end do

    call SYSTEM_CLOCK(mesh%tree%timestamp)

end subroutine

subroutine mesh_dealloc(mesh)

    use p4est_interfaces_mod, only: p4est_ghost_destroy

    type(mesh_t), intent(inout) :: mesh

    deallocate(mesh%quads)
    deallocate(mesh%sides)
    deallocate(mesh%verts)

    deallocate(mesh%mirroridx)
    deallocate(mesh%mirrorptr)
    deallocate(mesh%mpi_sendbuf)
    deallocate(mesh%mpi_recvbuf)

    deallocate(mesh%gfq)
    deallocate(mesh%proc_offsets)
    deallocate(mesh%mirror_proc_offsets)
    deallocate(mesh%mirror_proc_mirrors)

    deallocate(mesh%proc_count)
    deallocate(mesh%mirror_proc_count)

    mesh%nquads = 0
    mesh%nsides = 0
    mesh%ntotal = 0

    mesh%nghosts = 0
    mesh%nmirrors = 0

    call p4est_ghost_destroy(mesh%ghosts)

end subroutine

subroutine mesh_exchange(mesh)

    use iso_c_binding, only: c_loc
    use share_mod, only: rt => runtime
    use types_mod, only: payload_size, payload_dummy
    use mpi_f08, only: mpi_ialltoallv
    use mpi_f08, only: mpi_byte
    use mpi_f08, only: mpi_request
    use mpi_f08, only: mpi_status
    use mpi_f08, only: mpi_wait

    type(mesh_t), intent(inout), target :: mesh
    
    integer(kind=8) :: q
    character(len=1) :: char_dummy

    integer, parameter :: byte_size = STORAGE_SIZE(char_dummy) / 8

    integer :: ierr

    type(MPI_Request) :: req
    type(MPI_Status) :: stat

    if (rt%mpi%size > 1) then
# if DO_PROFILING
        call profiler_start(prof_mesh_exchange)
# endif

        !$OMP PARALLEL DO default(shared) private(q)
        do q = 1,size(mesh%mirror_proc_mirrors)
            mesh%mpi_sendbuf(1+(q-1)*payload_size:q*payload_size) = transfer(&
                SOURCE=mesh%quads(mesh%mirroridx(mesh%mirror_proc_mirrors(q) + 1))%pld,MOLD=char_dummy,SIZE=payload_size)
        end do
        !$OMP end parallel do

        call MPI_IAlltoallv(&
            mesh%mpi_sendbuf,INT(payload_size)*INT(mesh%mirror_proc_count),INT(payload_size)*INT(mesh%mirror_proc_offsets(1:rt%mpi%size)), MPI_BYTE, &
            mesh%mpi_recvbuf,INT(payload_size)*INT(mesh%proc_count),INT(payload_size)*INT(mesh%proc_offsets(1:rt%mpi%size)), MPI_BYTE, rt%mpi%comm, req, ierr)
        call mpi_wait(req,stat,ierr)

        !$OMP PARALLEL DO default(shared) private(q)
        do q = 1,mesh%nghosts
            mesh%quads(mesh%nquads+q)%pld = transfer(&
                SOURCE=mesh%mpi_recvbuf(1+(q-1)*payload_size:q*payload_size),MOLD=payload_dummy)
        end do
        !$OMP end parallel do

# if DO_PROFILING
        call profiler_stop(prof_mesh_exchange)
# endif
    end if

end subroutine

subroutine mesh_iterate_quads_without_flags(mesh,iterate_cb)

    type(mesh_t), intent(inout), target :: mesh
    procedure(mesh_iterate_quads_cb_t)  :: iterate_cb

    call mesh_iterate_quads_with_flags(mesh,iterate_cb,flags=0)

end subroutine

subroutine mesh_iterate_quads_with_flags(mesh,iterate_cb,flags)

# ifdef _OPENMP
    USE OMP_LIB, only: omp_get_thread_num
    USE OMP_LIB, only: omp_get_num_threads
# endif

    type(mesh_t), intent(inout), target :: mesh
    procedure(mesh_iterate_quads_cb_t)  :: iterate_cb
    integer                             :: flags

    integer(kind=8) :: i

    if (iand(flags,MESH_ITERATE_WITH_GHOSTS) > 0) then
        !$OMP parallel do private(i) shared(mesh) schedule(dynamic)
        do i = 1,mesh%ntotal
            call iterate_cb(mesh%quads(i))
        end do
        !$OMP END parallel DO

    else if (iand(flags,MESH_ITERATE_ONLY_INNER) > 0) then
        !$OMP parallel do private(i) shared(mesh) schedule(dynamic)
        do i = 1,mesh%nquads
            ! if (i == 1) write (*,*) 'iter only inner', omp_get_thread_num(), omp_get_num_threads()
            !!$OMP task default(shared) firstprivate(i)
            call iterate_cb(mesh%quads(i))
            !!$OMP END TASK
        end do
        !$OMP END parallel DO

    else if (iand(flags,MESH_ITERATE_ONLY_GHOSTS) > 0) then
        !$OMP parallel do private(i) shared(mesh) schedule(dynamic)
        do i = mesh%nquads+1,mesh%ntotal
            call iterate_cb(mesh%quads(i))
        end do
        !$OMP END parallel DO

    else if (iand(flags,MESH_ITERATE_ONLY_INNER_WITHOUT_MIRRORS) > 0) then
        !$OMP parallel do private(i) shared(mesh) schedule(dynamic)
        do i = 1,mesh%nquads
            if (.not.mesh%quads(i)%ismirror) then
                call iterate_cb(mesh%quads(i))
            end if
        end do
        !$OMP END parallel DO

    else
        !$OMP parallel do private(i) shared(mesh) schedule(dynamic)
        do i = 1,mesh%nquads
            call iterate_cb(mesh%quads(i))
        end do
        !$OMP END parallel DO
    end if

end subroutine

subroutine mesh_iterate_sides_without_flags(mesh,iterate_cb)

    type(mesh_t), intent(inout), target :: mesh
    procedure(mesh_iterate_sides_cb_t)  :: iterate_cb

    call mesh_iterate_sides_with_flags(mesh,iterate_cb,flags=0)

end subroutine

subroutine mesh_iterate_sides_with_flags(mesh,iterate_cb,flags)

# ifdef _OPENMP
    USE OMP_LIB, only: omp_get_thread_num
    USE OMP_LIB, only: omp_get_num_threads
# endif

    type(mesh_t), intent(inout), target :: mesh
    procedure(mesh_iterate_sides_cb_t)  :: iterate_cb
    integer                             :: flags

    integer(kind=8) :: i

    if (iand(flags,MESH_ITERATE_ONLY_INNER) > 0) then
        !$OMP PARALLEL DO private(i) shared(mesh) schedule(dynamic)
        do i = 1,mesh%nsides
            if (.not.mesh%sides(i)%isatzoneborder) then
                !!$OMP task default(shared) firstprivate(i)
                call iterate_cb(mesh%sides(i))
                !write (*,*) 'task', i, omp_get_thread_num(), omp_get_num_threads()
                !!$OMP END TASK
                !write (*,*) 'iter only inner', i, omp_get_thread_num(), omp_get_num_threads()
            end if
        end do
        !$OMP END PARALLEL DO

    else if (iand(flags,MESH_ITERATE_ONLY_ZONE_BOUNDARY) > 0) then
        !$OMP PARALLEL DO private(i) shared(mesh) schedule(dynamic)
        do i = 1,mesh%nsides
            if (mesh%sides(i)%isatzoneborder) then
                call iterate_cb(mesh%sides(i))
            end if
        end do
        !$OMP END PARALLEL DO

    else
        !$OMP PARALLEL DO private(i) shared(mesh) schedule(dynamic)
        do i = 1,mesh%nsides
            call iterate_cb(mesh%sides(i))
        end do
        !$OMP END PARALLEL DO
    end if

end subroutine

function dummy_probe_refine_cb(quad) result(ok)

    type(quad_t), intent(in)    :: quad
    logical                     :: ok

    ok = .false.

end function

function dummy_probe_coarsen_cb(family) result(ok)

    type(quad_t), intent(in)    :: family(MESH_CHILDREN)
    logical                     :: ok

    ok = .false.

end function

subroutine dummy_replace_refine_cb(incoming, outgoing)

    type(quad_t), intent(in)    :: incoming
    type(quad_t), intent(inout) :: outgoing(MESH_CHILDREN)

end subroutine

subroutine dummy_replace_coarsen_cb(incoming, outgoing)

    type(quad_t), intent(in)    :: incoming(MESH_CHILDREN)
    type(quad_t), intent(inout) :: outgoing

end subroutine

subroutine mesh_refine_without_replace(mesh,probe_refine_cb)

    type(mesh_t), intent(inout), target     :: mesh
    procedure(mesh_probe_refine_cb_t)       :: probe_refine_cb

    call mesh_adapt(mesh,probe_refine_cb,dummy_replace_refine_cb,dummy_probe_coarsen_cb,dummy_replace_coarsen_cb)

end subroutine

subroutine mesh_refine_with_replace(mesh,probe_refine_cb,replace_refine_cb)

    type(mesh_t), intent(inout), target     :: mesh
    procedure(mesh_probe_refine_cb_t)       :: probe_refine_cb
    procedure(mesh_replace_refine_cb_t)     :: replace_refine_cb

    call mesh_adapt(mesh,probe_refine_cb,replace_refine_cb,dummy_probe_coarsen_cb,dummy_replace_coarsen_cb)

end subroutine

subroutine mesh_coarsen(mesh,probe_coarsen_cb,replace_coarsen_cb)

    type(mesh_t), intent(inout), target     :: mesh
    procedure(mesh_probe_coarsen_cb_t)      :: probe_coarsen_cb
    procedure(mesh_replace_coarsen_cb_t)    :: replace_coarsen_cb

    call mesh_adapt(mesh,dummy_probe_refine_cb,dummy_replace_refine_cb,probe_coarsen_cb,replace_coarsen_cb)

end subroutine

subroutine mesh_adapt(mesh,probe_refine_cb,replace_refine_cb,probe_coarsen_cb,replace_coarsen_cb)

    type(mesh_t), intent(inout), target     :: mesh

    procedure(mesh_probe_refine_cb_t)       :: probe_refine_cb
    procedure(mesh_replace_refine_cb_t)     :: replace_refine_cb

    procedure(mesh_probe_coarsen_cb_t)      :: probe_coarsen_cb
    procedure(mesh_replace_coarsen_cb_t)    :: replace_coarsen_cb

# if DO_PROFILING
    call profiler_start(prof_mesh_adapt)
# endif

    call mesh_adapt_wrapped(mesh,probe_refine_cb,replace_refine_cb,probe_coarsen_cb,replace_coarsen_cb)

# if DO_PROFILING
    call profiler_stop(prof_mesh_adapt)
# endif

end subroutine

subroutine mesh_adapt_wrapped(mesh,probe_refine_cb,replace_refine_cb,probe_coarsen_cb,replace_coarsen_cb)

    use quad_utils_mod, only: quad_is_equal
    use iso_c_binding, only: c_associated
    use iso_c_binding, only: c_int8_t, c_loc, c_null_ptr

    use p4wrap_interfaces_mod, only: p4wrap_refine
    use p4wrap_interfaces_mod, only: p4wrap_coarsen
    use p4wrap_interfaces_mod, only: p4wrap_balance
    use p4est_interfaces_mod, only: p4est_quadrant_is_family

    use quad_utils_mod, only: quad_is_family
    use share_mod, only: mesh_ptr => mesh, meshes, mesh_active, rt => runtime

    use mpi_f08, only: MPI_Allreduce,MPI_SUM,MPI_INTEGER

    type(mesh_t), intent(inout), target     :: mesh

    procedure(mesh_probe_refine_cb_t)       :: probe_refine_cb
    procedure(mesh_replace_refine_cb_t)     :: replace_refine_cb

    procedure(mesh_probe_coarsen_cb_t)      :: probe_coarsen_cb
    procedure(mesh_replace_coarsen_cb_t)    :: replace_coarsen_cb

    type(mesh_t), pointer                   :: old,new

    integer(c_int8_t), allocatable, target :: flags(:)

    integer(kind=8) :: i,inew,iold
    integer         :: locstat,glostat !! local and global status code

    if (mesh%minlevel == mesh%maxlevel) return

    allocate(flags(mesh%nquads))
    flags = 0
    locstat = 0
    glostat = 0

    !$OMP PARALLEL DO default(shared) private(i) reduction(+:locstat)
    do i = 1,mesh%nquads
        if (mesh%quads(i)%level < mesh%maxlevel) then
            if (probe_refine_cb(mesh%quads(i))) then
                flags(i) = MESH_FLAG_REFINE
                locstat = 1
                cycle
            end if
        end if

        if (mesh%quads(i)%level <= mesh%minlevel) cycle
        if (mesh%quads(i)%locid + MESH_CHILDREN - 1 > mesh%nquads) cycle

        if (quad_is_family(mesh%quads(i:i+MESH_CHILDREN-1))) then
            if (probe_coarsen_cb(mesh%quads(i:i+MESH_CHILDREN-1))) then
                flags(i) = MESH_FLAG_COARSEN
                locstat = 1
            end if
        end if
    end do
    !$OMP END PARALLEL DO

    if (rt%mpi%size > 1) then
        call MPI_Allreduce(locstat,glostat,1,MPI_INTEGER,MPI_SUM,rt%mpi%comm)
    else
        glostat = locstat
    end if

    if (glostat == 0) then
        deallocate(flags)
        return
    end if

    !! ---------------------------------------------------------------------- !!

    old => meshes(mesh_active)

    mesh_active = 3-mesh_active
    mesh_ptr => meshes(mesh_active)

    new => meshes(mesh_active)

    new%p4est = old%p4est
    new%conn  = old%conn
    new%minlevel = old%minlevel
    new%maxlevel = old%maxlevel

    !! side effects on 'new'
    call p4wrap_coarsen(new%p4est,c_loc(flags))
    call p4wrap_refine(new%p4est,c_loc(flags))
    call p4wrap_balance(new%p4est)

    deallocate(flags)

    call mesh_alloc(new)

    !! ---------------------------------------------------------------------- !!

    iold = 1
    inew = 1

    !$OMP PARALLEL default(shared)
    !$OMP single
    do while (iold <= old%nquads .and. inew <= new%nquads)
        if (old%quads(iold)%level < new%quads(inew)%level) then
            !$OMP task default(shared) firstprivate(iold,inew)
            call replace_refine_cb(old%quads(iold),new%quads(inew:inew+MESH_CHILDREN-1))
            !$OMP end task
            iold = iold + 1
            inew = inew + MESH_CHILDREN

        else if (old%quads(iold)%level > new%quads(inew)%level) then
            !$OMP task default(shared) firstprivate(iold,inew)
            call replace_coarsen_cb(old%quads(iold:iold+MESH_CHILDREN-1),new%quads(inew))
            !$OMP end task
            iold = iold + MESH_CHILDREN
            inew = inew + 1

        else
            new%quads(inew)%pld = old%quads(iold)%pld
            new%quads(inew)%aux = old%quads(iold)%aux
            iold = iold + 1
            inew = inew + 1
        end if
    end do
    !$OMP end single
    !$OMP end parallel

    !! ---------------------------------------------------------------------- !!

    call mesh_dealloc(old)

end subroutine

!! TODO: implement weighting function support
!! TODO: Not tested yet.

subroutine mesh_partition(mesh)

    use iso_c_binding, only: c_int64_t, c_loc

    use types_mod, only: payload_dummy, payload_size

    use p4est_interfaces_mod, only: p4est_transfer_fixed
    use p4wrap_interfaces_mod, only: p4wrap_partition
    use p4wrap_interfaces_mod, only: p4wrap_copy_gfq

    use share_mod, only: rt => runtime

    type(mesh_t), intent(inout), target :: mesh

    integer(kind=8), save :: timestamp = 0

    character(len=1)                        :: char_dummy
    character(len=1), allocatable, target   :: old_payloads(:)
    character(len=1), allocatable, target   :: new_payloads(:)

    integer(c_int64_t), allocatable, target :: old_gfq(:)

    integer(kind=8) :: i

    if (rt%mpi%size < 2) return
    if (mesh%timestamp <= timestamp) return

# if DO_PROFILING
    call profiler_start(prof_mesh_partition)
# endif

    old_gfq = mesh%gfq
    allocate(old_payloads(payload_size*mesh%nquads))
    
    !$OMP PARALLEL DO default(shared) private(i)
    do i = 1,mesh%nquads
        old_payloads(1+(i-1)*payload_size:i*payload_size) = transfer(&
            SOURCE=mesh%quads(i)%pld,MOLD=char_dummy,SIZE=payload_size)
    end do
    !$OMP end parallel do

    call mesh_dealloc(mesh)
    call p4wrap_partition(mesh%p4est)
    call mesh_alloc(mesh)

    timestamp = mesh%timestamp

    allocate(new_payloads(payload_size*mesh%nquads))
    call p4est_transfer_fixed(c_loc(mesh%gfq),c_loc(old_gfq),& 
        rt%mpi%ccomm,42,c_loc(new_payloads),c_loc(old_payloads),payload_size)

    !$OMP PARALLEL DO default(shared) private(i)
    do i = 1,mesh%nquads
        mesh%quads(i)%pld = transfer(&
            SOURCE=new_payloads(1+(i-1)*payload_size:i*payload_size),MOLD=payload_dummy)
    end do
    !$OMP end parallel do
 
    deallocate(old_payloads)
    deallocate(new_payloads)

# if DO_PROFILING
    call profiler_stop(prof_mesh_partition)
# endif

end subroutine

end module
