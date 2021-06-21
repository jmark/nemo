!# define MODULE_SELF_GRAVITY 1
!# define MODULE_SOURCETERM 1

module gravity_fmm_mod

use constants_mod
use mesh_const_mod, only: MESH_VERTS

real(dp), parameter :: fmm3dprec = 1.0e-6

real(dp), allocatable, save :: charges(:)
real(dp), allocatable, save :: sources(:,:)
real(dp), allocatable, save :: targets(:,:)
real(dp), allocatable, save :: potential(:)
real(dp), allocatable, save :: gradients(:,:)
real(dp), allocatable, save :: coefficients(:,:,:)

integer, allocatable, save :: recvcounts_1(:)
integer, allocatable, save :: recvcounts_3(:)

integer, allocatable, save :: recvdispls_1(:)
integer, allocatable, save :: recvdispls_3(:)

integer(kind=8), save :: timestamp = 0
integer(kind=8), save :: nsources,ntargets

real(dp), save :: ccweights(8,8,8,8,8)
real(dp), save :: vcweights(8,8,N_NODES,N_NODES,N_NODES)

# if DO_PROFILING
    integer :: prof_self_gravity
# endif

contains

subroutine hook_init_900

    use cubic_interpol_mod
    use share_mod, only: mesh
    use mesh_const_mod, only: MESH_VERTS
    use mesh_utils_mod, only: mesh_iterate_quads

    use interpol_mod, only: LegendreGaussNodesAndWeights
    use interpol_mod, only: LegendreGaussLobattoNodesAndWeights

    use quad_utils_mod, only: nodes

# if DO_PROFILING
    use profiler_mod, only: profiler_init
# endif

    integer :: i,j,k, iv,jv

    real(dp) :: xnodes(N_NODES), ynodes(N_NODES), znodes(N_NODES)
    real(dp) :: weights(N_NODES)
    real(dp) :: x,y,z

    !call LegendreGaussNodesAndWeights(N_NODES-1, xnodes, weights)
    xnodes = nodes

    do i = 1,N_NODES
        xnodes(i) = (0.5*(1.0_dp+xnodes(i)))
    end do
    ynodes = xnodes
    znodes = xnodes

    do k = 1,N_NODES; do j = 1,N_NODES; do i = 1,N_NODES;
        vcweights(:,:,i,j,k) = tricubic_mat(xnodes(i),ynodes(j),znodes(k))
    end do; end do; end do;

    do iv = 1,MESH_VERTS
        do jv = 1,MESH_VERTS
            x = 0.5_dp*(iand(ishft(iv-1, 0),1) + iand(ishft(jv-1, 0),1))
            y = 0.5_dp*(iand(ishft(iv-1,-1),1) + iand(ishft(jv-1,-1),1))
            z = 0.5_dp*(iand(ishft(iv-1,-2),1) + iand(ishft(jv-1,-2),1))

            ccweights(:,:,1,iv,jv) = tricubic_mat      (x,y,z)
            ccweights(:,:,2,iv,jv) = tricubic_x_mat    (x,y,z)
            ccweights(:,:,3,iv,jv) = tricubic_y_mat    (x,y,z)
            ccweights(:,:,4,iv,jv) = tricubic_z_mat    (x,y,z)
            ccweights(:,:,5,iv,jv) = tricubic_xy_mat   (x,y,z)
            ccweights(:,:,6,iv,jv) = tricubic_xz_mat   (x,y,z)
            ccweights(:,:,7,iv,jv) = tricubic_yz_mat   (x,y,z)
            ccweights(:,:,8,iv,jv) = tricubic_xyz_mat  (x,y,z)
        end do
    end do

    call calc_self_gravity

# if DO_PROFILING
    prof_self_gravity = profiler_init('self_gravity')
# endif

end subroutine

subroutine hook_pre_timestep_605

    call calc_self_gravity

end subroutine

subroutine calc_self_gravity

    use share_mod, only: rt => runtime

    use mpi_f08, only: MPI_BYTE
    use mpi_f08, only: MPI_IN_PLACE
    use mpi_f08, only: MPI_INTEGER1
    use mpi_f08, only: MPI_Allgatherv
    use mpi_f08, only: MPI_DATATYPE_NULL
    use mpi_f08, only: MPI_DOUBLE_PRECISION

    use share_mod, only: mesh
    use FMM3D_mod, only: lfmm3d_t_c_g
    use setup_config_mod, only: plummer, diameter

    use quad_utils_mod, only: quad_get_volumes
    use quad_utils_mod, only: vert_get_coords

    use mesh_const_mod, only: MESH_CORNERS, MESH_VERTS, MESH_FACES, MESH_CHILDREN
    use mesh_const_mod, only: map_faceid_to_hanging_vertid_quads
    use mesh_const_mod, only: map_faceid_to_hanging_vertid_verts
    use mesh_const_mod, only: map_to_opposite_vertid
    use mesh_const_mod, only: map_faceid_to_opposite
    use mesh_const_mod, only: map_faceid_to_vertids

    use setup_config_mod, only: grav,plummer
    use equations_const_mod, only: DENS_VAR

# if DO_PROFILING
    use profiler_mod, only: profiler_start,profiler_stop
# endif

    integer(kind=8) :: i

# if DO_PROFILING
    call profiler_start(prof_self_gravity)
# endif

    if (mesh%timestamp > timestamp) then
        if (allocated(sources)) deallocate(sources)
        if (allocated(charges)) deallocate(charges)
        if (allocated(targets)) deallocate(targets)
        if (allocated(potential)) deallocate(potential)
        if (allocated(gradients)) deallocate(gradients)
        if (allocated(coefficients)) deallocate(coefficients)

        if (allocated(recvcounts_1)) deallocate(recvcounts_1)
        if (allocated(recvcounts_3)) deallocate(recvcounts_3)

        if (allocated(recvdispls_1)) deallocate(recvdispls_1)
        if (allocated(recvdispls_3)) deallocate(recvdispls_3)

        nsources = mesh%nglobal
        ntargets = mesh%nverts_with_ghosts*(1+MESH_VERTS)

        allocate(charges(  nsources))
        allocate(sources(3,nsources))

        allocate(targets(3,  ntargets))
        allocate(potential(  ntargets))
        allocate(gradients(3,ntargets))
        allocate(coefficients(MESH_VERTS,0:N_DIMS,mesh%nverts_with_ghosts))

        allocate(recvcounts_1(size(mesh%gfq)-1))
        allocate(recvcounts_3(size(mesh%gfq)-1))
        do i = 1,size(recvcounts_1)
            recvcounts_1(i) = INT(mesh%gfq(i+1)-mesh%gfq(i),kind=4)
        end do

        recvcounts_3 = 3*recvcounts_1

        allocate(recvdispls_1(size(mesh%gfq)))
        allocate(recvdispls_3(size(mesh%gfq)))

        recvdispls_1 = INT(mesh%gfq,kind=4)
        recvdispls_3 = 3*recvdispls_1

        !$omp parallel do default(shared) private(i)
        do i = 1,mesh%nquads
            sources(:,mesh%gfq(rt%mpi%rank+1)+i) = mesh%quads(i)%center
        end do
        !$omp end parallel do

        if (rt%mpi%size > 1) then
            call MPI_Allgatherv(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,&
                sources,recvcounts_3,recvdispls_3,MPI_DOUBLE_PRECISION,rt%mpi%comm)
        end if

        !$omp parallel do default(shared) private(i)
        do i = 1,mesh%nverts_with_ghosts
            call calc_targets(i)
        end do
        !$omp end parallel do

        timestamp = mesh%timestamp
    end if

    !$omp parallel do default(shared) private(i)
    do i = 1,mesh%nquads
        charges(mesh%gfq(rt%mpi%rank+1)+i) = &
            -grav*SUM(quad_get_volumes(mesh%quads(i))*mesh%quads(i)%pld%hydro%state(:,:,:,DENS_VAR))
    end do
    !$omp end parallel do

    if (rt%mpi%size > 1) then
        call MPI_Allgatherv(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,&
            charges,recvcounts_1,recvdispls_1,MPI_DOUBLE_PRECISION,rt%mpi%comm)
    end if

    call lfmm3d_t_c_g(fmm3dprec,INT(nsources,kind=4),sources,charges,INT(ntargets,kind=4),targets,potential,gradients,plummer)

    if (any(isnan(gradients))) then
        do i = 1,size(targets,dim=2)
            if (any(isnan(gradients(:,i)))) then
                write (*,'(a12,i7,99(ES12.4))') 'targets', i, targets(:,i),potential(i),gradients(:,i)
            end if
        end do
    end if

# if 0
    !! debugging
    do i = 1,ntargets
        associate(quad => mesh%quads(i))
        !potential(  i) = targets(3,i)
        !gradients(:,i) = targets(3,i)
        potential(  i) = targets(1,i) * targets(2,i) * targets(3,i)
        gradients(:,i) = targets(1,i) * targets(2,i) * targets(3,i)
        end associate
    end do
# endif

    !$omp parallel do default(shared) private(i)
    do i = 1,mesh%nverts_with_ghosts
        call calc_coefficients(i)
    end do
    !$omp end parallel do

    !$omp parallel do default(shared) private(i)
    do i = 1,mesh%nquads
        call calc_potential(i)
    end do
    !$omp end parallel do

# if DO_PROFILING
    call profiler_stop(prof_self_gravity)
# endif

end subroutine

subroutine calc_targets(i)

    use share_mod, only: mesh
    use mesh_const_mod, only: MESH_VERTS
    use quad_utils_mod, only: vert_get_coords

    integer(kind=8), intent(in) :: i !! vertex id

    integer(kind=8) :: j
    integer :: iv

    real(dp) :: coords(N_DIMS),delta(N_DIMS),dx
    real(dp) :: x,y,z

    coords = vert_get_coords(mesh%verts(i))

    !dx = HUGE(1.0_dp)
    dx = -1.0_dp
    do iv = 1,MESH_VERTS
        if (mesh%verts(i)%iquads(iv) > 0) then
            !dx = min(dx,minval(mesh%quads(mesh%verts(i)%iquads(iv))%delta))
            dx = max(dx,maxval(mesh%quads(mesh%verts(i)%iquads(iv))%delta))
        end if
    end do

    j = 1 + (i-1)*(1+MESH_VERTS)
    targets(:,j) = coords

    do iv = 1,MESH_VERTS
        x = 2.0_dp*(iand(ishft(iv-1, 0),1)-0.5_dp)
        y = 2.0_dp*(iand(ishft(iv-1,-1),1)-0.5_dp)
        z = 2.0_dp*(iand(ishft(iv-1,-2),1)-0.5_dp)

        targets(:,j+iv) = coords + dx*(/x,y,z/)
    end do

end subroutine

pure function finite_differences(fs) result(qc)
        
    real(dp), intent(in) :: fs(9)
    real(dp)             :: qc(8)

    !! function value
    qc(1) = fs(1)

    !! dx
    qc(2) = 0.25*(fs(1+SWF) - fs(1+NWF) + fs(1+SEF) - fs(1+NEF) &
                + fs(1+SWB) - fs(1+NWB) + fs(1+SEB) - fs(1+NEB))

    !! dy
    qc(3) = 0.25*(fs(1+NEF) - fs(1+NWF) + fs(1+SEF) - fs(1+SWF) &
                + fs(1+NEB) - fs(1+NWB) + fs(1+SEB) - fs(1+SWB))

    !! dz
    qc(4) = 0.25*(fs(1+NWB) - fs(1+NWF) + fs(1+SWB) - fs(1+SWF) &
                + fs(1+NEB) - fs(1+NEF) + fs(1+SEB) - fs(1+SEF))

    !! dxdy
    qc(5) = 0.5*((fs(1+SEF) - fs(1+NEF) + fs(1+SEB) - fs(1+NEB)) &
               - (fs(1+SWF) - fs(1+NWF) + fs(1+SWB) - fs(1+NWB)))

    !! dxdz
    qc(6) = 0.5*((fs(1+SWB) - fs(1+NWB) + fs(1+SEB) - fs(1+NEB)) &
               - (fs(1+SWF) - fs(1+NWF) + fs(1+SEF) - fs(1+NEF)))

    !! dydz
    qc(7) = 0.5*((fs(1+NEB) - fs(1+NWB) + fs(1+SEB) - fs(1+SWB)) &
               - (fs(1+NEF) - fs(1+NWF) + fs(1+SEF) - fs(1+SWF)))

    !! dxdydz
    qc(8) = (((fs(1+SEB) - fs(1+NEB)) - (fs(1+SWB) - fs(1+NWB))) &
          - ((fs(1+SEF) - fs(1+NEF)) - (fs(1+SWF) - fs(1+NWF))))

end function

subroutine calc_coefficients(i)

    use mesh_const_mod, only: MESH_VERTS

    integer(kind=8), intent(in) :: i !! vertex id

    integer(kind=8) :: j
    real(dp)        :: delta(N_DIMS),sd(8)

    j = 1 + (i-1)*(1+MESH_VERTS)

    delta = targets(:,j+SEB)-targets(:,j+NWF)

    sd(1) = 1.0_dp
    sd(2) = 1.0_dp/delta(1)
    sd(3) = 1.0_dp/delta(2)
    sd(4) = 1.0_dp/delta(3)
    sd(5) = 1.0_dp/delta(1)/delta(2)
    sd(6) = 1.0_dp/delta(1)/delta(3)
    sd(7) = 1.0_dp/delta(2)/delta(3)
    sd(8) = 1.0_dp/product(delta)

    coefficients(:,0,i) = sd*finite_differences(potential(  j:j+MESH_VERTS))
    coefficients(:,1,i) = sd*finite_differences(gradients(1,j:j+MESH_VERTS))
    coefficients(:,2,i) = sd*finite_differences(gradients(2,j:j+MESH_VERTS))
    coefficients(:,3,i) = sd*finite_differences(gradients(3,j:j+MESH_VERTS))

    block
    integer :: ii
    if (any(isnan(coefficients(:,:,i)))) then
        write (*,'(a12,1(i7),99(ES12.4))') 'coefficients', i, delta
        do ii = 0,N_DIMS
            write (*,'(99(ES12.4))') coefficients(:,ii,i)
        end do
        write (*,*)
    end if
    end block

end subroutine

pure function flipcount(a,b,nbits) result(c)

    integer, intent(in) :: a,b,nbits
    integer             :: n,i
    integer             :: c

    c = 0
    n = 1
    
    do i = 1,nbits
        if (iand(a,n) /= iand(b,n)) c = c + 1
        n = 2*n
    end do

end function

subroutine calc_potential(i)

    use share_mod, only: mesh
    use mesh_const_mod, only: MESH_CHILDREN,MESH_VERTS
    use quad_utils_mod, only: quad_is_family
    use p8wrap_interfaces_mod, only: p8wrap_quadrant_child_id

    use tensor_mod, only: blockmulXYZ
    use kernel_const_mod, only: DGFV_QMAX
    use kernel_share_mod, only: dg2fvMat

    integer(kind=8), intent(in) :: i !! vertex id

    integer :: iv,jv,kv, ii,jj,kk

    real(dp) :: coords(N_NODES,N_NODES,N_NODES,N_DIMS)
    real(dp) :: qcoeffs(MESH_VERTS,MESH_VERTS,0:N_DIMS)
    real(dp) :: ccoeffs(MESH_VERTS,MESH_VERTS,0:N_DIMS)

    real(dp) :: delta(MESH_VERTS),sd(MESH_VERTS)
    real(dp) :: xelta(MESH_VERTS),xd(MESH_VERTS)

    integer :: child,ivA,ivB
    integer :: flips,idir,face

    integer(kind=8) :: iq,iquad

    associate(quad => mesh%quads(i))
    delta(1) = 1.0_dp
    delta(2:4) = quad%delta
    delta(5) = delta(2)*delta(3)
    delta(6) = delta(2)*delta(4)
    delta(7) = delta(3)*delta(4)
    delta(8) = product(delta(2:4))
    sd       = 1.0_dp/delta

    if (any(quad%iverts < 1)) then
        xelta(1) = 1.0_dp
        xelta(2:4) = 2.0_dp*quad%delta
        xelta(5) = xelta(2)*xelta(3)
        xelta(6) = xelta(2)*xelta(4)
        xelta(7) = xelta(3)*xelta(4)
        xelta(8) = product(xelta(2:4))
        xd       = 1.0_dp/xelta
    end if

    do iv = 1,MESH_VERTS
        if (quad%iverts(iv) > 0) then
            do ii = 0,N_DIMS;
                qcoeffs(:,iv,ii) = coefficients(:,ii,quad%iverts(iv))
            end do

            if (any(isnan(qcoeffs(:,iv,:)))) then
                write (*,'(a12,99(i7))') 'conforming', quad%iverts, iv
                do ii = 0,N_DIMS
                    write (*,'(99(ES12.4))') qcoeffs(:,iv,ii)
                end do
                write (*,*)
            end if
        else
            child = p8wrap_quadrant_child_id(quad%level,quad%morton(1),quad%morton(2),quad%morton(3))+1
            flips = flipcount(child-1,iv-1,N_DIMS)

            if (flips == 1) then !! hanging node along edge
                idir = xor(child-1,iv-1)
                face = iand(iv-1,idir)

                do jv = 1,MESH_VERTS
                    iquad = mesh%verts(quad%iverts(child))%iquads(jv)
                    if (iquad < 1) cycle
                    if (mesh%quads(iquad)%level >= quad%level) cycle
                    if (iand(jv-1,idir) /= face) cycle

                    ivA = 1 +     iand(not(jv-1),2**N_DIMS-1)
                    ivB = 1 + xor(iand(not(jv-1),2**N_DIMS-1),idir)

                    if (mesh%quads(iquad)%iverts(ivA) > 0 .and. mesh%quads(iquad)%iverts(ivB) > 0) exit
                end do

                ccoeffs = 0.0_dp
                do ii = 0,N_DIMS
                    ccoeffs(:,ivA,ii) = coefficients(:,ii,mesh%quads(iquad)%iverts(ivA))
                    ccoeffs(:,ivB,ii) = coefficients(:,ii,mesh%quads(iquad)%iverts(ivB))
                end do

                if (any(isnan(ccoeffs)) .or. ivA < 1 .or. ivB < 1 .or. ivA > MESH_VERTS .or. ivB > MESH_VERTS) then
                    write (*,'(a18,99(i7))')    'hanging vert 1', mesh%quads(iquad)%iverts, jv, ivA, ivB
                    write (*,'(a18,99(i7))')    'hanging vert 2', mesh%quads(iquad)%iverts(ivA), mesh%quads(iquad)%iverts(ivB)

                    do ii = 0,N_DIMS
                        write (*,'(a18,99(ES12.4))')'ivA           ', coefficients(:,ii,mesh%quads(iquad)%iverts(ivA))
                    end do
                    write (*,*)
                    do ii = 0,N_DIMS
                        write (*,'(a18,99(ES12.4))')'ivB           ', coefficients(:,ii,mesh%quads(iquad)%iverts(ivB))
                    end do
                    write (*,*)

                    do ii = 0,N_DIMS
                    do jv = 1,MESH_VERTS
                        write (*,'(99(ES12.4))') ccoeffs(:,jv,ii)
                    end do
                    write (*,*)
                    end do
                    write (*,*)
                end if

            else !! hanging node at a face center
                ! idir = ~(corner ^ hanging) & ((1 << N_DIMS)-1)
                ! face =  (corner ^ idir) & idir
                idir = iand(not(xor(child-1,iv-1)),2**N_DIMS-1)
                face = iand(ior(child-1,idir),idir)

                iquad = mesh%verts(quad%iverts(child))%iquads(iv)

                ivA = xor(child-1,idir)+1
                ivB = xor(iv-1,idir)+1

                ccoeffs = 0.0_dp
                do jv = 1,MESH_VERTS
                    if (mesh%quads(iquad)%iverts(jv) > 0) then
                        do ii = 0,N_DIMS
                            ccoeffs(:,jv,ii) = coefficients(:,ii,mesh%quads(iquad)%iverts(jv))
                        end do
                    end if
                end do

                if (any(isnan(ccoeffs))) then
                    write (*,'(a12,99(i7))') 'hanging face', quad%iverts, iv
                    do ii = 0,N_DIMS
                    do jv = 1,MESH_VERTS
                        write (*,'(99(ES12.4))') ccoeffs(:,jv,ii)
                    end do
                    write (*,*)
                    end do
                    write (*,*)
                end if

            end if

            do ii = 0,N_DIMS; do kv = 1,MESH_VERTS; do jv = 1,MESH_VERTS;
                ccoeffs(jv,kv,ii) = xelta(jv)*ccoeffs(jv,kv,ii)
            end do; end do; end do;
       
            do ii = 0,N_DIMS; do jv = 1,MESH_VERTS;
                qcoeffs(jv,iv,ii) = xd(jv)*sum(ccweights(:,:,jv,ivA,ivB)*ccoeffs(:,:,ii))
            end do; end do;
        end if
    end do

    do ii = 0,N_DIMS; do iv = 1,MESH_VERTS; do jv = 1,MESH_VERTS;
        qcoeffs(jv,iv,ii) = delta(jv)*qcoeffs(jv,iv,ii)
    end do; end do; end do;

    ! qcoeffs(2,:,:) = delta(1)*qcoeffs(2,:,:)
    ! qcoeffs(3,:,:) = delta(2)*qcoeffs(3,:,:)
    ! qcoeffs(4,:,:) = delta(3)*qcoeffs(4,:,:)

    ! qcoeffs(5,:,:) = delta(4)*qcoeffs(5,:,:)
    ! qcoeffs(6,:,:) = delta(5)*qcoeffs(6,:,:)
    ! qcoeffs(7,:,:) = delta(6)*qcoeffs(7,:,:)

    ! qcoeffs(8,:,:) = delta(7)*qcoeffs(8,:,:)

    do kk = 1,N_NODES; do jj = 1,N_NODES; do ii = 1,N_NODES;
        quad%aux%gravity%pot(ii,jj,kk  ) =  sum(vcweights(:,:,ii,jj,kk)*qcoeffs(:,:,0))
        quad%aux%gravity%acc(ii,jj,kk,1) = -sum(vcweights(:,:,ii,jj,kk)*qcoeffs(:,:,1))
        quad%aux%gravity%acc(ii,jj,kk,2) = -sum(vcweights(:,:,ii,jj,kk)*qcoeffs(:,:,2))
        quad%aux%gravity%acc(ii,jj,kk,3) = -sum(vcweights(:,:,ii,jj,kk)*qcoeffs(:,:,3))
    end do; end do; end do;

    !if (any(abs(quad%aux%gravity%pot*0.0_dp) /= 0.0_dp)) then
    if (any(isnan(quad%aux%gravity%pot))) then
        write (*,'(a12,99(i7))')     'NaN gpot', 0, i, quad%iverts
        write (*,'(99(ES12.4))') sum(quad%aux%gravity%pot)
        write (*,'(99(ES12.4))') ccoeffs(:,:,0)
        write (*,'(99(ES12.4))') qcoeffs(:,:,0)
        write (*,*)
    end if

    do ii = 1,N_DIMS
        if (any(isnan(quad%aux%gravity%acc(:,:,:,ii)))) then
            write (*,'(a12,99(i7))') 'NaN acc', ii, i, quad%iverts
            ! write (*,'(a12,99(ES12.4))') 'accs   ', sum(quad%aux%gravity%acc(:,:,:,ii))
            do iv = 1,MESH_VERTS
                write (*,'(a12,99(ES12.4))') 'ccoeffs', ccoeffs(:,iv,ii)
            end do
            write (*,*)
            do iv = 1,MESH_VERTS
                write (*,'(a12,99(ES12.4))') 'qcoeffs', qcoeffs(:,iv,ii)
            end do
            write (*,*)
        end if
    end do

    ! quad%aux%gravity%pot          = blockmulXYZ(dg2fvMat(:,:,DGFV_QMAX),quad%aux%gravity%pot)
    ! quad%aux%gravity%acc(:,:,:,1) = blockmulXYZ(dg2fvMat(:,:,DGFV_QMAX),quad%aux%gravity%acc(:,:,:,1))
    ! quad%aux%gravity%acc(:,:,:,2) = blockmulXYZ(dg2fvMat(:,:,DGFV_QMAX),quad%aux%gravity%acc(:,:,:,2))
    ! quad%aux%gravity%acc(:,:,:,3) = blockmulXYZ(dg2fvMat(:,:,DGFV_QMAX),quad%aux%gravity%acc(:,:,:,3))
    end associate

end subroutine

pure subroutine calc_sourceterm(quad,rhs)

    use equations_const_mod
    use mesh_types_mod, only: quad_t

    type(quad_t), intent(inout) :: quad
    real(dp), intent(out)       :: rhs(N_NODES,N_NODES,N_NODES,N_VARS)

    associate(u => quad%pld%hydro%state, acc => quad%aux%gravity%acc)
    rhs(:,:,:,DENS_VAR) = 0.0_dp
    rhs(:,:,:,MOMX_VAR) = acc(:,:,:,X_DIR)*u(:,:,:,DENS_VAR)
    rhs(:,:,:,MOMY_VAR) = acc(:,:,:,Y_DIR)*u(:,:,:,DENS_VAR)
    rhs(:,:,:,MOMZ_VAR) = acc(:,:,:,Z_DIR)*u(:,:,:,DENS_VAR)
    rhs(:,:,:,ENER_VAR) = acc(:,:,:,X_DIR)*u(:,:,:,MOMX_VAR) &
                        + acc(:,:,:,Y_DIR)*u(:,:,:,MOMY_VAR) &
                        + acc(:,:,:,Z_DIR)*u(:,:,:,MOMZ_VAR)
    end associate

end subroutine

subroutine calc_timestep(quad,dt)

    use mesh_types_mod, only: quad_t
    use equations_mod, only: soundspeed
    use setup_config_mod, only: kappa

    type(quad_t), intent(in)    :: quad
    real(dp), intent(inout)     :: dt

    real(dp) :: snds(N_NODES,N_NODES,N_NODES)

    snds = soundspeed(quad%pld%hydro%state)

    associate(u => quad%pld%hydro%state, acc => quad%aux%gravity%acc)
    dt = minval(snds/SQRT(2*kappa*(kappa-1))/SQRT(acc(:,:,:,1)**2+acc(:,:,:,2)**2+acc(:,:,:,3)))
    end associate

end subroutine

end module
