module setup_sourceterm_mod

use constants_mod

real(dp), parameter :: fmm3dprec = 1.0e-6

real(dp), allocatable, save :: charges(:)
real(dp), allocatable, save :: sources(:,:)
real(dp), allocatable, save :: targets(:,:)
real(dp), allocatable, save :: potential(:)
real(dp), allocatable, save :: gradients(:,:)
real(dp), allocatable, save :: coeffs(:,:,:)

integer(kind=1), allocatable, save :: alllevels(:)
integer(kind=8), allocatable, save :: recvdispls(:)
integer(kind=8), allocatable, save :: recvcounts(:)

integer(kind=8), save :: timestamp = 0
integer(kind=8), save :: nsources,ntargets

real(dp), save :: ccweights(8,8,8,8,8)
real(dp), save :: vcweights(8,8,N_NODES,N_NODES,N_NODES)

contains

subroutine hook_init_900

    use quad_utils_mod, only: nodes

    use cubic_interpol_mod
    use share_mod, only: mesh
    use mesh_const_mod, only: MESH_VERTS
    use mesh_utils_mod, only: mesh_iterate_quads

    integer :: i,j,k, iv,jv

    real(dp) :: xnodes(N_NODES), ynodes(N_NODES), znodes(N_NODES)
    real(dp) :: x,y,z

    do i = 1,N_NODES
        xnodes(i) = (0.5*(1.0_dp+nodes(i)))
    end do
    ynodes = xnodes
    znodes = xnodes

    forall (k = 1:N_NODES, j = 1:N_NODES, i = 1:N_NODES)
        vcweights(:,:,i,j,k) = tricubic_mat(xnodes(i),ynodes(j),znodes(k))
    end forall

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

    call hook_pre_timestep_605

end subroutine

subroutine hook_pre_timestep_605

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

    integer(kind=8) :: i

    if (mesh%timestamp > timestamp) then
        if (allocated(alllevels)) deallocate(alllevels)
        if (allocated(sources)) deallocate(sources)
        if (allocated(charges)) deallocate(charges)
        if (allocated(targets)) deallocate(targets)
        if (allocated(potential)) deallocate(potential)
        if (allocated(gradients)) deallocate(gradients)
        if (allocated(coeffs)) deallocate(coeffs)
        if (allocated(recvcounts)) deallocate(recvcounts)

        nsources = mesh%nglobal
        ntargets = mesh%nverts*(1+MESH_VERTS)

        allocate(alllevels(nsources))
        allocate(charges(  nsources))
        allocate(sources(3,nsources))
        allocate(targets(3,  ntargets))
        allocate(potential(  ntargets))
        allocate(gradients(3,ntargets))
        allocate(coeffs(MESH_VERTS,0:N_DIMS,mesh%nverts))

        allocate(recvcounts(size(mesh%gfq)-1))
        do i = 1,size(recvcounts)
            recvcounts(i) = mesh%gfq(i+1)-mesh%gfq(i)
        end do

        !$omp parallel do default(shared) private(i)
        do i = 1,mesh%nquads
            block
            integer(kind=8) :: iq
            iq = mesh%gfq(rt%mpi%rank+1)+i
            alllevels(iq) = mesh%quads(i)%level
            sources(:,iq) = mesh%quads(i)%center
            end block
        end do
        !$omp end parallel do

        if (rt%mpi%size > 1) then
            call MPI_Allgatherv(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,&
                alllevels,INT(recvcounts,kind=4),INT(mesh%gfq,kind=4),MPI_INTEGER1,rt%mpi%comm)

            call MPI_Allgatherv(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,&
                sources,INT(3*recvcounts,kind=4),3*INT(mesh%gfq,kind=4),MPI_DOUBLE_PRECISION,rt%mpi%comm)
        end if

        !$omp parallel do default(shared) private(i)
        do i = 1,mesh%nverts
            call calc_targets(i)
        end do
        !$omp end parallel do

        timestamp = mesh%timestamp
    end if

    !$omp parallel do default(shared) private(i)
    do i = 1,mesh%nquads
        associate(quad => mesh%quads(i))
        charges(mesh%gfq(rt%mpi%rank+1)+i) = -grav*SUM(quad_get_volumes(quad)*quad%pld%hydro%state(:,:,:,DENS_VAR))
        end associate
    end do
    !$omp end parallel do

    if (rt%mpi%size > 1) then
        call MPI_Allgatherv(MPI_IN_PLACE,0,MPI_DOUBLE_PRECISION,&
            charges,INT(recvcounts,kind=4),INT(mesh%gfq,kind=4),MPI_DOUBLE_PRECISION,rt%mpi%comm)
    end if

    call lfmm3d_t_c_g(fmm3dprec,INT(nsources,kind=4),sources,charges,INT(ntargets,kind=4),targets,potential,gradients,plummer)

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
    do i = 1,mesh%nverts
        call calc_coefficients(i)
    end do
    !$omp end parallel do

    !$omp parallel do default(shared) private(i)
    do i = 1,mesh%nquads
        call calc_potential(i)
    end do
    !$omp end parallel do

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

    dx = HUGE(1.0_dp)
    do iv = 1,MESH_VERTS
        if (mesh%verts(i)%iquads(iv) > 0) then
            dx = min(dx,minval(mesh%quads(mesh%verts(i)%iquads(iv))%delta))
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

    coeffs(:,0,i) = sd*finite_differences(potential(  j:j+MESH_VERTS))
    coeffs(:,1,i) = sd*finite_differences(gradients(1,j:j+MESH_VERTS))
    coeffs(:,2,i) = sd*finite_differences(gradients(2,j:j+MESH_VERTS))
    coeffs(:,3,i) = sd*finite_differences(gradients(3,j:j+MESH_VERTS))

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
    use quad_utils_mod, only: quad_get_coords
    use quad_utils_mod, only: quad_is_family
    use p8wrap_interfaces_mod, only: p8wrap_quadrant_child_id

    integer(kind=8), intent(in) :: i !! vertex id

    integer :: iv,jv,kv, ii,jj,kk

    real(dp) :: coords(N_NODES,N_NODES,N_NODES,N_DIMS)
    real(dp) :: qcoeffs(MESH_VERTS,MESH_VERTS,0:N_DIMS)
    real(dp) :: ccoeffs(MESH_VERTS,MESH_VERTS,0:N_DIMS)

    real(dp) :: delta(8),sd(8)
    real(dp) :: xelta(8),xd(8)

    integer :: child,ivA,ivB
    integer :: flips,idir,face

    integer(kind=8) :: iq,iquad

    associate(quad => mesh%quads(i))
    delta(1:N_DIMS) = quad%delta
    delta(4) = delta(1)*delta(2)
    delta(5) = delta(1)*delta(3)
    delta(6) = delta(2)*delta(3)
    delta(7) = product(delta(1:N_DIMS))

    sd = 1.0_dp/delta

    if (any(quad%iverts < 1)) then
        xelta(1:N_DIMS) = 2*quad%delta
        xelta(4) = xelta(1)*xelta(2)
        xelta(5) = xelta(1)*xelta(3)
        xelta(6) = xelta(2)*xelta(3)
        xelta(7) = product(xelta(1:N_DIMS))
        xd = 1.0_dp/xelta
    end if

    do iv = 1,MESH_VERTS
        if (quad%iverts(iv) > 0) then
            qcoeffs(:,iv,:) = coeffs(:,:,quad%iverts(iv))
        else
            child = p8wrap_quadrant_child_id(quad%level,quad%morton(1),quad%morton(2),quad%morton(3))+1
            flips = flipcount(child-1,iv-1,N_DIMS)

            !! hanging node along edge
            if (flips == 1) then
                idir = xor(child-1,iv-1)
                face = iand(iv-1,idir)

                do jv = 1,MESH_VERTS
                    iquad = mesh%verts(quad%iverts(child))%iquads(jv)
                    if (iquad < 1) cycle
                    if (mesh%quads(iquad)%level >= quad%level) cycle
                    if (iand(jv-1,idir) /= face) cycle

                    ivA = 1 +     iand(not(jv-1),2**N_DIMS-1)
                    ivB = 1 + xor(iand(not(jv-1),2**N_DIMS-1),idir)

                    ccoeffs = 0.0_dp
                    ccoeffs(:,ivA,:) = coeffs(:,:,mesh%quads(iquad)%iverts(ivA))
                    ccoeffs(:,ivB,:) = coeffs(:,:,mesh%quads(iquad)%iverts(ivB))

                    exit
                end do
        
            !! hanging node at a face center
            else

                ! idir = ~(corner ^ hanging) & ((1 << N_DIMS)-1)
                ! face =  (corner ^ idir) & idir
                idir = iand(not(xor(child-1,iv-1)),2**N_DIMS-1)
                face = iand(ior(child-1,idir),idir)

                iquad = mesh%verts(quad%iverts(child))%iquads(iv)

                ivA = xor(child-1,idir)+1
                ivB = xor(iv-1,idir)+1

                ccoeffs = 0.0_dp
                do jv = 1,MESH_VERTS
                    if (mesh%quads(iquad)%iverts(jv) < 1) cycle
                    ccoeffs(:,jv,:) = coeffs(:,:,mesh%quads(iquad)%iverts(jv))
                end do
            end if

            ccoeffs(2,:,:) = xelta(1)*ccoeffs(2,:,:)
            ccoeffs(3,:,:) = xelta(2)*ccoeffs(3,:,:)
            ccoeffs(4,:,:) = xelta(3)*ccoeffs(4,:,:)
            ccoeffs(5,:,:) = xelta(4)*ccoeffs(5,:,:)
            ccoeffs(6,:,:) = xelta(5)*ccoeffs(6,:,:)
            ccoeffs(7,:,:) = xelta(6)*ccoeffs(7,:,:)
            ccoeffs(8,:,:) = xelta(7)*ccoeffs(8,:,:)
        
            forall (ii = 0:N_DIMS)
                qcoeffs(1,iv,ii) =       sum(ccweights(:,:,1,ivA,ivB)*ccoeffs(:,:,ii))

                qcoeffs(2,iv,ii) = xd(1)*sum(ccweights(:,:,2,ivA,ivB)*ccoeffs(:,:,ii))
                qcoeffs(3,iv,ii) = xd(2)*sum(ccweights(:,:,3,ivA,ivB)*ccoeffs(:,:,ii))
                qcoeffs(4,iv,ii) = xd(3)*sum(ccweights(:,:,4,ivA,ivB)*ccoeffs(:,:,ii))

                qcoeffs(5,iv,ii) = xd(4)*sum(ccweights(:,:,5,ivA,ivB)*ccoeffs(:,:,ii))
                qcoeffs(6,iv,ii) = xd(5)*sum(ccweights(:,:,6,ivA,ivB)*ccoeffs(:,:,ii))
                qcoeffs(7,iv,ii) = xd(6)*sum(ccweights(:,:,7,ivA,ivB)*ccoeffs(:,:,ii))

                qcoeffs(8,iv,ii) = xd(7)*sum(ccweights(:,:,8,ivA,ivB)*ccoeffs(:,:,ii))
            end forall
        end if
    end do

    qcoeffs(2,:,:) = delta(1)*qcoeffs(2,:,:)
    qcoeffs(3,:,:) = delta(2)*qcoeffs(3,:,:)
    qcoeffs(4,:,:) = delta(3)*qcoeffs(4,:,:)

    qcoeffs(5,:,:) = delta(4)*qcoeffs(5,:,:)
    qcoeffs(6,:,:) = delta(5)*qcoeffs(6,:,:)
    qcoeffs(7,:,:) = delta(6)*qcoeffs(7,:,:)

    qcoeffs(8,:,:) = delta(7)*qcoeffs(8,:,:)

    forall (kk = 1:N_NODES, jj = 1:N_NODES, ii = 1:N_NODES)
        mesh%quads(i)%aux%gravity%pot(ii,jj,kk  ) =  sum(vcweights(:,:,ii,jj,kk)*qcoeffs(:,:,0))
        mesh%quads(i)%aux%gravity%acc(ii,jj,kk,1) = -sum(vcweights(:,:,ii,jj,kk)*qcoeffs(:,:,1))
        mesh%quads(i)%aux%gravity%acc(ii,jj,kk,2) = -sum(vcweights(:,:,ii,jj,kk)*qcoeffs(:,:,2))
        mesh%quads(i)%aux%gravity%acc(ii,jj,kk,3) = -sum(vcweights(:,:,ii,jj,kk)*qcoeffs(:,:,3))
    end forall
    end associate

end subroutine

subroutine setup_sourceterm(quad,rhs)

    use equations_const_mod
    use mesh_types_mod, only: quad_t

    type(quad_t), intent(inout) :: quad
    real(dp), intent(inout)     :: rhs(N_NODES,N_NODES,N_NODES,N_VARS)

    rhs = 0.0_dp

    associate(u => quad%pld%hydro%state, acc => quad%aux%gravity%acc)
    rhs(:,:,:,MOMX_VAR) = rhs(:,:,:,MOMX_VAR) + acc(:,:,:,1)*u(:,:,:,DENS_VAR)
    rhs(:,:,:,MOMY_VAR) = rhs(:,:,:,MOMY_VAR) + acc(:,:,:,2)*u(:,:,:,DENS_VAR)
    rhs(:,:,:,MOMZ_VAR) = rhs(:,:,:,MOMZ_VAR) + acc(:,:,:,3)*u(:,:,:,DENS_VAR)
    rhs(:,:,:,ENER_VAR) = rhs(:,:,:,ENER_VAR) + acc(:,:,:,1)*u(:,:,:,MOMX_VAR) &
                                              + acc(:,:,:,2)*u(:,:,:,MOMY_VAR) &
                                              + acc(:,:,:,3)*u(:,:,:,MOMZ_VAR)
    end associate

end subroutine

subroutine setup_sourceterm_dt(quad,dt)

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

subroutine modify_limiter_cb(quad)

    use mesh_types_mod, only: quad_t
    use setup_config_mod

    use quad_utils_mod, only: quad_get_coords
    use utils_mod, only: vmin

    type(quad_t), intent(inout) :: quad
    logical                     :: ok

    real(dp) :: dw

    real(dp) :: coords(N_NODES,N_NODES,N_DIMS)

    return

    !if (SQRT(SUM(quad%center**2)) < 1.2*cloudradius) then
    !    return 
    !end if

    if (quad%center(1) > xlim(1) + quad%delta(1) .and. quad%center(1) < xlim(2) - quad%delta(1)) then
        if (quad%center(2) > ylim(1) + quad%delta(2) .and. quad%center(2) < ylim(2) - quad%delta(2)) then
            if (quad%center(3) > zlim(1) + quad%delta(3) .and. quad%center(3) < zlim(2) - quad%delta(3)) return
        end if
    end if

    quad%aux%hydro%subcell = .true.

    !quad%aux%hydro%qmin = 1
    !quad%aux%hydro%qmax = 0

    !quad%aux%hydro%blend = 0.0_dp
    quad%aux%hydro%blend = vmin(0.5_dp,quad%aux%hydro%blend)

end subroutine

!! subroutine setup_sourceterm_eos_cb(quad)
!! 
!!     use mesh_types_mod, only: quad_t
!!     use equations_const_mod
!!     use setup_config_mod
!! 
!!     type(quad_t), intent(inout) :: quad
!! 
!!     associate(u => quad%pld%hydro%state)
!!     !where (u(:,:,DENS_VAR) < dens_min)
!!     !    u(:,:,DENS_VAR) = dens_min
!!     !end where
!! 
!!     where (u(:,:,DENS_VAR) < dens_critical)
!!         u(:,:,ENER_VAR) = 0.5/u(:,:,DENS_VAR)*(u(:,:,MOMX_VAR)**2 + u(:,:,MOMY_VAR)**2) + K_isothermal*u(:,:,DENS_VAR)/(kappa-1)
!!     else where
!!         u(:,:,ENER_VAR) = 0.5/u(:,:,DENS_VAR)*(u(:,:,MOMX_VAR)**2 + u(:,:,MOMY_VAR)**2) + K_adiabatic*u(:,:,DENS_VAR)**kappa/(kappa-1)
!!     end where
!!     !u(:,:,ENER_VAR) = 0.5/u(:,:,DENS_VAR)*(u(:,:,MOMX_VAR)**2 + u(:,:,MOMY_VAR)**2) + K_isothermal*u(:,:,DENS_VAR)/(kappa-1)
!!     end associate
!! 
!! end subroutine
!! 
!! subroutine hook_post_rk_timestep_200
!! 
!!     use share_mod, only: mesh
!!     use mesh_utils_mod, only: mesh_iterate_quads
!! 
!!     call mesh_iterate_quads(mesh,setup_sourceterm_eos_cb)
!! 
!! end subroutine

end module
