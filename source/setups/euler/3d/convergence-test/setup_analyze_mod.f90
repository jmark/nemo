# define ANGULAR_MOMENTUM 0

module setup_analyze_mod

use constants_mod

character(len=255) :: fpath = 'analysis.dat'

integer, parameter :: idens      = 1
integer, parameter :: imomx      = 2
integer, parameter :: imomy      = 3
integer, parameter :: imomz      = 4
integer, parameter :: ietot      = 5

integer, parameter :: idl2norm   = 7
integer, parameter :: idl8norm   = 8

integer, parameter :: iel2norm   = 9
integer, parameter :: iel8norm   = 10

integer, parameter :: sums_len  = 10

real(dp), save :: sums(1:SUMS_LEN) !! accumulators
real(dp), save :: nextanalyze = 0.0

real(dp), save :: nodes(N_NODES),ws(N_NODES)

contains

subroutine hook_init_830

    use interpol_mod, only: LegendreGaussNodesAndWeights

    call LegendreGaussNodesAndWeights(N_NODES-1,nodes,ws)
    nodes = 0.5_dp*nodes
    ws = 0.5_dp*ws

end subroutine

subroutine analyze_cb(quad)

    use setup_config_mod
    use equations_const_mod

    use equations_mod, only: pressure
    use mesh_types_mod, only: quad_t
    use quad_utils_mod, only: quad_get_volumes
    use quad_utils_mod, only: quad_get_coords
    use setup_manufacsol_mod
    use share_mod, only: rt => runtime

    type(quad_t), intent(inout) :: quad
    
    real(dp) :: coords(N_NODES,N_NODES,N_NODES,N_DIMS)

    real(dp), dimension(N_NODES,N_NODES,N_NODES) :: dens, momx, momy, momz, etot, ekin, eint, epot
    real(dp) :: vols(N_NODES,N_NODES,N_NODES)

    real(dp) :: u(N_NODES,N_NODES,N_NODES,N_VARS)

    real(dp) :: xs(N_NODES)
    real(dp) :: ys(N_NODES)
    real(dp) :: zs(N_NODES)

    real(dp) :: time

    integer :: i,j,k, ii,jj,kk

    vols = quad_get_volumes(quad)
    coords = quad_get_coords(quad)

    time = rt%simtime

    u = 0.0_dp
    do k = 1,N_NODES; do j = 1,N_NODES; do i = 1,N_NODES;
        xs = coords(i,j,k,1) + quad%delta(1)/REAL(N_NODES,dp)*nodes
        ys = coords(i,j,k,2) + quad%delta(2)/REAL(N_NODES,dp)*nodes
        zs = coords(i,j,k,3) + quad%delta(3)/REAL(N_NODES,dp)*nodes

        do kk = 1,N_NODES; do jj = 1,N_NODES; do ii = 1,N_NODES;
            u(i,j,k,:) = u(i,j,k,:) + ws(ii)*ws(jj)*ws(kk)*init_manufacsol(time,xs(ii),ys(jj),zs(kk))
        end do; end do; end do
    end do; end do; end do

    dens = quad%pld%hydro%state(:,:,:,DENS_VAR)
    momx = quad%pld%hydro%state(:,:,:,MOMX_VAR)
    momy = quad%pld%hydro%state(:,:,:,MOMY_VAR)
    momz = quad%pld%hydro%state(:,:,:,MOMZ_VAR)
    etot = quad%pld%hydro%state(:,:,:,ENER_VAR)

    !$omp critical
    sums(idens) = sums(idens) + sum(vols * dens)
    sums(imomx) = sums(imomx) + sum(vols * momx)
    sums(imomy) = sums(imomy) + sum(vols * momy)
    sums(imomz) = sums(imomz) + sum(vols * momz)
    sums(ietot) = sums(ietot) + sum(vols * etot)

    sums(idl2norm) = sums(idl2norm) + sum(vols*(u(:,:,:,DENS_VAR)-quad%pld%hydro%state(:,:,:,DENS_VAR))**2)
    sums(idl8norm) = max(sums(idl8norm),maxval(abs(u(:,:,:,DENS_VAR)-quad%pld%hydro%state(:,:,:,DENS_VAR))))

    sums(iel2norm) = sums(iel2norm) + sum(vols*(u(:,:,:,ENER_VAR)-quad%pld%hydro%state(:,:,:,ENER_VAR))**2)
    sums(iel8norm) = max(sums(iel8norm),maxval(abs(u(:,:,:,ENER_VAR)-quad%pld%hydro%state(:,:,:,ENER_VAR))))
    !$omp end critical

end subroutine

subroutine hook_analyze_500

    use mesh_utils_mod, only: mesh_iterate_quads
    use mpi_f08, only: mpi_reduce,MPI_SUM,MPI_DOUBLE_PRECISION

    use utils_mod, only: incr, write2file
    use share_mod, only: mesh, rt => runtime

    use setup_config_mod, only: dt_analyze

    integer :: ierr
    integer :: cnt = 0
    real(dp) :: resBuf(99)
    real(dp) :: recv(1:SUMS_LEN) !! mpi accumulator

    if (rt%simtime < nextanalyze) return

    nextanalyze = rt%simtime + dt_analyze

    sums = 0.0
 
    call mesh_iterate_quads(mesh,analyze_cb)

    call mpi_reduce(sums,recv,SUMS_LEN, MPI_DOUBLE_PRECISION, MPI_SUM, 0, rt%mpi%comm, ierr)

    if (rt%mpi%root) then
        cnt = 0

        resBuf(incr(cnt)) = rt%simtime
        resBuf(incr(cnt)) = rt%maxtimestep
        resBuf(incr(cnt)) = recv(idens)
        resBuf(incr(cnt)) = recv(ietot)

        resBuf(incr(cnt)) = sums(idl8norm) / REAL(rt%mpi%size,dp)
        resBuf(incr(cnt)) = SQRT(sums(idl2norm))

        resBuf(incr(cnt)) = sums(iel8norm) / REAL(rt%mpi%size,dp)
        resBuf(incr(cnt)) = SQRT(sums(iel2norm))

        call write2file(fpath, cnt, resBuf)
    end if

end subroutine

end module
