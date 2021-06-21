# define ANGULAR_MOMENTUM 0

module setup_analyze_mod

use constants_mod

character(len=255) :: fpath = 'analysis.dat'

integer, parameter :: idens      = 1
integer, parameter :: imomx      = 2
integer, parameter :: imomy      = 3
integer, parameter :: ietot      = 4

integer, parameter :: idl2norm   = 7
integer, parameter :: idl8norm   = 8

integer, parameter :: iel2norm   = 9
integer, parameter :: iel8norm   = 10

integer, parameter :: sums_len  = 10

real(dp), save :: sums(1:SUMS_LEN) !! accumulators
real(dp), save :: nextanalyze = 0.0_dp

contains

subroutine analyze_cb(quad)

    use setup_config_mod
    use equations_const_mod

    use equations_mod, only: pressure
    use mesh_types_mod, only: quad_t
    use quad_utils_mod, only: quad_get_volumes
    use quad_utils_mod, only: quad_get_coords

    type(quad_t), intent(inout) :: quad
    
    real(dp) :: coords(N_NODES,N_NODES,N_DIMS)

    real(dp), dimension(N_NODES,N_NODES) :: dens, momx, momy, momz, etot, ekin, eint, epot
    real(dp) :: vols(N_NODES,N_NODES)

    vols = quad_get_volumes(quad)
    coords = quad_get_coords(quad)

    dens = quad%pld%hydro%state(:,:,DENS_VAR)
    momx = quad%pld%hydro%state(:,:,MOMX_VAR)
    momy = quad%pld%hydro%state(:,:,MOMY_VAR)
    etot = quad%pld%hydro%state(:,:,ENER_VAR)

    !$omp critical
    sums(idens) = sums(idens) + sum(vols * dens)
    sums(imomx) = sums(imomx) + sum(vols * momx)
    sums(imomy) = sums(imomy) + sum(vols * momy)
    sums(ietot) = sums(ietot) + sum(vols * etot)

    sums(idl2norm) = sums(idl2norm) + sum(vols*(quad%aux%manufacsol%u0(:,:,DENS_VAR)-quad%pld%hydro%state(:,:,DENS_VAR))**2)
    sums(idl8norm) = max(sums(idl8norm),maxval(abs(quad%aux%manufacsol%u0(:,:,DENS_VAR)-quad%pld%hydro%state(:,:,DENS_VAR))))

    sums(iel2norm) = sums(iel2norm) + sum(vols*(quad%aux%manufacsol%u0(:,:,ENER_VAR)-quad%pld%hydro%state(:,:,ENER_VAR))**2)
    sums(iel8norm) = max(sums(iel8norm),maxval(abs(quad%aux%manufacsol%u0(:,:,ENER_VAR)-quad%pld%hydro%state(:,:,ENER_VAR))))
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
