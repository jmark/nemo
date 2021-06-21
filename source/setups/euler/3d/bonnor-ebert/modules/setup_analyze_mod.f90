module setup_analyze_mod

use constants_mod

character(len=255) :: fpath = 'analysis.dat'

integer, parameter :: idens      = 1
integer, parameter :: imomx      = 2
integer, parameter :: imomy      = 3
integer, parameter :: imomz      = 4
integer, parameter :: iekin      = 5
integer, parameter :: iener      = 6
integer, parameter :: sums_len  = 10

real(dp), save :: sums(1:SUMS_LEN) !! accumulators
real(dp), save :: nextanalyze = 0.0

contains

subroutine analyze_cb(quad)

    use setup_config_mod
    use equations_const_mod

    use mesh_types_mod, only: quad_t
    use quad_utils_mod, only: quad_get_volumes

    type(quad_t), intent(inout) :: quad

    real(dp), dimension(N_NODES,N_NODES,N_NODES) :: dens, momx, momy, momz, ener, ekin
    real(dp) :: vols(N_NODES,N_NODES,N_NODES)

    vols = quad_get_volumes(quad)

    dens = quad%pld%hydro%state(:,:,:,DENS_VAR)
    momx = quad%pld%hydro%state(:,:,:,MOMX_VAR)
    momy = quad%pld%hydro%state(:,:,:,MOMY_VAR)
    momz = quad%pld%hydro%state(:,:,:,MOMZ_VAR)
    ener = quad%pld%hydro%state(:,:,:,ENER_VAR)
    ekin = 0.5 * (momx**2 + momy**2 + momz**2)/dens

    !$omp critical
    sums(idens) = sums(idens) + sum(vols * dens)
    sums(imomx) = sums(imomx) + sum(vols * momx)
    sums(imomy) = sums(imomy) + sum(vols * momy)
    sums(imomz) = sums(imomz) + sum(vols * momz)
    sums(iekin) = sums(iekin) + sum(vols * ekin)
    sums(iener) = sums(iener) + sum(vols * ener)
    !$omp end critical

end subroutine

subroutine hook_analyze_500

    use mesh_utils_mod, only: mesh_iterate_quads
    use mpi_f08, only: mpi_reduce,MPI_SUM,MPI_DOUBLE_PRECISION

    use utils_mod, only: incr, write2file
    use share_mod, only: mesh, rt => runtime

    use setup_config_mod, only: dt_analyze

    integer :: ierr
    integer :: count = 0
    real(dp) :: resBuf(10)
    real(dp) :: recv(1:SUMS_LEN) !! mpi accumulator

    if (rt%simtime < nextanalyze) return

    nextanalyze = rt%simtime + dt_analyze

    sums = 0.0
 
    call mesh_iterate_quads(mesh,analyze_cb)

    call mpi_reduce(sums,recv,SUMS_LEN, MPI_DOUBLE_PRECISION, MPI_SUM, 0, rt%mpi%comm, ierr)

    if (rt%mpi%root) then
        count = 0
        resBuf(incr(count)) = rt%simtime
        resBuf(incr(count)) = rt%maxtimestep
        resBuf(incr(count)) = recv(idens)
        resBuf(incr(count)) = recv(iener)
        call write2file(fpath, count, resBuf)
    end if

end subroutine

end module
