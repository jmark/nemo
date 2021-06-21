# define ANGULAR_MOMENTUM 0

module setup_analyze_mod

use constants_mod

character(len=255) :: fpath = 'analysis.dat'

integer, parameter :: idens      = 1
integer, parameter :: imomx      = 2
integer, parameter :: imomy      = 3
integer, parameter :: imomz      = 4
integer, parameter :: ieint      = 5
integer, parameter :: iekin      = 6
integer, parameter :: iepot      = 7
integer, parameter :: ietot      = 8

integer, parameter :: iLtot      = 9    !! angular momentum
integer, parameter :: iHtot      = 10   !! specific angular momentum

integer, parameter :: sums_len  = 10

real(dp), save :: sums(1:SUMS_LEN) !! accumulators
real(dp), save :: nextanalyze = 0.0

contains

subroutine analyze_cb(quad)

    use setup_config_mod
    use equations_const_mod

    use equations_mod, only: pressure
    use mesh_types_mod, only: quad_t
    use quad_utils_mod, only: quad_get_volumes
    use quad_utils_mod, only: quad_get_coords

    type(quad_t), intent(inout) :: quad
    
    real(dp) :: coords(N_NODES,N_NODES,N_NODES,N_DIMS)

    real(dp), dimension(N_NODES,N_NODES,N_NODES) :: dens, momx, momy, momz, etot, ekin, eint, epot
    real(dp) :: vols(N_NODES,N_NODES,N_NODES)

# if ANGULAR_MOMENTUM
    real(dp) :: angx(N_NODES,N_NODES,N_NODES) !! x-angular momentum
    real(dp) :: angy(N_NODES,N_NODES,N_NODES) !! y-angular momentum
    real(dp) :: angz(N_NODES,N_NODES,N_NODES) !! z-angular momentum

    real(dp) :: angL(N_NODES,N_NODES,N_NODES) !! angular momentum
    real(dp) :: angH(N_NODES,N_NODES,N_NODES) !! specific angular momentum
# endif

    vols = quad_get_volumes(quad)
    coords = quad_get_coords(quad)

    dens = quad%pld%hydro%state(:,:,:,DENS_VAR)
    momx = quad%pld%hydro%state(:,:,:,MOMX_VAR)
    momy = quad%pld%hydro%state(:,:,:,MOMY_VAR)
    momz = quad%pld%hydro%state(:,:,:,MOMZ_VAR)
    etot = quad%pld%hydro%state(:,:,:,ENER_VAR)

    eint = pressure(quad%pld%hydro%state)/(kappa-1)
    ekin = 0.5 * (momx**2 + momy**2 + momz**2)/dens
# if MODULE_SELF_GRAVITY
    epot = 0.5 * dens*quad%aux%gravity%pot
# endif

# if ANGULAR_MOMENTUM
    angx = coords(:,:,:,Y_DIR)*momz - coords(:,:,:,Z_DIR)*momy
    angy = coords(:,:,:,Z_DIR)*momx - coords(:,:,:,X_DIR)*momz
    angz = coords(:,:,:,X_DIR)*momy - coords(:,:,:,Y_DIR)*momx

    angL = SQRT(angx**2 + angy**2 + angz**2)
    angH = angL/dens
# endif

    !$omp critical
    sums(idens) = sums(idens) + sum(vols * dens)
    sums(imomx) = sums(imomx) + sum(vols * momx)
    sums(imomy) = sums(imomy) + sum(vols * momy)
    sums(imomz) = sums(imomz) + sum(vols * momz)

    sums(ieint) = sums(ieint) + sum(vols * eint)
    sums(iekin) = sums(iekin) + sum(vols * ekin)
# if MODULE_SELF_GRAVITY
    sums(iepot) = sums(iepot) + sum(vols * epot)
# endif

    sums(ietot) = sums(ietot) + sum(vols * etot)

# if ANGULAR_MOMENTUM
    sums(iLtot) = sums(iLtot) + sum(vols * angL)
    sums(iHtot) = sums(iHtot) + sum(vols * angH)
# endif
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

        ! resBuf(incr(cnt)) = recv(ieint)
        ! resBuf(incr(cnt)) = recv(iekin)

! # if MODULE_SELF_GRAVITY
!         resBuf(incr(cnt)) = recv(iepot)
! # endif
! 
! # if ANGULAR_MOMENTUM
!         resBuf(incr(cnt)) = recv(iLtot)
!         resBuf(incr(cnt)) = recv(iHtot)
! # endif

        call write2file(fpath, cnt, resBuf)
    end if

end subroutine

end module
