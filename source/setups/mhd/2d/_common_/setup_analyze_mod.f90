# define ANGULAR_MOMENTUM 0

module setup_analyze_mod

use constants_mod

character(len=255) :: fpath = 'analysis.dat'

! integer, parameter :: idens      = 1
! 
! integer, parameter :: imomx      = 2
! integer, parameter :: imomy      = 3
! integer, parameter :: imomz      = 4
! 
! integer, parameter :: ieint      = 5
! integer, parameter :: iekin      = 6
! integer, parameter :: iepot      = 7
! integer, parameter :: ietot      = 8
! 
! integer, parameter :: idivb      = 9
! 
! integer, parameter :: sums_len  = 30

enum, bind(c)

    enumerator :: izero = 0

    enumerator :: idens

    enumerator :: imomx
    enumerator :: imomy
    enumerator :: imomz

    enumerator :: ieint
    enumerator :: iekin
    enumerator :: iepot
    enumerator :: ietot

    enumerator :: itmag
    enumerator :: idivb

    enumerator :: ilast

endenum

integer, parameter :: SUMS_LEN = ilast - 1

real(dp), save :: sums(1:SUMS_LEN) !! accumulators
real(dp), save :: nextanalyze = 0.0

contains

subroutine analyze_fill_halo(side)

    use mesh_types_mod, only: quad_t, side_t
    use mesh_const_mod, only: MESH_HALF, MESH_FACES

    use setup_boundary_mod, only: setup_side_boundary

    integer, parameter :: N = N_NODES
    integer, parameter :: V = N_VARS

    type(side_t), intent(inout) :: side

    integer :: iside, iquad

    real(dp), parameter :: dummy(N_NODES,N_VARS) = 0.0_dp

    if (side%boundary > 0) then

        call setup_side_boundary(side,dummy,side%quads(1,side%inner)%p%aux%analyze%halo(:,:,side%faceid(side%inner)))

    else

        !! Only support for uniform grid now!
        do iside = 1,2
            select case(side%faceid(3-iside)) !! opposite side
                case (NOR); do iquad = 1,1; side%quads(iquad,iside)%p%aux%analyze%halo(:,:,side%faceid(iside)) = side%quads(iquad,3-iside)%p%pld%hydro%state(1,:,:); end do
                case (SOU); do iquad = 1,1; side%quads(iquad,iside)%p%aux%analyze%halo(:,:,side%faceid(iside)) = side%quads(iquad,3-iside)%p%pld%hydro%state(N,:,:); end do
                case (WES); do iquad = 1,1; side%quads(iquad,iside)%p%aux%analyze%halo(:,:,side%faceid(iside)) = side%quads(iquad,3-iside)%p%pld%hydro%state(:,1,:); end do
                case (EAS); do iquad = 1,1; side%quads(iquad,iside)%p%aux%analyze%halo(:,:,side%faceid(iside)) = side%quads(iquad,3-iside)%p%pld%hydro%state(:,N,:); end do
            end select
        end do
    end if

end subroutine

subroutine analyze_cb(quad)

    use setup_config_mod
    use equations_const_mod

    use equations_mod, only: pressure
    use mesh_types_mod, only: quad_t
    use quad_utils_mod, only: quad_get_volumes
    use quad_utils_mod, only: quad_get_coords

    type(quad_t), intent(inout) :: quad
    
    integer, parameter :: N = N_NODES
    integer, parameter :: P = N_NODES+1
    integer, parameter :: V = N_VARS

    real(dp) :: coords(N_NODES,N_NODES,N_DIMS)

    real(dp) :: dens(N_NODES,N_NODES)

    real(dp) :: momx(N_NODES,N_NODES)
    real(dp) :: momy(N_NODES,N_NODES)
    real(dp) :: momz(N_NODES,N_NODES)

    real(dp) :: etot(N_NODES,N_NODES)
    real(dp) :: ekin(N_NODES,N_NODES)
    real(dp) :: eint(N_NODES,N_NODES)
    real(dp) :: epot(N_NODES,N_NODES)

    real(dp) :: tmag(N_NODES,N_NODES)

    real(dp) :: magx(0:N_NODES+1,0:N_NODES+1)
    real(dp) :: magy(0:N_NODES+1,0:N_NODES+1)
    real(dp) :: magz(0:N_NODES+1,0:N_NODES+1)

    real(dp) :: glmp(N_NODES,N_NODES)

    real(dp) :: divB(N_NODES,N_NODES)

    real(dp) :: vols(N_NODES,N_NODES)

    real(dp) :: sdx,sdy

    integer :: i,j

    vols = quad_get_volumes(quad)
    coords = quad_get_coords(quad)

    dens = quad%pld%hydro%state(:,:,DENS_VAR)
    momx = quad%pld%hydro%state(:,:,MOMX_VAR)
    momy = quad%pld%hydro%state(:,:,MOMY_VAR)
    etot = quad%pld%hydro%state(:,:,ENER_VAR)

    magx(1:N,1:N) = quad%pld%hydro%state(:,:,MAGX_VAR)
    magy(1:N,1:N) = quad%pld%hydro%state(:,:,MAGY_VAR)
    magz(1:N,1:N) = quad%pld%hydro%state(:,:,MAGZ_VAR)

    magx(0,1:N) = quad%aux%analyze%halo(:,MAGX_VAR,NOR)
    magx(P,1:N) = quad%aux%analyze%halo(:,MAGX_VAR,SOU)
    magx(1:N,0) = quad%aux%analyze%halo(:,MAGX_VAR,WES)
    magx(1:N,P) = quad%aux%analyze%halo(:,MAGX_VAR,EAS)

    magy(0,1:N) = quad%aux%analyze%halo(:,MAGY_VAR,NOR)
    magy(P,1:N) = quad%aux%analyze%halo(:,MAGY_VAR,SOU)
    magy(1:N,0) = quad%aux%analyze%halo(:,MAGY_VAR,WES)
    magy(1:N,P) = quad%aux%analyze%halo(:,MAGY_VAR,EAS)

    magz(0,1:N) = quad%aux%analyze%halo(:,MAGZ_VAR,NOR)
    magz(P,1:N) = quad%aux%analyze%halo(:,MAGZ_VAR,SOU)
    magz(1:N,0) = quad%aux%analyze%halo(:,MAGZ_VAR,WES)
    magz(1:N,P) = quad%aux%analyze%halo(:,MAGZ_VAR,EAS)

    sdx = 0.5_dp*N_NODES/quad%delta(1)
    sdy = 0.5_dp*N_NODES/quad%delta(2)

    do j = 1,N; do i = 1,N;
        divB(i,j) = sdx*(magx(i+1,j) - magx(i-1,j)) + sdy*(magy(i,j+1) - magy(i,j-1))
    end do; end do;

    tmag = sqrt(magx(1:N,1:N)**2 + magy(1:N,1:N)**2 + magz(1:N,1:N)**2)

    eint = pressure(quad%pld%hydro%state)/(kappa-1.0_dp)
    ekin = 0.5 * (momx**2 + momy**2)/dens
# if MODULE_SELF_GRAVITY
    epot = 0.5 * dens*quad%aux%gravity%pot
# endif

    !$omp critical
    sums(idens) = sums(idens) + sum(vols * dens)
    sums(imomx) = sums(imomx) + sum(vols * momx)
    sums(imomy) = sums(imomy) + sum(vols * momy)

    sums(ieint) = sums(ieint) + sum(vols * eint)
    sums(iekin) = sums(iekin) + sum(vols * ekin)

    sums(iekin) = sums(iekin) + sum(vols * ekin)

# if MODULE_SELF_GRAVITY
    sums(iepot) = sums(iepot) + sum(vols * epot)
# endif

    sums(ietot) = sums(ietot) + sum(vols * etot)

    sums(itmag) = sums(itmag) + sum(vols * tmag)
    sums(idivb) = sums(idivb) + sum(vols * abs(divb))
    !$omp end critical

end subroutine

subroutine hook_analyze_500

    use share_mod, only: mesh, rt => runtime

    use mesh_utils_mod, only: mesh_exchange
    use mesh_utils_mod, only: mesh_iterate_quads
    use mesh_utils_mod, only: mesh_iterate_sides

    use mpi_f08, only: mpi_reduce,MPI_SUM,MPI_DOUBLE_PRECISION

    use utils_mod, only: incr, write2file

    use setup_config_mod, only: dt_analyze

    use mesh_const_mod, only: MESH_ITERATE_WITH_GHOSTS

    use equations_share_mod, only: glm_ch

    integer :: ierr
    integer :: cnt = 0
    real(dp) :: resBuf(99)
    real(dp) :: recv(1:SUMS_LEN) !! mpi accumulator

    if (rt%simtime < nextanalyze) return

    nextanalyze = rt%simtime + dt_analyze

    call mesh_exchange(mesh)
    call mesh_iterate_sides(mesh,analyze_fill_halo)

    sums = 0.0
 
    call mesh_iterate_quads(mesh,analyze_cb)

    call mpi_reduce(sums,recv,SUMS_LEN,MPI_DOUBLE_PRECISION,MPI_SUM,0,rt%mpi%comm,ierr)

    if (rt%mpi%root) then
        cnt = 0

        resBuf(incr(cnt)) = rt%simtime
        resBuf(incr(cnt)) = rt%maxtimestep
        resBuf(incr(cnt)) = recv(idens)
        resBuf(incr(cnt)) = recv(ietot)
        resBuf(incr(cnt)) = recv(itmag)
        resBuf(incr(cnt)) = recv(idivb)
        resBuf(incr(cnt)) = glm_ch

        call write2file(fpath,cnt,resBuf)
    end if

end subroutine

end module
