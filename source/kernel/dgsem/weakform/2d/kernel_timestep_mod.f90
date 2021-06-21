module kernel_timestep_mod

use constants_mod
use kernel_share_mod, only: dtmin
use kernel_share_mod, only: lmax

integer(kind=8) :: iquad_dtmin

# if PP_EQN_MHD
!! global maximum wave speed
real(dp), save :: lmax_global

!! global maximum fluid speed
real(dp), save :: vmax_global
# endif

contains

subroutine hook_calc_timestep_200

    use setup_config_mod, only: CFL
    use share_mod, only: mesh, rt => runtime
    use mesh_utils_mod, only: mesh_iterate_quads
    use mpi_f08, only: mpi_allreduce, MPI_MIN, MPI_MAX, MPI_DOUBLE_PRECISION

# if PP_EQN_MHD
    use equations_share_mod, only: glm_ch
    use equations_share_mod, only: glm_alpha

    real(dp) :: lmax
    real(dp) :: vmax
# endif

    real(dp) :: maxtimestep
    integer :: ierr

    dtmin = HUGE(1.0_dp)
    lmax = 0.0_dp

    iquad_dtmin = 0

# if PP_EQN_MHD
    lmax_global = 0.0_dp
    vmax_global = 0.0_dp
# endif

    call mesh_iterate_quads(mesh,timestep_cb)

    maxtimestep = 0.5_dp* CFL / REAL(2*N_NODES-1,dp) * dtmin

    call mpi_allreduce(maxtimestep,rt%maxtimestep,1,MPI_DOUBLE_PRECISION,MPI_MIN,rt%mpi%comm,ierr)

    ! mesh%quads(iquad_dtmin)%aux%hydro%dtmin_count = mesh%quads(iquad_dtmin)%aux%hydro%dtmin_count + 1

# if PP_EQN_MHD
    call mpi_allreduce(lmax_global,lmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,rt%mpi%comm,ierr)
    call mpi_allreduce(vmax_global,vmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,rt%mpi%comm,ierr)

    glm_ch = 2.0_dp*abs(lmax - vmax)
    glm_alpha = glm_ch / 0.18_dp
# endif

end subroutine

subroutine timestep_cb(quad)

    use mesh_types_mod, only: quad_t
    use equations_mod, only: xlambdaMax
    use equations_mod, only: ylambdaMax
    use equations_const_mod

    ! use setup_sourceterm_mod, only: setup_sourceterm_dt

    type(quad_t), intent(inout) :: quad

    real(dp) :: lmaxF ! (N_NODES,N_NODES)
    real(dp) :: lmaxG ! (N_NODES,N_NODES)
    real(dp) :: dt

    real(dp) :: means(N_VARS)
    real(dp) :: pmins(N_VARS)
    real(dp) :: prim(N_NODES,N_NODES,N_VARS)
    real(dp) :: cons(N_NODES,N_NODES,N_VARS)

# if PP_EQN_MHD
    real(dp) :: lMax,vmax
# endif

    integer :: i,j,h

    lmaxF = maxval(xlambdaMax(quad%pld%hydro%state))
    lmaxG = maxval(ylambdaMax(quad%pld%hydro%state))

    dt = HUGE(1.0_dp)
    !call setup_sourceterm_dt(quad,dt)

    dt = min(dt,min(quad%delta(X_DIR)/lmaxF,quad%delta(Y_DIR)/lmaxG))

    ! write (*,*) quad%delta, maxval(lmaxF), maxval(lmaxG), dt

# if PP_EQN_MHD
    associate(u => quad%pld%hydro%state)
    lmax = max(lmaxF,lmaxG)
    vmax = max(maxval(abs(u(:,:,MOMX_VAR)/u(:,:,DENS_VAR))),maxval(abs(u(:,:,MOMY_VAR)/u(:,:,DENS_VAR))))
    end associate
# endif

    !$omp critical
    if (dt < dtmin) then
        iquad_dtmin = quad%locid
    end if
    dtmin = min(dtmin,dt)
    lmax = max(lmax,lmaxF,lmaxG)

# if PP_EQN_MHD
    lmax_global = max(lmax_global,lmax)
    vmax_global = max(vmax_global,vmax)
# endif
    !$omp end critical

end subroutine

end module
