module kernel_timestep_mod

use constants_mod

real(dp), save :: dtMin

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

    integer, parameter :: N = N_NODES           !! # nodes of the block

    dtmin = HUGE(1.0_dp)

# if PP_EQN_MHD
    lmax_global = 0.0_dp
    vmax_global = 0.0_dp
# endif

    call mesh_iterate_quads(mesh,timestep_cb)

    maxtimestep = 0.5_dp*CFL/REAL(N,dp) * dtmin

    call mpi_allreduce(maxtimestep, rt%maxtimestep, 1, MPI_DOUBLE_PRECISION, MPI_MIN, rt%mpi%comm, ierr)

# if PP_EQN_MHD
    call mpi_allreduce(lmax_global,lmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,rt%mpi%comm,ierr)
    call mpi_allreduce(vmax_global,vmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,rt%mpi%comm,ierr)

    glm_ch = 2.0_dp*abs(lmax - vmax)
    !glm_ch = 0.0_dp
    glm_alpha = GLM_ch / 0.18_dp
# endif

end subroutine

subroutine timestep_cb(quad)

    use mesh_types_mod, only: quad_t
    use equations_mod, only: xlambdaMax
    use equations_mod, only: ylambdaMax
    use equations_const_mod

    ! use setup_sourceterm_mod, only: setup_sourceterm_dt

    type(quad_t), intent(inout) :: quad

# if PP_EQN_MHD
    real(dp) :: lMax,vmax
# endif

    real(dp) :: lmaxF(N_NODES,N_NODES)
    real(dp) :: lmaxG(N_NODES,N_NODES)
    real(dp) :: dt

    lmaxF = xlambdaMax(quad%pld%hydro%state)
    lmaxG = ylambdaMax(quad%pld%hydro%state)

# if PP_EQN_MHD
    associate(u => quad%pld%hydro%state)
    lmax = max(maxval(lmaxF),maxval(lmaxG))
    vmax = max(maxval(abs(u(:,:,MOMX_VAR)/u(:,:,DENS_VAR))),maxval(abs(u(:,:,MOMY_VAR)/u(:,:,DENS_VAR))))
    end associate
# endif

    dt = HUGE(1.0_dp)
    !call setup_sourceterm_dt(quad,dt)
    dt = min(dt,min(quad%delta(1)/maxval(lmaxF),quad%delta(2)/maxval(lmaxG)))

    ! write (*,*) quad%gloid, dt

    !$omp critical
    dtmin = min(dtmin,dt)

# if PP_EQN_MHD
    lmax_global = max(lmax_global,lmax)
    vmax_global = max(vmax_global,vmax)
# endif
    !$omp end critical

end subroutine

end module
