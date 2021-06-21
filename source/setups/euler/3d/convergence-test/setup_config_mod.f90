module setup_config_mod

use constants_mod

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: init_time = 0.0
real(dp), parameter :: stop_time = 0.01

real(dp), parameter :: CFL = 0.8

real(dp), parameter :: beta = 0.8

real(dp), parameter :: dt_checkpoint = 1e-1

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: dt_analyze = 1e-3

real(dp), parameter :: chkpt_init_time = init_time
real(dp), parameter :: chkpt_stop_time = stop_time

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: diameter = 1.0

real(dp), parameter :: xlim(2) = diameter * ((/0.0,1.0/) - 0.0)
real(dp), parameter :: ylim(2) = diameter * ((/0.0,1.0/) - 0.0)
real(dp), parameter :: zlim(2) = diameter * ((/0.0,1.0/) - 0.0)

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: bell_mass  = 1.5
real(dp), parameter :: bell_sigma = 0.1

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: dens0 = 1e-0
real(dp), parameter :: pres0 = 1e-0

real(dp), parameter :: velx0 = 0.0_dp
real(dp), parameter :: vely0 = 0.0_dp
real(dp), parameter :: velz0 = 0.0_dp
 
real(dp), parameter :: kappa = 7.0_dp/5.0_dp

!! -------------------------------------------------------------------------- !!

integer, parameter :: mesh_minlevel = 4
integer, parameter :: mesh_maxlevel = 4

real(dp), parameter :: setup_amr_threshold_refine = -6
real(dp), parameter :: setup_amr_threshold_coarsen = -7

real(dp), parameter :: dens_refine_threshold = 1.3*dens0
real(dp), parameter :: dens_coarse_threshold = 1.2*dens0

# if DO_PROFILING
integer, parameter :: prof_report_rate = 64
# endif

!! -------------------------------------------------------------------------- !!

! real(dp), parameter :: d0 = 1.0, d1 = 1.0 * 0.05, d2 = 1.0 * 0.04, d3 = 0.0 * 0.03, d4 = 0.0 ! 44;
! real(dp), parameter :: u0 = 0.1, u1 = 0.0 * 0.03, u2 = 0.0 * 0.02, u3 = 0.0 * 0.01, u4 = 0.0 ! 31;
! real(dp), parameter :: v0 = 0.2, v1 = 0.0 * 0.04, v2 = 0.0 * 0.03, v3 = 0.0 * 0.06, v4 = 0.0 ! 24;
! real(dp), parameter :: p0 = 1.0, p1 = 1.0 * 0.03, p2 = 1.0 * 0.09, p3 = 0.0 * 0.02, p4 = 0.0 ! 23;

! !! amplitudes
! real(dp), parameter :: d0 = 1.0_dp, d1 = 1.0_dp * 0.05_dp, d2 = 1.0_dp * 0.04_dp, d3 = 1.0_dp * 0.05_dp, d4 = 0.0_dp ! 44;
! real(dp), parameter :: u0 = 0.1_dp, u1 = 0.0_dp * 0.03_dp, u2 = 0.0_dp * 0.02_dp, u3 = 0.0_dp * 0.01_dp, u4 = 0.0_dp ! 31;
! real(dp), parameter :: v0 = 0.2_dp, v1 = 0.0_dp * 0.04_dp, v2 = 0.0_dp * 0.03_dp, v3 = 0.0_dp * 0.04_dp, v4 = 0.0_dp ! 24;
! real(dp), parameter :: w0 = 0.3_dp, w1 = 0.0_dp * 0.01_dp, w2 = 0.0_dp * 0.01_dp, w3 = 0.0_dp * 0.02_dp, w4 = 0.0_dp ! 12;
! real(dp), parameter :: p0 = 1.0_dp, p1 = 1.0_dp * 0.07_dp, p2 = 1.0_dp * 0.09_dp, p3 = 1.0_dp * 0.02_dp, p4 = 0.0_dp ! 23;

!! amplitudes
real(dp), parameter :: d0 = 1.0_dp, d1 = 1.0_dp * 0.35_dp, d2 = 1.0_dp * 0.24_dp, d3 = 1.0_dp * 0.12_dp, d4 = 0.0_dp ! 44;
real(dp), parameter :: u0 = 0.1_dp ! , u1 = 0.0_dp * 0.03_dp, u2 = 0.0_dp * 0.02_dp, u3 = 0.0_dp * 0.01_dp, u4 = 0.0_dp ! 31;
real(dp), parameter :: v0 = 0.2_dp ! , v1 = 0.0_dp * 0.04_dp, v2 = 0.0_dp * 0.03_dp, v3 = 0.0_dp * 0.04_dp, v4 = 0.0_dp ! 24;
real(dp), parameter :: w0 = 0.3_dp ! , w1 = 0.0_dp * 0.01_dp, w2 = 0.0_dp * 0.01_dp, w3 = 0.0_dp * 0.02_dp, w4 = 0.0_dp ! 12;
real(dp), parameter :: p0 = 1.0_dp, p1 = 1.0_dp * 0.23_dp, p2 = 1.0_dp * 0.19_dp, p3 = 1.0_dp * 0.21_dp, p4 = 0.0_dp ! 23;

!! frequencies
! real(dp), parameter :: df1 = 1.0_dp, df2 = 1.0_dp, df3 = 1.0_dp, df4 = 1.0_dp;
! real(dp), parameter :: uf1 = 1.0_dp, uf2 = 1.0_dp, uf3 = 1.0_dp, uf4 = 1.0_dp;
! real(dp), parameter :: vf1 = 1.0_dp, vf2 = 1.0_dp, vf3 = 1.0_dp, vf4 = 1.0_dp;
! real(dp), parameter :: wf1 = 1.0_dp, wf2 = 1.0_dp, wf3 = 1.0_dp, wf4 = 1.0_dp;
! real(dp), parameter :: pf1 = 1.0_dp, pf2 = 1.0_dp, pf3 = 1.0_dp, pf4 = 1.0_dp;

real(dp), parameter :: L = 1.0_dp;

end module
