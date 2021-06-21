module setup_config_mod

use constants_mod

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: init_time = 0.0_dp
real(dp), parameter :: stop_time = 1.78_dp

real(dp), parameter :: CFL = 0.9_dp

real(dp), parameter :: dt_checkpoint = 0.4_dp

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: dt_analyze = dt_checkpoint

real(dp), parameter :: chkpt_init_time = init_time
real(dp), parameter :: chkpt_stop_time = stop_time

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: diameter = 9.0_dp

real(dp), parameter :: xlim(2) = diameter * ((/0.0_dp,1.0_dp/) - 0.5_dp)
real(dp), parameter :: ylim(2) = diameter * ((/0.0_dp,1.0_dp/) - 0.5_dp)

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: epsil0 = 0.2_dp
real(dp), parameter :: defreq = 5.0_dp

real(dp), parameter :: pivot = -4.0_dp

real(dp), parameter :: densL = 3.8567143_dp
real(dp), parameter :: presL = 10.33333_dp

real(dp), parameter :: densR = 1.0_dp + epsil0*sin(defreq*xlim(2))
real(dp), parameter :: presR = 1.0_dp

real(dp), parameter :: velxL = 2.629369_dp
real(dp), parameter :: velxR = 0.0_dp

real(dp), parameter :: vely0 = 0.0_dp

real(dp), parameter :: kappa = 1.4_dp

!! -------------------------------------------------------------------------- !!

integer, parameter :: mesh_minlevel = 5
integer, parameter :: mesh_maxlevel = 5

real(dp), parameter :: setup_amr_threshold_refine = -6
real(dp), parameter :: setup_amr_threshold_coarsen = -7

# if DO_PROFILING
integer, parameter :: prof_report_rate = 64
# endif

end module
