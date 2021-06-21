module setup_config_mod

use constants_mod

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: init_time = 0.0_dp
real(dp), parameter :: stop_time = 0.15_dp
!real(dp), parameter :: stop_time = 0.0_dp

real(dp), parameter :: CFL = 0.8_dp

real(dp), parameter :: dt_checkpoint = 0.05_dp

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: dt_analyze = dt_checkpoint

real(dp), parameter :: chkpt_init_time = init_time
real(dp), parameter :: chkpt_stop_time = stop_time

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: diameter = 1.0_dp

real(dp), parameter :: xlim(2) = diameter * ((/0.0_dp,1.0_dp/) - 0.0_dp)
real(dp), parameter :: ylim(2) = diameter * ((/0.0_dp,1.0_dp/) - 0.0_dp)

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: densL = 0.445_dp
real(dp), parameter :: velxL = 0.698_dp
real(dp), parameter :: presL = 3.528_dp

real(dp), parameter :: densR = 0.5_dp
real(dp), parameter :: velxR = 0.0_dp
real(dp), parameter :: presR = 0.571_dp

real(dp), parameter :: velyL = 0.0_dp
real(dp), parameter :: velyR = 0.0_dp

real(dp), parameter :: kappa = 1.4_dp

real(dp), parameter :: pivot = 0.5_dp

!! -------------------------------------------------------------------------- !!

integer, parameter :: mesh_minlevel = 7
integer, parameter :: mesh_maxlevel = 7

real(dp), parameter :: setup_amr_threshold_refine = -6
real(dp), parameter :: setup_amr_threshold_coarsen = -7

# if DO_PROFILING
integer, parameter :: prof_report_rate = 64
# endif

end module
