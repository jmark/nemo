module setup_config_mod

use constants_mod

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: init_time = 0.00_dp
real(dp), parameter :: stop_time = 0.06_dp

real(dp), parameter :: CFL = 0.9_dp

real(dp), parameter :: dt_checkpoint = 1e-2_dp

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: dt_analyze = dt_checkpoint

real(dp), parameter :: chkpt_init_time = init_time
real(dp), parameter :: chkpt_stop_time = stop_time

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: diameter = 1.0_dp

real(dp), parameter :: xlim(2) = diameter * ((/0.0,1.0/) - 0.0)
real(dp), parameter :: ylim(2) = diameter * ((/0.0,1.0/) - 0.0)
real(dp), parameter :: zlim(2) = diameter * ((/0.0,1.0/) - 0.0)

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: mu0 = 1.0_dp
real(dp), parameter :: kappa = 5.0_dp/3.0_dp

!! -------------------------------------------------------------------------- !!

integer, parameter :: amr_each_ntimesteps = 8

integer, parameter :: mesh_minlevel = 2
integer, parameter :: mesh_maxlevel = 6

# if DO_PROFILING
integer, parameter :: prof_report_rate = 64
# endif

end module
