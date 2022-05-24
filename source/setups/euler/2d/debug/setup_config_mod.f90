module setup_config_mod

use constants_mod

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: init_time = 0.0_dp
!real(dp), parameter :: stop_time = 2.00_dp
!real(dp), parameter :: stop_time = 0.25_dp
real(dp), parameter :: stop_time = 3.7_dp
!real(dp), parameter :: stop_time = 0.0_dp

real(dp), parameter :: CFL = 0.9_dp

real(dp), parameter :: dt_checkpoint = 0.5_dp

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: dt_analyze = dt_checkpoint

real(dp), parameter :: chkpt_init_time = init_time
real(dp), parameter :: chkpt_stop_time = stop_time

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: diameter = 1.0_dp

real(dp), parameter :: xlim(2) = diameter * ((/0.0_dp,1.0_dp/) - 0.0_dp)
real(dp), parameter :: ylim(2) = diameter * ((/0.0_dp,1.0_dp/) - 0.0_dp)

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: dens0 = 0.5_dp
real(dp), parameter :: dens1 = 1.5_dp
real(dp), parameter :: pres0 = 1.0_dp
real(dp), parameter :: velx0 = 0.5_dp
real(dp), parameter :: vely0 = 0.1_dp
real(dp), parameter :: slope = 15.0_dp

real(dp), parameter :: kappa = 1.4_dp
! real(dp), parameter :: mu0 = 1.0_dp

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: refine_radius = 0.75_dp

integer, parameter :: mesh_minlevel = 1
integer, parameter :: mesh_maxlevel = mesh_minlevel

# if DO_PROFILING
integer, parameter :: prof_report_rate = 64
# endif

end module
