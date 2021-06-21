module setup_config_mod

use constants_mod

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: init_time = 0.0_dp
real(dp), parameter :: stop_time = 0.2

real(dp), parameter :: CFL = 0.8_dp

real(dp), parameter :: beta = 0.8_dp

real(dp), parameter :: dt_checkpoint = 0.05_dp
!real(dp), parameter :: dt_checkpoint = 0.001_dp

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: dt_analyze = dt_checkpoint

real(dp), parameter :: chkpt_init_time = init_time
real(dp), parameter :: chkpt_stop_time = stop_time

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: diameter = 4.0_dp

real(dp), parameter :: xlim(2) = diameter * ((/0.0_dp,1.0_dp/) - 0.0_dp)
real(dp), parameter :: ylim(2) = diameter * ((/0.0_dp,1.0_dp/) - 0.0_dp)

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: dens0 = 1.0_dp
real(dp), parameter :: pres0 = 1.0_dp

real(dp), parameter :: velx0 = 0.0_dp
real(dp), parameter :: vely0 = 0.0_dp

real(dp), parameter :: dens_L =    8.0_dp
real(dp), parameter :: momx_L =   57.1597_dp
real(dp), parameter :: momy_L =  -33.0012_dp
real(dp), parameter :: ener_L =  563.544_dp

real(dp), parameter :: dens_R = 1.4_dp
real(dp), parameter :: momx_R = 0.0_dp
real(dp), parameter :: momy_R = 0.0_dp
real(dp), parameter :: ener_R = 2.5_dp

real(dp), parameter :: kappa = 1.4_dp

real(dp), parameter :: pivot = 1.0/6.0_dp
real(dp), parameter :: angle = 30.0_dp

!! -------------------------------------------------------------------------- !!

integer, parameter :: mesh_minlevel = 6
integer, parameter :: mesh_maxlevel = 6

real(dp), parameter :: setup_amr_threshold_refine = -6
real(dp), parameter :: setup_amr_threshold_coarsen = -7

real(dp), parameter :: dens_refine_threshold = 1.3_dp*dens0
real(dp), parameter :: dens_coarse_threshold = 1.2_dp*dens0

# if DO_PROFILING
integer, parameter :: prof_report_rate = 64
# endif

end module
