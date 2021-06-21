module setup_config_mod

use constants_mod

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: init_time = 0.0
real(dp), parameter :: stop_time = 0.05
!real(dp), parameter :: stop_time = 0.1
!real(dp), parameter :: stop_time = 0.0

real(dp), parameter :: CFL = 0.8

real(dp), parameter :: beta = 0.8

real(dp), parameter :: dt_checkpoint = 1e-1

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: dt_analyze = dt_checkpoint

real(dp), parameter :: chkpt_init_time = init_time
real(dp), parameter :: chkpt_stop_time = stop_time

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: diameter = 1.0

real(dp), parameter :: xlim(2) = diameter * ((/0.0,1.0/) - 0.5)
real(dp), parameter :: ylim(2) = diameter * ((/0.0,1.0/) - 0.5)
real(dp), parameter :: zlim(2) = diameter * ((/0.0,1.0/) - 0.5)

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: bell_mass  = 0.5_dp
real(dp), parameter :: bell_sigma = 0.1_dp

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: dens0 = 1e-0
real(dp), parameter :: pres0 = 1e-0

real(dp), parameter :: velx0 = 0.5_dp
real(dp), parameter :: vely0 = 0.5_dp
real(dp), parameter :: velz0 = 0.5_dp
 
real(dp), parameter :: magx0 = 0.0_dp
real(dp), parameter :: magy0 = 0.0_dp
real(dp), parameter :: magz0 = 0.0_dp

real(dp), parameter :: glmp0 = 0.0_dp
 
real(dp), parameter :: kappa = 7.0_dp/5.0_dp

real(dp), parameter :: mu0 = 1.0_dp

!! -------------------------------------------------------------------------- !!

integer, parameter :: amr_each_ntimesteps = 8

integer, parameter :: mesh_minlevel = 6
integer, parameter :: mesh_maxlevel = mesh_minlevel

real(dp), parameter :: setup_amr_threshold_refine = -6
real(dp), parameter :: setup_amr_threshold_coarsen = -7

real(dp), parameter :: dens_refine_threshold = 1.3*dens0
real(dp), parameter :: dens_coarse_threshold = 1.2*dens0

# if DO_PROFILING
integer, parameter :: prof_report_rate = 64
# endif

end module
