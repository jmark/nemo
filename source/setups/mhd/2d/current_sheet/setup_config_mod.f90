module setup_config_mod

use constants_mod

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: init_time = 0.0
real(dp), parameter :: stop_time = 5.0

real(dp), parameter :: CFL = 0.9

real(dp), parameter :: beta = 0.8

real(dp), parameter :: dt_checkpoint = 0.5

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: dt_analyze = 0.1

real(dp), parameter :: chkpt_init_time = init_time
real(dp), parameter :: chkpt_stop_time = stop_time

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: diameter = 2.0

real(dp), parameter :: xlim(2) = diameter * ((/0.0,1.0/) - 0.5)
real(dp), parameter :: ylim(2) = diameter * ((/0.0,1.0/) - 0.5)

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: kappa = 1.4_dp

real(dp), parameter :: dens0 = 1.0_dp

real(dp), parameter :: pres0 = 0.1_dp

real(dp), parameter :: velx0 = 0.1_dp
real(dp), parameter :: vely0 = 0.0_dp
real(dp), parameter :: velz0 = 0.0_dp

real(dp), parameter :: magx0 = 0.0_dp
real(dp), parameter :: magy0 = 1.0_dp
real(dp), parameter :: magz0 = 0.0_dp

real(dp), parameter :: glmp0 = 0.0_dp
 
real(dp), parameter :: mu0 = 1.0_dp
! real(dp), parameter :: GLM_ch = 1.0_dp
! real(dp), parameter :: GLM_alpha = GLM_ch / 0.18_dp

!! -------------------------------------------------------------------------- !!

integer, parameter :: mesh_minlevel = 5
integer, parameter :: mesh_maxlevel = 5

real(dp), parameter :: setup_amr_threshold_refine = -6
real(dp), parameter :: setup_amr_threshold_coarsen = -7

real(dp), parameter :: dens_refine_threshold = 1.3*dens0
real(dp), parameter :: dens_coarse_threshold = 1.2*dens0

# if DO_PROFILING
integer, parameter :: prof_report_rate = 64
# endif

end module
