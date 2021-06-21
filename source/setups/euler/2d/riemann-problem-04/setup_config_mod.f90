module setup_config_mod

use constants_mod

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: init_time = 0.0_dp
real(dp), parameter :: stop_time = 0.25_dp
!real(dp), parameter :: stop_time = 1.5_dp
!real(dp), parameter :: stop_time = 0.0

real(dp), parameter :: CFL = 0.8_dp

real(dp), parameter :: beta = 0.8_dp

real(dp), parameter :: dt_checkpoint = 0.05_dp
!real(dp), parameter :: dt_checkpoint = 0.001_dp

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: dt_analyze = dt_checkpoint

real(dp), parameter :: chkpt_init_time = init_time
real(dp), parameter :: chkpt_stop_time = stop_time

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: diameter = 2.0_dp

real(dp), parameter :: xlim(2) = diameter * ((/0.0_dp,1.0_dp/) - 0.5_dp)
real(dp), parameter :: ylim(2) = diameter * ((/0.0_dp,1.0_dp/) - 0.5_dp)

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: dens0 = 1.0_dp
real(dp), parameter :: pres0 = 1.0_dp

real(dp), parameter :: velx0 = 0.0_dp
real(dp), parameter :: vely0 = 0.0_dp

real(dp), parameter :: dens_UL =  1.1_dp
real(dp), parameter :: velx_UL =  0.8939_dp
real(dp), parameter :: vely_UL =  0.8939_dp
real(dp), parameter :: pres_UL =  1.1_dp

real(dp), parameter :: dens_LL =  0.5065_dp
real(dp), parameter :: velx_LL =  0.0_dp
real(dp), parameter :: vely_LL =  0.8939_dp
real(dp), parameter :: pres_LL =  0.35_dp

real(dp), parameter :: dens_UR =  0.5065_dp
real(dp), parameter :: velx_UR =  0.8939_dp
real(dp), parameter :: vely_UR =  0.0_dp
real(dp), parameter :: pres_UR =  0.35_dp

real(dp), parameter :: dens_LR =  1.1_dp
real(dp), parameter :: velx_LR =  0.0_dp
real(dp), parameter :: vely_LR =  0.0_dp
real(dp), parameter :: pres_LR =  1.1_dp

real(dp), parameter :: kappa = 1.4_dp

real(dp), parameter :: pivot = 0.0_dp

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: refine_radius = 0.5_dp

integer, parameter :: mesh_minlevel = 6
integer, parameter :: mesh_maxlevel = 6

# if DO_PROFILING
integer, parameter :: prof_report_rate = 64
# endif

end module
