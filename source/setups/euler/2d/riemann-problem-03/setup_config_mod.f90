module setup_config_mod

use constants_mod

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: init_time = 0.0_dp
real(dp), parameter :: stop_time = 0.3_dp
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


!! weak problem
! real(dp), parameter :: dens_UL = 0.838_dp
! real(dp), parameter :: velx_UL = 0.206_dp
! real(dp), parameter :: vely_UL = 0.206_dp
! real(dp), parameter :: pres_UL = 0.629_dp
! 
! real(dp), parameter :: dens_LL = 0.9323_dp
! real(dp), parameter :: velx_LL = 0.0_dp
! real(dp), parameter :: vely_LL = 0.206_dp
! real(dp), parameter :: pres_LL = 0.8_dp
! 
! real(dp), parameter :: dens_UR = 0.9323_dp
! real(dp), parameter :: velx_UR = 0.206_dp
! real(dp), parameter :: vely_UR = 0.0_dp
! real(dp), parameter :: pres_UR = 0.8_dp
! 
! real(dp), parameter :: dens_LR = 1.25_dp
! real(dp), parameter :: velx_LR = 0.0_dp
! real(dp), parameter :: vely_LR = 0.0_dp
! real(dp), parameter :: pres_LR = 1.0_dp

real(dp), parameter :: dens_UL = 0.138_dp
real(dp), parameter :: velx_UL = 1.206_dp
real(dp), parameter :: vely_UL = 1.206_dp
real(dp), parameter :: pres_UL = 0.029_dp

real(dp), parameter :: dens_LL = 0.5323_dp
real(dp), parameter :: velx_LL = 0.0_dp
real(dp), parameter :: vely_LL = 1.206_dp
real(dp), parameter :: pres_LL = 0.3_dp

real(dp), parameter :: dens_UR = 0.5323_dp
real(dp), parameter :: velx_UR = 1.206_dp
real(dp), parameter :: vely_UR = 0.0_dp
real(dp), parameter :: pres_UR = 0.3_dp

real(dp), parameter :: dens_LR = 1.5_dp
real(dp), parameter :: velx_LR = 0.0_dp
real(dp), parameter :: vely_LR = 0.0_dp
real(dp), parameter :: pres_LR = 1.5_dp

real(dp), parameter :: kappa = 1.4_dp

real(dp), parameter :: pivot = 0.0_dp

!! -------------------------------------------------------------------------- !!

integer, parameter :: mesh_minlevel = 5
integer, parameter :: mesh_maxlevel = 5

real(dp), parameter :: setup_amr_threshold_refine = -6
real(dp), parameter :: setup_amr_threshold_coarsen = -7

real(dp), parameter :: dens_refine_threshold = 1.3_dp*dens0
real(dp), parameter :: dens_coarse_threshold = 1.2_dp*dens0

# if DO_PROFILING
integer, parameter :: prof_report_rate = 64
# endif

end module
