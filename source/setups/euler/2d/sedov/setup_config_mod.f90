module setup_config_mod

use constants_mod

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: init_time = 0.0_dp
real(dp), parameter :: stop_time = 0.05_dp

real(dp), parameter :: CFL = 0.8_dp

real(dp), parameter :: beta = 0.2_dp

real(dp), parameter :: dt_checkpoint = 1e-2_dp

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: dt_analyze = dt_checkpoint

real(dp), parameter :: chkpt_init_time = init_time
real(dp), parameter :: chkpt_stop_time = stop_time

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: diameter = 0.5_dp

real(dp), parameter :: xlim(2) = diameter * ((/0.0_dp,1.0_dp/) - 0.5_dp)
real(dp), parameter :: ylim(2) = diameter * ((/0.0_dp,1.0_dp/) - 0.5_dp)

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: blast_energy = 1.0_dp
real(dp), parameter :: blast_sigma = 5e-3_dp

real(dp), parameter :: blast_ener_normalized = blast_energy/SQRT((2*pi)**N_DIMS)/blast_sigma**N_DIMS

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: dens0 = 1e-0_dp
!real(dp), parameter :: pres0 = 1e-5_dp
real(dp), parameter :: pres0 = 1e-14_dp

real(dp), parameter :: velx0 = 0.0_dp
real(dp), parameter :: vely0 = 0.0_dp
 
real(dp), parameter :: kappa = 7.0_dp/5.0_dp

!! -------------------------------------------------------------------------- !!

integer, parameter :: mesh_minlevel = 5
integer, parameter :: mesh_maxlevel = 5

real(dp), parameter :: setup_amr_threshold_refine = -6
real(dp), parameter :: setup_amr_threshold_coarsen = -7

real(dp), parameter :: dens_refine_threshold = 1.3*dens0
real(dp), parameter :: dens_coarse_threshold = 1.2*dens0

real(dp), parameter :: init_amr_min_diameter = 0.2* 0.5_dp*diameter

# if DO_PROFILING
integer, parameter :: prof_report_rate = 64
# endif

end module
