module setup_config_mod

use constants_mod

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: init_time = 0.0
real(dp), parameter :: stop_time = 1000.0
!real(dp), parameter :: stop_time = 0.0

real(dp), parameter :: CFL = 0.8

real(dp), parameter :: beta = 0.8

real(dp), parameter :: dt_checkpoint = 1e2

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: dt_analyze = 1e-2

real(dp), parameter :: chkpt_init_time = init_time
real(dp), parameter :: chkpt_stop_time = stop_time

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: diameter = 1.0

real(dp), parameter :: xlim(2) = diameter * ((/0.0,1.0/) - 0.5)
real(dp), parameter :: ylim(2) = diameter * ((/0.0,1.0/) - 0.5)

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: bell_mass  = 1.5
real(dp), parameter :: bell_sigma = 0.1

real(dp), parameter :: blast_energy = 1.0
real(dp), parameter :: blast_sigma  = 0.25

real(dp), parameter :: bell_mass_normalized = bell_mass/SQRT((2*pi)**N_DIMS)/bell_sigma**N_DIMS
real(dp), parameter :: blast_ener_normalized = blast_energy/SQRT((2*pi)**N_DIMS)/blast_sigma**N_DIMS

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: dens0 = 1e-0
real(dp), parameter :: pres0 = 1e-0

real(dp), parameter :: velx0 = 0.0_dp
real(dp), parameter :: vely0 = 0.0_dp
 
real(dp), parameter :: kappa = 7.0_dp/5.0_dp

!! -------------------------------------------------------------------------- !!

integer, parameter :: mesh_minlevel = 3
integer, parameter :: mesh_maxlevel = 4

real(dp), parameter :: setup_amr_threshold_refine = -6
real(dp), parameter :: setup_amr_threshold_coarsen = -7

real(dp), parameter :: dens_refine_threshold = 1.3*dens0
real(dp), parameter :: dens_coarse_threshold = 1.2*dens0

# if DO_PROFILING
integer, parameter :: prof_report_rate = 64
# endif

end module
