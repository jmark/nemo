module setup_config_mod

use constants_mod

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: init_time =  10.0
real(dp), parameter :: stop_time = 905.0
!real(dp), parameter :: stop_time = 705.0
!real(dp), parameter :: stop_time = 520.0

real(dp), parameter :: CFL = 0.8_dp

real(dp), parameter :: beta = 0.8_dp

real(dp), parameter :: dt_checkpoint = 50.0_dp

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: dt_analyze = dt_checkpoint

real(dp), parameter :: chkpt_init_time = init_time
real(dp), parameter :: chkpt_stop_time = stop_time

!! -------------------------------------------------------------------------- !!

!! Ferrand (2012)
real(dp), parameter :: blast_energy = 5.25166181297794e-05
real(dp), parameter :: ejecta_mass = 1.4
real(dp), parameter :: dens0 = 0.0024535904882592957
real(dp), parameter :: pres0 = 2.1309232076296036e-13
real(dp), parameter :: rc = 0.07906945269055088
real(dp), parameter :: dc = 386.34489645765456

real(dp), parameter :: velx0 = 0.0
real(dp), parameter :: vely0 = 0.0
 
real(dp), parameter :: kappa = 5.0_dp/3.0_dp

real(dp), parameter :: blast_sigma = 0.5_dp*rc

real(dp), parameter :: blast_ener_normalized = blast_energy/SQRT((2*pi)**N_DIMS)/blast_sigma**N_DIMS

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: diameter = 5.0

real(dp), parameter :: xlim(2) = diameter * ((/0.0,1.0/) - 0.0)
real(dp), parameter :: ylim(2) = diameter * ((/0.0,1.0/) - 0.0)

!! -------------------------------------------------------------------------- !!

integer, parameter :: mesh_minlevel = 6
integer, parameter :: mesh_maxlevel = 6

real(dp), parameter :: setup_amr_threshold_refine = -6
real(dp), parameter :: setup_amr_threshold_coarsen = -7

real(dp), parameter :: dens_refine_threshold = 1.3*dens0
real(dp), parameter :: dens_coarse_threshold = 1.2*dens0

# if DO_PROFILING
integer, parameter :: prof_report_rate = 64
# endif

end module
