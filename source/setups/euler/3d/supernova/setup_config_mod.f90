module setup_config_mod

use constants_mod

!! -------------------------------------------------------------------------- !!

!! !! Fraschetti (2010)
!! real(dp), parameter :: blast_energy = 8.402658900764704e-05
!! real(dp), parameter :: ejecta_mass = 5.0
!! real(dp), parameter :: dens0 = 0.010305080050689041
!! real(dp), parameter :: pres0 = 5.369926483226601e-13
!! real(dp), parameter :: rc = 0.052923400308888896
!! real(dp), parameter :: dc = 4601.507427022113

!! Ferrand (2012)
real(dp), parameter :: blast_energy = 5.25166181297794e-05_dp
real(dp), parameter :: ejecta_mass = 1.4_dp
real(dp), parameter :: dens0 = 0.0024535904882592957_dp
real(dp), parameter :: pres0 = 2.1309232076296036e-13_dp
real(dp), parameter :: rc = 0.07906945269055088_dp
real(dp), parameter :: dc = 386.34489645765456_dp

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: init_time =   0.0_dp !! years
real(dp), parameter :: stop_time = 500.0_dp !! years

real(dp), parameter :: chkpt_init_time = 0.0_dp !! years
real(dp), parameter :: chkpt_stop_time = stop_time !! years

!real(dp), parameter :: dt_checkpoint = (chkpt_stop_time-chkpt_init_time)/25
real(dp), parameter :: dt_checkpoint = 10.0_dp

real(dp), parameter :: dt_analyze = dt_checkpoint/10

integer, parameter :: nt_analyze = 16

real(dp), parameter :: CFL = 0.8_dp

real(dp), parameter :: beta = 0.8_dp

real(dp), parameter :: switch_time = 50.0_dp

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: diameter = 5.0 !! parsec

real(dp), parameter :: xlim(2) = diameter * ((/0.0,1.0/) - 0.0)
real(dp), parameter :: ylim(2) = diameter * ((/0.0,1.0/) - 0.0)
real(dp), parameter :: zlim(2) = diameter * ((/0.0,1.0/) - 0.0)

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: blast_sigma = 0.75_dp*rc
real(dp), parameter :: blast_ener_normalized = blast_energy/SQRT((2*pi)**N_DIMS)/blast_sigma**N_DIMS

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: velx0 = 0.0_dp
real(dp), parameter :: vely0 = 0.0_dp
real(dp), parameter :: velz0 = 0.0_dp
 
real(dp), parameter :: kappa = 5.0_dp/3.0_dp

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: r_offset = 0.1_dp*diameter

integer, parameter :: amr_each_ntimesteps = 16

integer, parameter :: mesh_minlevel = 2
integer, parameter :: mesh_maxlevel = 6

real(dp), parameter :: setup_amr_threshold_refine = -6
real(dp), parameter :: setup_amr_threshold_coarsen = -7

# if DO_PROFILING
integer, parameter :: prof_report_rate = 64
# endif

end module
