module setup_config_mod

use constants_mod

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: init_time = 0.0
real(dp), parameter :: stop_time = 0.0

real(dp), parameter :: CFL = 0.9

real(dp), parameter :: beta = 0.8

real(dp), parameter :: dt_checkpoint = 1e-1

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: dt_analyze = dt_checkpoint

real(dp), parameter :: chkpt_init_time = init_time
real(dp), parameter :: chkpt_stop_time = stop_time

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: diameter = 1 * 0.2
!real(dp), parameter :: diameter = 2

real(dp), parameter :: xlim(2) = diameter * ((/0.0,1.0/) - 0.5)
real(dp), parameter :: ylim(2) = diameter * ((/0.0,1.0/) - 0.5)
real(dp), parameter :: zlim(2) = diameter * ((/0.0,1.0/) - 0.5)

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: dens0 = 1.2558834406799149
real(dp), parameter :: pres0 = 54.608515091988004

real(dp), parameter :: velx0 = 0.0
real(dp), parameter :: vely0 = 0.0
real(dp), parameter :: velz0 = 0.0
 
real(dp), parameter :: kappa = 5.0_dp/3.0_dp

!! -------------------------------------------------------------------------- !!

integer, parameter :: mesh_minlevel = 3
integer, parameter :: mesh_maxlevel = 7

real(dp), parameter :: setup_amr_threshold_refine = -6
real(dp), parameter :: setup_amr_threshold_coarsen = -7

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: grav = 0.004498512523396364
real(dp), parameter :: plummer = 1e-16_dp * diameter

end module
