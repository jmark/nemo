!# define PP_N_NODES        8
!# define PP_KERNEL_DGFV_QMAX           2
!# define PP_MESH_PERIODIC  0
!# define PP_SPLIT_FORM     0

!# define DO_PROFILING      0

# define PP_RIEMANN_INNER     riemann_hlle_mod
# define PP_RIEMANN_BOUNDARY  riemann_roe_entropy_fix_mod

module setup_config_mod

use constants_mod

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: init_time = 0.0
real(dp), parameter :: stop_time = 3.0

real(dp), parameter :: CFL = 0.8

real(dp), parameter :: beta = 0.8

real(dp), parameter :: dt_checkpoint = 5e-1

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: dt_analyze = 0.05*dt_checkpoint

real(dp), parameter :: chkpt_init_time = init_time
real(dp), parameter :: chkpt_stop_time = stop_time

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: diameter = 6

real(dp), parameter :: xlim(2) = diameter * ((/0.0,1.0/) - 0.5)
real(dp), parameter :: ylim(2) = diameter * ((/0.0,1.0/) - 0.5)
real(dp), parameter :: zlim(2) = diameter * ((/0.0,1.0/) - 0.5)

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: kappa = 5.0_dp/3.0_dp

real(dp), parameter :: grav  = 1.0
real(dp), parameter :: mass0 = 1.0
real(dp), parameter :: radi0 = 1.0

real(dp), parameter :: dens0 = 1e-5
real(dp), parameter :: ener0 = 5e-2

real(dp), parameter :: velx0 = 0.0
real(dp), parameter :: vely0 = 0.0
real(dp), parameter :: velz0 = 0.0
 
!! -------------------------------------------------------------------------- !!

integer, parameter :: mesh_minlevel = 3
integer, parameter :: mesh_maxlevel = 6

real(dp), parameter :: setup_amr_threshold_refine = -6
real(dp), parameter :: setup_amr_threshold_coarsen = -7

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: plummer = 1e-16_dp * diameter

# if DO_PROFILING
integer, parameter :: prof_report_rate = 16
# endif

end module
