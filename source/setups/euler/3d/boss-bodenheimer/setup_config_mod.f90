module setup_config_mod

use constants_mod

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: init_time = 0.0
!real(dp), parameter :: stop_time = 0.08 !! Myr
real(dp), parameter :: stop_time = 0.0 !! Myr

real(dp), parameter :: CFL = 0.8

real(dp), parameter :: beta = 0.1

real(dp), parameter :: dt_checkpoint = 2e-4_dp

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: dt_analyze = dt_checkpoint

real(dp), parameter :: chkpt_init_time = init_time
real(dp), parameter :: chkpt_stop_time = stop_time

!! -------------------------------------------------------------------------- !!

# if 0
!! xmax    =  6.25e+16
!! xmin    = -6.25e+16
!! ymax    =  6.25e+16
!! ymin    = -6.25e+16
!! zmax    =  6.25e+16
!! zmin    = -6.25e+16
!! sim_rho0 = 3.82e-18
!! sim_rc = 5.e+16
!! sim_xCenter = 0.0
!! sim_yCenter = 0.0
!! sim_zCenter = 0.0
!! sim_omega = 7.2e-13
!! # Mabye just stick with 2.0e+4
!! sim_cs = 1.66e+4
!! sim_densFluc = 0.1

real(dp), parameter :: cloudradius = 5.e+16
real(dp), parameter :: diameter = 2 * 6.25e+16
# endif

real(dp), parameter :: cloudradius = 0.031625014819641306 !! lyr

!real(dp), parameter :: diameter = 0.031625014819641306
real(dp), parameter :: diameter = 2.0_dp * 1.2_dp * cloudradius
!real(dp), parameter :: diameter = 1.0

real(dp), parameter :: xlim(2) = diameter * ((/0.0,1.0/) - 0.5)
real(dp), parameter :: ylim(2) = diameter * ((/0.0,1.0/) - 0.5)
real(dp), parameter :: zlim(2) = diameter * ((/0.0,1.0/) - 0.5)

real(dp), parameter :: domsize(N_DIMS) = (/xlim(2)-xlim(1),ylim(2)-ylim(1),zlim(2)-zlim(1)/)
!! -------------------------------------------------------------------------- !!

!real(dp), parameter :: kappa = 5.0_dp/3.0_dp
real(dp), parameter :: kappa = 7.0_dp/5.0_dp

! real(dp), parameter :: xofs = diameter*0.25
! real(dp), parameter :: yofs = diameter*0.25
! real(dp), parameter :: zofs = diameter*0.25

real(dp), parameter :: bell_sigma = 0.25*cloudradius
real(dp), parameter :: xofs = 0.75*cloudradius

!! -------------------------------------------------------------------------- !!

!real(dp), parameter :: dens_cloud = 3.82e-18
real(dp), parameter :: dens_cloud = 6813.539303377171

real(dp), parameter :: dens0 = dens_cloud * 1e-5
!real(dp), parameter :: dens0 = 1.0_dp

real(dp), parameter :: dens_min = dens0 * 1e-2

real(dp), parameter :: angl0 = 1.25*49.229856

real(dp), parameter :: velx0 = 0.0_dp
real(dp), parameter :: vely0 = 0.0_dp
real(dp), parameter :: velz0 = 0.0_dp

!real(dp), parameter :: R_specific = kB/(mu*mH)
!real(dp), parameter :: dens_critical = 42584620.646107316

real(dp), parameter :: dens_critical = dens_cloud * 100.0_dp
real(dp), parameter :: temp_min = 10
!real(dp), parameter :: snds_min = 0.9618257275276481 !! lightyear/Myr
real(dp), parameter :: snds_min = 0.6801134942544845 !! lightyear/Myr

real(dp), parameter :: K_isothermal = snds_min**2
real(dp), parameter :: K_adiabatic = snds_min**2 * dens_critical**(1-kappa)

!! -------------------------------------------------------------------------- !!

integer, parameter :: mesh_minlevel = 2
integer, parameter :: mesh_maxlevel = 8

real(dp), parameter :: setup_amr_threshold_refine = -6
real(dp), parameter :: setup_amr_threshold_coarsen = -7

!! -------------------------------------------------------------------------- !!

real(dp), parameter :: grav = 0.15607787865870767
real(dp), parameter :: plummer = 1e-9_dp*diameter/2**(mesh_maxlevel)

!! -------------------------------------------------------------------------- !!
!! AMR

integer, parameter :: amr_each_ntimesteps = 3

!! Jeans criterion
integer, parameter :: ncells_refine = 32
integer, parameter :: ncells_coarse = ncells_refine

integer, parameter :: prof_report_rate = 128

end module
