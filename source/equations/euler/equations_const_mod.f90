!# define PP_N_VARS 5

module equations_const_mod

use constants_mod

integer, parameter :: DENS_VAR = 1
integer, parameter :: MOMX_VAR = 2
integer, parameter :: MOMY_VAR = 3
integer, parameter :: MOMZ_VAR = 4
integer, parameter :: ENER_VAR = 5

integer, parameter :: VELX_VAR = 2
integer, parameter :: VELY_VAR = 3
integer, parameter :: VELZ_VAR = 4
integer, parameter :: PRES_VAR = 5

integer, parameter :: posvars(2) = (/DENS_VAR,PRES_VAR/)

real(dp), parameter :: lowfactor = 1e-3_dp

end module
