!# define PP_N_VARS 9
!# define PP_EQN_MHD 1

module equations_const_mod

!! Conservative variables.
integer, parameter :: DENS_VAR = 1

integer, parameter :: MOMX_VAR = 2
integer, parameter :: MOMY_VAR = 3
integer, parameter :: MOMZ_VAR = 4

integer, parameter :: ENER_VAR = 5

integer, parameter :: MAGX_VAR = 6
integer, parameter :: MAGY_VAR = 7
integer, parameter :: MAGZ_VAR = 8

integer, parameter :: GLMP_VAR = 9

!! Primitive variables.
integer, parameter :: VELX_VAR = 2
integer, parameter :: VELY_VAR = 3
integer, parameter :: VELZ_VAR = 4

integer, parameter :: PRES_VAR = 5

!! Positive variables.
integer, parameter :: posvars(2) = (/DENS_VAR,PRES_VAR/)

end module
