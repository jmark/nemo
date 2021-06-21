module timedisc_const_mod

use constants_mod

integer, parameter :: nstages = 3

real(dp), parameter :: dt_coeffs(nstages) = (/0.0_dp, 1.0_dp, 0.5_dp/)

end module
