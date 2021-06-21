module timedisc_const_mod

use constants_mod

integer, parameter :: nstages = 1

real(dp), parameter :: dt_coeffs(nstages) = (/1.0_dp/)

end module
