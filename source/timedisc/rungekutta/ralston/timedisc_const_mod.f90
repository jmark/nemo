module timedisc_const_mod

use constants_mod

integer, parameter :: nstages = 2

real(dp), parameter :: dt_coeffs(nstages) = (/ 0.0_dp, 2.0_dp/3.0_dp /)

end module
