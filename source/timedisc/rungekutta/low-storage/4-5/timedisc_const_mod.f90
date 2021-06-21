module timedisc_const_mod

use constants_mod

real(dp), parameter :: a_1 = REAL(  0.0_dp                                , dp)
real(dp), parameter :: a_2 = REAL( -567301805773.0_dp /  1357537059087.0_dp,dp)
real(dp), parameter :: a_3 = REAL(-2404267990393.0_dp /  2016746695238.0_dp,dp)
real(dp), parameter :: a_4 = REAL(-3550918686646.0_dp /  2091501179385.0_dp,dp)
real(dp), parameter :: a_5 = REAL(-1275806237668.0_dp /   842570457699.0_dp,dp)

real(dp), parameter :: b_1 = REAL( 0.0_dp                                  ,dp)
real(dp), parameter :: b_2 = REAL( 1432997174477.0_dp /  9575080441755.0_dp,dp)
real(dp), parameter :: b_3 = REAL( 2526269341429.0_dp /  6820363962896.0_dp,dp)
real(dp), parameter :: b_4 = REAL( 2006345519317.0_dp /  3224310063776.0_dp,dp)
real(dp), parameter :: b_5 = REAL( 2802321613138.0_dp /  2924317926251.0_dp,dp)

real(dp), parameter :: c_1 = REAL( 1432997174477.0_dp /  9575080441755.0_dp,dp)
real(dp), parameter :: c_2 = REAL( 5161836677717.0_dp / 13612068292357.0_dp,dp)
real(dp), parameter :: c_3 = REAL( 1720146321549.0_dp /  2090206949498.0_dp,dp)
real(dp), parameter :: c_4 = REAL( 3134564353537.0_dp /  4481467310338.0_dp,dp)
real(dp), parameter :: c_5 = REAL( 2277821191437.0_dp / 14882151754819.0_dp,dp)

integer, parameter :: nstages = 5

real(dp), parameter :: a_rk(nstages) = (/ a_1,a_2,a_3,a_4,a_5 /)
real(dp), parameter :: b_rk(nstages) = (/ b_1,b_2,b_3,b_4,b_5 /)
real(dp), parameter :: c_rk(nstages) = (/ c_1,c_2,c_3,c_4,c_5 /)

real(dp), parameter :: dt_coeffs(nstages) = b_rk

end module
