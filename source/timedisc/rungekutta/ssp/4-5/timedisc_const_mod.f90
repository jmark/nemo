module timedisc_const_mod

use constants_mod

real(dp), parameter :: A0 = 1.00000000000000_dp
real(dp), parameter :: A1 = 1.00000000000000_dp
real(dp), parameter :: A2 = 0.00000000000000_dp
real(dp), parameter :: A3 = 0.39175222700392_dp

real(dp), parameter :: B0 = 0.39175222700392_dp
real(dp), parameter :: B1 = 0.44437049406734_dp
real(dp), parameter :: B2 = 0.55562950593266_dp
real(dp), parameter :: B3 = 0.36841059262959_dp

real(dp), parameter :: C0 = 0.58607968896780_dp
real(dp), parameter :: C1 = 0.62010185138540_dp
real(dp), parameter :: C2 = 0.37989814861460_dp
real(dp), parameter :: C3 = 0.25189177424738_dp

real(dp), parameter :: D0 = 0.47454236302687_dp
real(dp), parameter :: D1 = 0.17807995410773_dp
real(dp), parameter :: D2 = 0.82192004589227_dp
real(dp), parameter :: D3 = 0.54497475021237_dp

real(dp), parameter :: E0 = 0.93501063100924_dp
real(dp), parameter :: E1 = 0.00683325884039_dp
real(dp), parameter :: E2 = 0.34833675773694_dp
real(dp), parameter :: E3 = 0.22600748319395_dp

real(dp), parameter :: E4 = 0.51723167208978_dp
real(dp), parameter :: E5 = 0.12759831133288_dp
real(dp), parameter :: E6 = 0.08460416338212_dp

integer, parameter :: nstages = 5

real(dp), parameter :: dt_coeffs(nstages) = (/A0, B0, C0, D0, E0/)

end module
