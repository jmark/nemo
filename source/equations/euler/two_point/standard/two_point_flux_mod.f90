!# define PP_SPLIT_FORM 1

module two_point_flux_mod

use constants_mod

private

public :: xflux, yflux, zflux

contains

pure function xflux(u,v) result(r)

    use equations_mod, only: f => xflux

    real(dp), intent(in)    :: u(N_VARS)
    real(dp), intent(in)    :: v(N_VARS)
    real(dp)                :: r(N_VARS)

    r = 0.5_dp*(f(u) + f(v))

end function

pure function yflux(u,v) result(r)

    use equations_mod, only: g => yflux

    real(dp), intent(in)    :: u(N_VARS)
    real(dp), intent(in)    :: v(N_VARS)

    real(dp)                :: r(N_VARS)

    r = 0.5_dp*(g(u) + g(v))

end function

pure function zflux(u,v) result(r)

    use equations_mod, only: h => zflux

    real(dp), intent(in)    :: u(N_VARS)
    real(dp), intent(in)    :: v(N_VARS)

    real(dp)                :: r(N_VARS)

    r = 0.5_dp*(h(u) + h(v))

end function

end module
