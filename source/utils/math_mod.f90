module math_mod

use constants_mod

interface minmod
    module procedure minmod2
    module procedure minmod3_foo
    module procedure minmod4
    module procedure minmod5
end interface

contains

pure elemental function sigmoid(x) result(r)

    real(dp), intent(in) :: x
    real(dp)             :: r

    !r = 0.5_dp + 0.5_dp*tanh(0.5_dp*x)
    r = 1.0_dp/(1.0_dp + exp(-x))

end function

pure elemental function minmod2(a,b) result(r)

    real(dp), intent(in) :: a,b
    real(dp)             :: r

    r = 0.5_dp*(sign(1.0_dp,a)+sign(1.0_dp,b))*min(abs(a),abs(b))

end function

pure elemental function minmod3_foo(a,b,c) result(r)

    real(dp), intent(in) :: a,b,c
    real(dp)             :: r

    r = minmod2(minmod2(a,b),c)

end function

pure elemental function minmod4(a,b,c,d) result(r)

    real(dp), intent(in) :: a,b,c,d
    real(dp)             :: r

    r = minmod2(minmod2(minmod2(a,b),c),d)

end function

pure elemental function minmod5(a,b,c,d,e) result(r)

    real(dp), intent(in) :: a,b,c,d,e
    real(dp)             :: r

    r = minmod2(minmod2(minmod2(minmod2(a,b),c),d),e)

end function

pure elemental function median(l,c,r) result(m)
    
    real(dp), intent(in) :: l,c,r
    real(dp)             :: m

    m = c + minmod2(c-l,c-r)
    
end function

pure elemental function curvature(l,c,r) result(m)
    
    real(dp), intent(in) :: l,c,r
    real(dp)             :: m

    m = l - 2.0_dp*c + r
    
end function

pure elemental function minmod3(a1,a2,a3) result(mm)

    real(dp), intent(in) :: a1, a2, a3
    real(dp)             :: mm

    if (((sign(1.0_dp,a1)*sign(1.0_dp,a2)) > 0.0_dp) .and. ((sign(1.0_dp,a2)*sign(1.0_dp,a3)) > 0.0_dp)) then
        mm = sign(min(abs(a1),abs(a2),abs(a3)),a1)
    else
        mm = 0.0_dp
    end if

end function

end module
