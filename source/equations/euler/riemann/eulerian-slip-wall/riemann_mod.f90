module riemann_mod

use constants_mod

contains

!! Eulerian Slippery Wall
pure function eulerian_slip_wall(inner,outer) result(F)

    use equations_const_mod
    use equations_mod, only: pressure
    use setup_config_mod, only: kappa

    real(dp), intent(in) :: inner(N_VARS)
    real(dp), intent(in) :: outer(N_VARS)
    real(dp)             :: F(N_VARS)

    real(dp) :: v,p,c, a,b
    
    v = inner(MOMX_VAR)/inner(DENS_VAR)
    p = pressure(inner)

    if (v <= 0.0_dp) then

        c = SQRT(kappa*p/inner(DENS_VAR))
        p = p*max((1.0_dp + 0.5_dp*(kappa-1.0_dp)*v/c),1e-5_dp)**(2.0_dp*kappa/(kappa-1.0_dp))

    else

        a = 2.0_dp/(kappa+1.0_dp)/inner(DENS_VAR)
        b = (kappa-1.0_dp)/(kappa+1.0_dp)*p
        p = p + v/a*0.5_dp*(v + SQRT(v**2 + 4.0_dp*a*(p+b)))

    end if

    F = 0.0_dp
    F(MOMX_VAR) = p

end function

#:if PP_MODULE == 'riemann_north_mod'
pure function xriemann(U_L,U_R) result(F)

    use equations_const_mod

    real(dp), intent(in) :: U_L(N_VARS)
    real(dp), intent(in) :: U_R(N_VARS)
    real(dp)             :: F(N_VARS)

    real(dp) :: inner(N_VARS), outer(N_VARS)

    inner(DENS_VAR) =  U_R(DENS_VAR)
    inner(MOMX_VAR) = -U_R(MOMX_VAR)
    inner(MOMY_VAR) =  U_R(MOMY_VAR)
    inner(MOMZ_VAR) =  U_R(MOMZ_VAR)
    inner(ENER_VAR) =  U_R(ENER_VAR)
    
    outer(DENS_VAR) =  U_L(DENS_VAR)
    outer(MOMX_VAR) = -U_L(MOMX_VAR)
    outer(MOMY_VAR) =  U_L(MOMY_VAR)
    inner(MOMZ_VAR) =  U_R(MOMZ_VAR)
    outer(ENER_VAR) =  U_L(ENER_VAR)

    f = eulerian_slip_wall(inner,outer)

end function
#:endif

#:if PP_MODULE == 'riemann_south_mod'
pure function xriemann(U_L,U_R) result(F)

    real(dp), intent(in) :: U_L(N_VARS)
    real(dp), intent(in) :: U_R(N_VARS)
    real(dp)             :: F(N_VARS)

    f = eulerian_slip_wall(U_L,U_R)

end function
#:endif

#:if PP_MODULE == 'riemann_west_mod'
pure function yriemann(U_L,U_R) result(F)

    use equations_const_mod

    real(dp), intent(in) :: U_L(N_VARS)
    real(dp), intent(in) :: U_R(N_VARS)
    real(dp)             :: F(N_VARS)

    real(dp) :: inner(N_VARS), outer(N_VARS), rflux(N_VARS)

    inner(DENS_VAR) =  U_R(DENS_VAR)
    inner(MOMX_VAR) = -U_R(MOMY_VAR)
    inner(MOMY_VAR) =  U_R(MOMZ_VAR)
    inner(MOMZ_VAR) =  U_R(MOMX_VAR)
    inner(ENER_VAR) =  U_R(ENER_VAR)
    
    outer(DENS_VAR) =  U_L(DENS_VAR)
    outer(MOMX_VAR) = -U_L(MOMY_VAR)
    outer(MOMY_VAR) =  U_L(MOMZ_VAR)
    outer(MOMZ_VAR) =  U_L(MOMX_VAR)
    outer(ENER_VAR) =  U_L(ENER_VAR)

    rflux = eulerian_slip_wall(inner,outer)

    f(DENS_VAR) = rflux(DENS_VAR)
    f(MOMX_VAR) = rflux(MOMZ_VAR)
    f(MOMY_VAR) = rflux(MOMX_VAR)
    f(MOMZ_VAR) = rflux(MOMY_VAR)
    f(ENER_VAR) = rflux(ENER_VAR)

 end function
#:endif

#:if PP_MODULE == 'riemann_east_mod'
pure function yriemann(U_L,U_R) result(F)

    use equations_const_mod

    real(dp), intent(in) :: U_L(N_VARS)
    real(dp), intent(in) :: U_R(N_VARS)
    real(dp)             :: F(N_VARS)

    real(dp) :: inner(N_VARS), outer(N_VARS), rflux(N_VARS)

    inner(DENS_VAR) = U_L(DENS_VAR)
    inner(MOMX_VAR) = U_L(MOMY_VAR)
    inner(MOMY_VAR) = U_L(MOMZ_VAR)
    inner(MOMZ_VAR) = U_L(MOMX_VAR)
    inner(ENER_VAR) = U_L(ENER_VAR)
    
    outer(DENS_VAR) = U_R(DENS_VAR)
    outer(MOMX_VAR) = U_R(MOMY_VAR)
    outer(MOMY_VAR) = U_R(MOMZ_VAR)
    outer(MOMZ_VAR) = U_R(MOMX_VAR)
    outer(ENER_VAR) = U_R(ENER_VAR)

    rflux = eulerian_slip_wall(inner,outer)

    f(DENS_VAR) = rflux(DENS_VAR)
    f(MOMX_VAR) = rflux(MOMZ_VAR)
    f(MOMY_VAR) = rflux(MOMX_VAR)
    f(MOMZ_VAR) = rflux(MOMY_VAR)
    f(ENER_VAR) = rflux(ENER_VAR)

end function
#:endif

#:if PP_MODULE == 'riemann_front_mod'
pure function zriemann(U_L,U_R) result(F)

    use equations_const_mod

    real(dp), intent(in) :: U_L(N_VARS)
    real(dp), intent(in) :: U_R(N_VARS)
    real(dp)             :: F(N_VARS)

    real(dp) :: inner(N_VARS), outer(N_VARS), rflux(N_VARS)

    inner(DENS_VAR) =  U_R(DENS_VAR)
    inner(MOMX_VAR) = -U_R(MOMZ_VAR)
    inner(MOMY_VAR) =  U_R(MOMX_VAR)
    inner(MOMZ_VAR) =  U_R(MOMY_VAR)
    inner(ENER_VAR) =  U_R(ENER_VAR)
    
    outer(DENS_VAR) =  U_L(DENS_VAR)
    outer(MOMX_VAR) = -U_L(MOMZ_VAR)
    outer(MOMY_VAR) =  U_L(MOMX_VAR)
    outer(MOMZ_VAR) =  U_L(MOMY_VAR)
    outer(ENER_VAR) =  U_L(ENER_VAR)

    rflux = eulerian_slip_wall(inner,outer)

    f(DENS_VAR) = rflux(DENS_VAR)
    f(MOMX_VAR) = rflux(MOMY_VAR)
    f(MOMY_VAR) = rflux(MOMZ_VAR)
    f(MOMZ_VAR) = rflux(MOMX_VAR)
    f(ENER_VAR) = rflux(ENER_VAR)

 end function
#:endif

#:if PP_MODULE == 'riemann_back_mod'
pure function zriemann(U_L,U_R) result(F)

    use equations_const_mod

    real(dp), intent(in) :: U_L(N_VARS)
    real(dp), intent(in) :: U_R(N_VARS)
    real(dp)             :: F(N_VARS)

    real(dp) :: inner(N_VARS), outer(N_VARS), rflux(N_VARS)

    inner(DENS_VAR) = U_L(DENS_VAR)
    inner(MOMX_VAR) = U_L(MOMZ_VAR)
    inner(MOMY_VAR) = U_L(MOMX_VAR)
    inner(MOMZ_VAR) = U_L(MOMY_VAR)
    inner(ENER_VAR) = U_L(ENER_VAR)
    
    outer(DENS_VAR) = U_R(DENS_VAR)
    outer(MOMX_VAR) = U_R(MOMZ_VAR)
    outer(MOMY_VAR) = U_R(MOMX_VAR)
    outer(MOMZ_VAR) = U_R(MOMY_VAR)
    outer(ENER_VAR) = U_R(ENER_VAR)

    rflux = eulerian_slip_wall(inner,outer)

    f(DENS_VAR) = rflux(DENS_VAR)
    f(MOMX_VAR) = rflux(MOMY_VAR)
    f(MOMY_VAR) = rflux(MOMZ_VAR)
    f(MOMZ_VAR) = rflux(MOMX_VAR)
    f(ENER_VAR) = rflux(ENER_VAR)

end function
#:endif

end module
