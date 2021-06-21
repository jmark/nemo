!# define PP_SPLIT_FORM 1

module two_point_flux_mod

use constants_mod

private

public :: xflux, yflux, zflux

contains

pure elemental function logmean(a,b) result(r)

    real(dp), intent(in)    :: a,b
    real(dp)                :: x,u,r

    real(dp),parameter  :: eps = 1e-4_dp
    real(dp), parameter :: c1 = 1.0_dp/6.0_dp
    real(dp), parameter :: c2 = 2.0_dp/45.0_dp
    real(dp), parameter :: c3 = 22.0_dp/945.0_dp

    x = a/b

    if (abs(a-b) < eps) then
        u  = (x*(x-2.0_dp)+1.0_dp)/(x*(x+2.0_dp)+1.0_dp)
        r = (a+b)*(0.5_dp - u*(c1 - u*(c2 - c3*u)))
    else
        r = (a-b)/log(x)
    endif

end function

pure function xflux(u,v) result(f)

    use equations_const_mod
    use setup_config_mod, only: kappa
    use equations_mod, only: pressure

    real(dp), intent(in)    :: u(N_VARS)
    real(dp), intent(in)    :: v(N_VARS)
    real(dp)                :: f(N_VARS)

    real(dp)                :: dens,pres,velx,vely,velz,enth
    real(dp)                :: bu,bv !! beta
    real(dp)                :: pu,pv !! pressure
    real(dp)                :: uu(3),vv(3) !! velocity

    pu = pressure(u)
    pv = pressure(v)

    bu = 0.5_dp*u(DENS_VAR)/pu
    bv = 0.5_dp*v(DENS_VAR)/pv

    uu(1) = u(MOMX_VAR)/u(DENS_VAR)
    uu(2) = u(MOMY_VAR)/u(DENS_VAR)
    uu(3) = u(MOMZ_VAR)/u(DENS_VAR)

    vv(1) = v(MOMX_VAR)/v(DENS_VAR)
    vv(2) = v(MOMY_VAR)/v(DENS_VAR)
    vv(3) = v(MOMZ_VAR)/v(DENS_VAR)

    dens = logmean(u(DENS_VAR),v(DENS_VAR))
    velx = 0.5_dp*(uu(1) + vv(1))
    vely = 0.5_dp*(uu(2) + vv(2))
    velz = 0.5_dp*(uu(3) + vv(3))
    pres = 0.5_dp*(u(DENS_VAR)+v(DENS_VAR))/(bu+bv)

    enth = 0.5_dp/(kappa-1.0_dp)/logmean(bu,bv) + pres/dens + velx*velx + vely*vely + velz*velz - 0.25_dp*sum(uu*uu + vv*vv)

    f(DENS_VAR) = dens * velx 
    f(MOMX_VAR) = f(DENS_VAR) * velx + pres
    f(MOMY_VAR) = f(DENS_VAR) * vely 
    f(MOMZ_VAR) = f(DENS_VAR) * velz
    f(ENER_VAR) = f(DENS_VAR) * enth 

end function

pure function yflux(uL,uR) result(f)

    use equations_const_mod

    real(dp), intent(in) :: uL(N_VARS), uR(N_VARS)
    real(dp)             :: f(N_VARS)

    real(dp) :: ruL(N_VARS), ruR(N_VARS), rf(N_VARS) !! rotated

    ruL(DENS_VAR) = uL(DENS_VAR)
    ruL(MOMX_VAR) = uL(MOMY_VAR)
    ruL(MOMY_VAR) = uL(MOMZ_VAR)
    ruL(MOMZ_VAR) = uL(MOMX_VAR)
    ruL(ENER_VAR) = uL(ENER_VAR)
    
    ruR(DENS_VAR) = uR(DENS_VAR)
    ruR(MOMX_VAR) = uR(MOMY_VAR)
    ruR(MOMY_VAR) = uR(MOMZ_VAR)
    ruR(MOMZ_VAR) = uR(MOMX_VAR)
    ruR(ENER_VAR) = uR(ENER_VAR)

    rf = xflux(ruL,ruR)

    f(DENS_VAR) = rf(DENS_VAR)
    f(MOMX_VAR) = rf(MOMZ_VAR)
    f(MOMY_VAR) = rf(MOMX_VAR)
    f(MOMZ_VAR) = rf(MOMY_VAR)
    f(ENER_VAR) = rf(ENER_VAR)
    
end function

pure function zflux(uL,uR) result(f)

    use equations_const_mod

    real(dp), intent(in) :: uL(N_VARS), uR(N_VARS)
    real(dp)             :: f(N_VARS)

    real(dp) :: ruL(N_VARS), ruR(N_VARS), rf(N_VARS) !! rotated

    ruL(DENS_VAR) = uL(DENS_VAR)
    ruL(MOMX_VAR) = uL(MOMZ_VAR)
    ruL(MOMY_VAR) = uL(MOMX_VAR)
    ruL(MOMZ_VAR) = uL(MOMY_VAR)
    ruL(ENER_VAR) = uL(ENER_VAR)
    
    ruR(DENS_VAR) = uR(DENS_VAR)
    ruR(MOMX_VAR) = uR(MOMZ_VAR)
    ruR(MOMY_VAR) = uR(MOMX_VAR)
    ruR(MOMZ_VAR) = uR(MOMY_VAR)
    ruR(ENER_VAR) = uR(ENER_VAR)

    rf = xflux(ruL,ruR)

    f(DENS_VAR) = rf(DENS_VAR)
    f(MOMX_VAR) = rf(MOMY_VAR)
    f(MOMY_VAR) = rf(MOMZ_VAR)
    f(MOMZ_VAR) = rf(MOMX_VAR)
    f(ENER_VAR) = rf(ENER_VAR)
    
end function

end module
