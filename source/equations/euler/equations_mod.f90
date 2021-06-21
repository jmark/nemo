module equations_mod

use constants_mod
use equations_const_mod

interface cons2entr
    procedure cons2entr_elemental
    procedure cons2entr_0d
    procedure cons2entr_1d
    procedure cons2entr_2d
    procedure cons2entr_3d
end interface

interface soundspeed
    procedure soundspeed_elemental
    procedure soundspeed_0d
    procedure soundspeed_1d
    procedure soundspeed_2d
    procedure soundspeed_3d
end interface

interface pressure
    procedure pressure_elemental
    procedure pressure_0d
    procedure pressure_1d
    procedure pressure_2d
    procedure pressure_3d
end interface

interface xlambdaMax
    procedure xlambdaMax_elemental
    procedure xlambdaMax_0d
    procedure xlambdaMax_1d
    procedure xlambdaMax_2d
    procedure xlambdaMax_3d
end interface

interface ylambdaMax
    procedure ylambdaMax_elemental
    procedure ylambdaMax_0d
    procedure ylambdaMax_1d
    procedure ylambdaMax_2d
    procedure ylambdaMax_3d
end interface

interface zlambdaMax
    procedure zlambdaMax_elemental
    procedure zlambdaMax_0d
    procedure zlambdaMax_1d
    procedure zlambdaMax_2d
    procedure zlambdaMax_3d
end interface

interface isvalid
    procedure isvalid_elemental
    procedure isvalid_0d
    procedure isvalid_1d
    procedure isvalid_2d
    procedure isvalid_3d
end interface

contains

pure function xflux(u) result(f)

    real(dp), intent(in)    :: u(N_VARS)
    real(dp)                :: f(N_VARS)
    real(dp)                :: pres,velx

    pres = pressure(u)
    velx = u(MOMX_VAR)/u(DENS_VAR)

    f(DENS_VAR) = u(MOMX_VAR)
    f(MOMX_VAR) = velx*u(MOMX_VAR) + pres
    f(MOMY_VAR) = velx*u(MOMY_VAR)
    f(MOMZ_VAR) = velx*u(MOMZ_VAR)
    f(ENER_VAR) = velx*(u(ENER_VAR) + pres)

end function

pure function yflux(u) result(f)

    real(dp), intent(in)    :: u(N_VARS)
    real(dp)                :: f(N_VARS)
    real(dp)                :: ru(N_VARS), rf(N_VARS) !! rotated

    !! rotating to x-direction
    ru(DENS_VAR) = u(DENS_VAR)
    ru(MOMX_VAR) = u(MOMY_VAR)
    ru(MOMY_VAR) = u(MOMZ_VAR)
    ru(MOMZ_VAR) = u(MOMX_VAR)
    ru(ENER_VAR) = u(ENER_VAR)
    
    rf = xflux(ru)

    !! rotating back
    f(DENS_VAR) = rf(DENS_VAR)
    f(MOMX_VAR) = rf(MOMZ_VAR)
    f(MOMY_VAR) = rf(MOMX_VAR)
    f(MOMZ_VAR) = rf(MOMY_VAR)
    f(ENER_VAR) = rf(ENER_VAR)
    
end function

pure function zflux(u) result(f)

    real(dp), intent(in)    :: u(N_VARS)
    real(dp)                :: f(N_VARS)
    real(dp)                :: ru(N_VARS), rf(N_VARS) !! rotated

    !! rotating to x-direction
    ru(DENS_VAR) = u(DENS_VAR)
    ru(MOMX_VAR) = u(MOMZ_VAR)
    ru(MOMY_VAR) = u(MOMX_VAR)
    ru(MOMZ_VAR) = u(MOMY_VAR)
    ru(ENER_VAR) = u(ENER_VAR)
    
    rf = xflux(ru)

    !! rotating back
    f(DENS_VAR) = rf(DENS_VAR)
    f(MOMX_VAR) = rf(MOMY_VAR)
    f(MOMY_VAR) = rf(MOMZ_VAR)
    f(MOMZ_VAR) = rf(MOMX_VAR)
    f(ENER_VAR) = rf(ENER_VAR)
    
end function

!! ========================================================================== !! 

pure elemental function pressure_elemental(dens, momx, momy, momz, ener) result(pres)

    use setup_config_mod, only: kappa

    real(dp),intent(in)    :: dens, momx, momy, momz, ener 
    real(dp)               :: pres

    pres = (kappa-1.0_dp) * (ener - 0.5_dp/dens * (momx*momx + momy*momy + momz*momz))

end function

pure function pressure_0d(u) result(res)

    use setup_config_mod, only: kappa

    real(dp),intent(in)    :: u(N_VARS)
    real(dp)               :: res

    res = pressure(u(DENS_VAR),u(MOMX_VAR),u(MOMY_VAR),u(MOMZ_VAR),u(ENER_VAR))

end function

pure function pressure_1d(u) result(p)

    real(dp),intent(in)    :: u(:,:)
    real(dp)               :: p(size(u,dim=1))

    p = pressure(u(:,DENS_VAR),u(:,MOMX_VAR),u(:,MOMY_VAR),u(:,MOMZ_VAR),u(:,ENER_VAR))

end function

pure function pressure_2d(u) result(p)

    real(dp),intent(in)    :: u(:,:,:)
    real(dp)               :: p(size(u,dim=1),size(u,dim=2))

    p = pressure(u(:,:,DENS_VAR),u(:,:,MOMX_VAR),u(:,:,MOMY_VAR),u(:,:,MOMZ_VAR),u(:,:,ENER_VAR))

end function

pure function pressure_3d(u) result(p)

    real(dp),intent(in)    :: u(:,:,:,:)
    real(dp)               :: p(size(u,dim=1),size(u,dim=2),size(u,dim=3))

    p = pressure(u(:,:,:,DENS_VAR),u(:,:,:,MOMX_VAR),u(:,:,:,MOMY_VAR),u(:,:,:,MOMZ_VAR),u(:,:,:,ENER_VAR))

end function

!! ========================================================================== !! 

pure function cons2prim(u) result(p)

    real(dp),intent(in)    :: u(N_VARS)
    real(dp)               :: p(N_VARS)

    p(DENS_VAR) = u(DENS_VAR)
    p(VELX_VAR) = u(MOMX_VAR)/u(DENS_VAR)
    p(VELY_VAR) = u(MOMY_VAR)/u(DENS_VAR)
    p(VELZ_VAR) = u(MOMZ_VAR)/u(DENS_VAR)
    p(PRES_VAR) = pressure(u)

end function

pure function prim2cons(p) result(u)
    
    use setup_config_mod, only: kappa

    real(dp),intent(in)    :: p(N_VARS)
    real(dp)               :: u(N_VARS)

    u(DENS_VAR) = p(DENS_VAR)
    u(MOMX_VAR) = p(VELX_VAR)*p(DENS_VAR)
    u(MOMY_VAR) = p(VELY_VAR)*p(DENS_VAR)
    u(MOMZ_VAR) = p(VELZ_VAR)*p(DENS_VAR)
    u(ENER_VAR) = p(PRES_VAR)/(kappa-1.0_dp) + 0.5_dp*p(DENS_VAR)*(p(VELX_VAR)**2 + p(VELY_VAR)**2 + p(VELZ_VAR)**2)

end function

!! ========================================================================== !! 

pure elemental function cons2entr_elemental(dens,momx,momy,momz,ener) result(s)

    use setup_config_mod, only: kappa

    real(dp),intent(in)    :: dens,momx,momy,momz,ener
    real(dp)               :: s

    !! entropy density
    s = log(pressure(dens,momx,momy,momz,ener)) - kappa*log(dens)
    s = -dens*s/(kappa-1.0_dp)

end function

pure function cons2entr_0d(u) result(s)

    real(dp),intent(in)    :: u(N_VARS)
    real(dp)               :: s

    s = cons2entr(u(DENS_VAR),u(MOMX_VAR),u(MOMY_VAR),u(MOMZ_VAR),u(ENER_VAR))

end function

pure function cons2entr_1d(u) result(s)

    real(dp),intent(in)    :: u(N_NODES,N_VARS)
    real(dp)               :: s(N_NODES)

    s = cons2entr(u(:,DENS_VAR),u(:,MOMX_VAR),u(:,MOMY_VAR),u(:,MOMZ_VAR),u(:,ENER_VAR))

end function

pure function cons2entr_2d(u) result(s)

    real(dp),intent(in)    :: u(N_NODES,N_NODES,N_VARS)
    real(dp)               :: s(N_NODES,N_NODES)

    s = cons2entr(u(:,:,DENS_VAR),u(:,:,MOMX_VAR),u(:,:,MOMY_VAR),u(:,:,MOMZ_VAR),u(:,:,ENER_VAR))

end function

pure function cons2entr_3d(u) result(s)

    real(dp),intent(in)    :: u(N_NODES,N_NODES,N_NODES,N_VARS)
    real(dp)               :: s(N_NODES,N_NODES,N_NODES)

    s = cons2entr(u(:,:,:,DENS_VAR),u(:,:,:,MOMX_VAR),u(:,:,:,MOMY_VAR),u(:,:,:,MOMZ_VAR),u(:,:,:,ENER_VAR))

end function

pure function cons2evec(u) result(w)

    use setup_config_mod, only: kappa

    real(dp),intent(in)    :: u(N_VARS)
    real(dp)               :: w(N_VARS)

    real(dp) :: p(N_VARS),s,b

    !! primitive variables
    p = cons2prim(u)

    !! compressible gas dynamics entropy
    s = log(p(PRES_VAR)) - kappa*log(p(DENS_VAR))

    !! beta
    b = p(DENS_VAR)/p(PRES_VAR)
        
    !! entropy variables
    w(DENS_VAR) = (kappa - s)/(kappa-1.0_dp) - 0.5_dp*b*(p(VELX_VAR)**2 + p(VELY_VAR)**2 + p(VELZ_VAR)**2)
    w(VELX_VAR) = b*p(VELX_VAR)
    w(VELY_VAR) = b*p(VELY_VAR)
    w(VELZ_VAR) = b*p(VELZ_VAR)
    w(PRES_VAR) = -b

end function

pure elemental function soundspeed_elemental(dens,momx,momy,momz,ener) result(c)

    use setup_config_mod, only: kappa

    real(dp),intent(in)    :: dens,momx,momy,momz,ener 
    real(dp)               :: c

    c = sqrt(kappa*pressure(dens,momx,momy,momz,ener)/dens)

end function

pure function soundspeed_0d(u) result(c)

    real(dp),intent(in)    :: u(N_VARS)
    real(dp)               :: c

    c = soundspeed(u(DENS_VAR),u(MOMX_VAR),u(MOMY_VAR),u(MOMZ_VAR),u(ENER_VAR))

end function

pure function soundspeed_1d(u) result(c)

    real(dp),intent(in)    :: u(:,:)
    real(dp)               :: c(size(u,dim=1))

    c = soundspeed(u(:,DENS_VAR),u(:,MOMX_VAR),u(:,MOMY_VAR),u(:,MOMZ_VAR),u(:,ENER_VAR))

end function

pure function soundspeed_2d(u) result(c)

    real(dp),intent(in)    :: u(:,:,:)
    real(dp)               :: c(size(u,dim=1),size(u,dim=2))

    c = soundspeed(u(:,:,DENS_VAR),u(:,:,MOMX_VAR),u(:,:,MOMY_VAR),u(:,:,MOMZ_VAR),u(:,:,ENER_VAR))

end function

pure function soundspeed_3d(u) result(c)

    real(dp),intent(in)    :: u(:,:,:,:)
    real(dp)               :: c(size(u,dim=1),size(u,dim=2),size(u,dim=3))

    c = soundspeed(u(:,:,:,DENS_VAR),u(:,:,:,MOMX_VAR),u(:,:,:,MOMY_VAR),u(:,:,:,MOMZ_VAR),u(:,:,:,ENER_VAR))

end function

!! ========================================================================== !! 

pure elemental function isvalid_elemental(dens,momx,momy,momz,ener) result(rest)

    real(dp),intent(in)     :: dens,momx,momy,momz,ener 
    logical                 :: rest

    rest = (dens > 0.0_dp) .and. (pressure(dens,momx,momy,momz,ener) > 0.0_dp)

end function

pure function isvalid_0d(u) result(ok)

    real(dp),intent(in) :: u(N_VARS)
    logical             :: ok

    ok = isvalid(u(DENS_VAR),u(MOMX_VAR),u(MOMY_VAR),u(MOMZ_VAR),u(ENER_VAR))

end function

pure function isvalid_1d(u) result(ok)

    real(dp),intent(in)     :: u(:,:)
    logical                 :: ok

    ok = all(isvalid(u(:,DENS_VAR),u(:,MOMX_VAR),u(:,MOMY_VAR),u(:,MOMZ_VAR),u(:,ENER_VAR)))

end function

pure function isvalid_2d(u) result(ok)

    real(dp),intent(in)     :: u(:,:,:)
    logical                 :: ok

    ok = all(isvalid(u(:,:,DENS_VAR),u(:,:,MOMX_VAR),u(:,:,MOMY_VAR),u(:,:,MOMZ_VAR),u(:,:,ENER_VAR)))

end function

pure function isvalid_3d(u) result(ok)

    real(dp),intent(in)     :: u(:,:,:,:)
    logical                 :: ok

    ok = all(isvalid(u(:,:,:,DENS_VAR),u(:,:,:,MOMX_VAR),u(:,:,:,MOMY_VAR),u(:,:,:,MOMZ_VAR),u(:,:,:,ENER_VAR)))

end function

!! ========================================================================== !! 

pure elemental function xlambdaMax_elemental(dens,momx,momy,momz,ener) result(res)

    use setup_config_mod, only: kappa

    real(dp),intent(in)    :: dens,momx,momy,momz,ener
    real(dp)               :: res 

    res = abs(momx/dens) + soundspeed(dens,momx,momy,momz,ener)

end function

pure function xlambdaMax_0d(u) result(res)

    use setup_config_mod, only: kappa

    real(dp),intent(in)    :: u(N_VARS)
    real(dp)               :: res 

    res = xlambdaMax(u(DENS_VAR),u(MOMX_VAR),u(MOMY_VAR),u(MOMZ_VAR),u(ENER_VAR))

end function

pure function xlambdaMax_1d(u) result(l)

    real(dp),intent(in)    :: u(:,:)
    real(dp)               :: l(size(u,dim=1))

    l = xlambdaMax(u(:,DENS_VAR),u(:,MOMX_VAR),u(:,MOMY_VAR),u(:,MOMZ_VAR),u(:,ENER_VAR))

end function

pure function xlambdaMax_2d(u) result(l)

    real(dp),intent(in)    :: u(:,:,:)
    real(dp)               :: l(size(u,dim=1),size(u,dim=2))

    l = xlambdaMax(u(:,:,DENS_VAR),u(:,:,MOMX_VAR),u(:,:,MOMY_VAR),u(:,:,MOMZ_VAR),u(:,:,ENER_VAR))

end function

pure function xlambdaMax_3d(u) result(l)

    real(dp),intent(in)    :: u(:,:,:,:)
    real(dp)               :: l(size(u,dim=1),size(u,dim=2),size(u,dim=3))

    l = xlambdaMax(u(:,:,:,DENS_VAR),u(:,:,:,MOMX_VAR),u(:,:,:,MOMY_VAR),u(:,:,:,MOMZ_VAR),u(:,:,:,ENER_VAR))

end function

!! ========================================================================== !! 

pure elemental function ylambdaMax_elemental(dens,momx,momy,momz,ener) result(res)

    real(dp),intent(in)    :: dens, momx, momy, momz, ener
    real(dp)               :: res

    res = xlambdaMax(dens,momy,momz,momx,ener)

end function

pure function ylambdaMax_0d(u) result(res)

    use setup_config_mod, only: kappa

    real(dp),intent(in)    :: u(N_VARS)
    real(dp)               :: res 

    res = ylambdaMax(u(DENS_VAR),u(MOMX_VAR),u(MOMY_VAR),u(MOMZ_VAR),u(ENER_VAR))

end function

pure function ylambdaMax_1d(u) result(l)

    real(dp),intent(in)    :: u(:,:)
    real(dp)               :: l(size(u,dim=1))

    l = ylambdaMax(u(:,DENS_VAR),u(:,MOMX_VAR),u(:,MOMY_VAR),u(:,MOMZ_VAR),u(:,ENER_VAR))

end function


pure function ylambdaMax_2d(u) result(l)

    real(dp),intent(in)    :: u(:,:,:)
    real(dp)               :: l(size(u,dim=1),size(u,dim=2))

    l = ylambdaMax(u(:,:,DENS_VAR),u(:,:,MOMX_VAR),u(:,:,MOMY_VAR),u(:,:,MOMZ_VAR),u(:,:,ENER_VAR))

end function

pure function ylambdaMax_3d(u) result(l)

    real(dp),intent(in)    :: u(:,:,:,:)
    real(dp)               :: l(size(u,dim=1),size(u,dim=2),size(u,dim=3))

    l = ylambdaMax(u(:,:,:,DENS_VAR),u(:,:,:,MOMX_VAR),u(:,:,:,MOMY_VAR),u(:,:,:,MOMZ_VAR),u(:,:,:,ENER_VAR))

end function

!! ========================================================================== !! 

pure elemental function zlambdaMax_elemental(dens, momx, momy, momz, ener) result(res)

    real(dp),intent(in)    :: dens, momx, momy, momz, ener
    real(dp)               :: res

    res = xlambdaMax(dens,momz,momx,momy,ener)

end function

pure function zlambdaMax_0d(u) result(res)

    use setup_config_mod, only: kappa

    real(dp),intent(in)    :: u(N_VARS)
    real(dp)               :: res 

    res = zlambdaMax(u(DENS_VAR),u(MOMX_VAR),u(MOMY_VAR),u(MOMZ_VAR),u(ENER_VAR))

end function

pure function zlambdaMax_1d(u) result(l)

    real(dp),intent(in)    :: u(:,:)
    real(dp)               :: l(size(u,dim=1))

    l = zlambdaMax(u(:,DENS_VAR),u(:,MOMX_VAR),u(:,MOMY_VAR),u(:,MOMZ_VAR),u(:,ENER_VAR))

end function

pure function zlambdaMax_2d(u) result(l)

    real(dp),intent(in)    :: u(:,:,:)
    real(dp)               :: l(size(u,dim=1),size(u,dim=2))

    l = zlambdaMax(u(:,:,DENS_VAR),u(:,:,MOMX_VAR),u(:,:,MOMY_VAR),u(:,:,MOMZ_VAR),u(:,:,ENER_VAR))

end function

pure function zlambdaMax_3d(u) result(l)

    real(dp),intent(in)    :: u(:,:,:,:)
    real(dp)               :: l(size(u,dim=1),size(u,dim=2),size(u,dim=3))

    l = zlambdaMax(u(:,:,:,DENS_VAR),u(:,:,:,MOMX_VAR),u(:,:,:,MOMY_VAR),u(:,:,:,MOMZ_VAR),u(:,:,:,ENER_VAR))

end function

end module
