module equations_mod

use constants_mod
use equations_const_mod

! interface cons2entr
!     procedure cons2entr_elemental
!     procedure cons2entr_0d
!     procedure cons2entr_1d
!     procedure cons2entr_2d
!     procedure cons2entr_3d
! end interface

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

    use setup_config_mod, only: kappa,mu0
    use equations_share_mod, only: glm_ch

    real(dp), intent(in)    :: u(N_VARS)
    real(dp)                :: f(N_VARS)
    real(dp)                :: pres,etot,vels(3)

    real(dp), parameter :: smu_0 = 1.0_dp/mu0
    real(dp), parameter :: s2mu_0 = 0.5_dp/mu0

    real(dp), parameter :: kappaM1 = kappa-1.0_dp
    real(dp), parameter :: kappaM2 = kappa-2.0_dp

    real(dp) :: srho,v1,v2,v3,v_2,bb2,vb,pt,Ep

    associate(rho   => U(1), &
            rhov1   => U(2), &
            rhov2   => U(3), &
            rhov3   => U(4), &
            Etotal  => U(5) - 0.5*smu_0*U(9)**2, &
            b1      => U(6), &
            b2      => U(7), &
            b3      => U(8))

    srho = 1.0_dp / rho
    v1   = rhov1*srho 
    v2   = rhov2*srho 
    v3   = rhov3*srho 
    v_2  = v1*v1+v2*v2+v3*v3 
    bb2  = (b1*b1+b2*b2+b3*b3)
    vb   = (b1*v1+b2*v2+b3*v3)

    pt   = kappaM1*(Etotal-0.5_dp*rho*(v_2)) - KappaM2*s2mu_0*bb2
    Ep   = (Etotal + pt)

    f(1) = rhov1                   
    f(2) = rhov1*v1+pt - smu_0*b1*b1
    f(3) = rhov1*v2    - smu_0*b1*b2
    f(4) = rhov1*v3    - smu_0*b1*b3
    f(5) = Ep*v1       - smu_0*b1*vb

    f(6) = 0.0_dp
    f(7) = v1*b2-b1*v2
    f(8) = v1*b3-b1*v3

    f(5) = f(5)+smu_0*GLM_ch*b1*U(9)
    f(6) = f(6)+GLM_ch*U(9)
    f(9) = GLM_ch*b1
    end associate

    !! There are some bugs in the following code...

    ! pres = pressure(u)
    ! vels = u(MOMX_VAR:MOMZ_VAR)/u(DENS_VAR)

    ! f(DENS_VAR) = u(MOMX_VAR)

    ! f(MOMX_VAR) = pres + 0.5_dp/mu0*sum(u(MAGX_VAR:MAGZ_VAR)**2) &
    !             + vels(X_DIR)*u(MOMX_VAR) - 1.0_dp/mu0*u(MAGX_VAR)*u(MAGX_VAR)
    ! f(MOMY_VAR) = vels(X_DIR)*u(MOMY_VAR) - 1.0_dp/mu0*u(MAGX_VAR)*u(MAGY_VAR)
    ! f(MOMZ_VAR) = vels(X_DIR)*u(MOMZ_VAR) - 1.0_dp/mu0*u(MAGX_VAR)*u(MAGZ_VAR)

    ! f(ENER_VAR) = vels(X_DIR)*(0.5_dp*u(DENS_VAR)*sum(vels**2) &
    !     & + 1.0_dp/mu0*sum(u(MAGX_VAR:MAGZ_VAR)**2) + kappa/(kappa-1.0_dp)*pres) &
    !     & - 1.0_dp/mu0*u(MAGX_VAR)*sum(vels*u(MAGX_VAR:MAGZ_VAR)) &
    !     & + GLM_ch/mu0*u(MAGX_VAR)*u(GLMP_VAR)

    ! f(MAGX_VAR) = GLM_ch*u(GLMP_VAR)
    ! f(MAGY_VAR) = u(VELX_VAR)*u(MAGY_VAR) - u(VELY_VAR)*u(MAGX_VAR)
    ! f(MAGZ_VAR) = u(VELX_VAR)*u(MAGZ_VAR) - u(VELZ_VAR)*u(MAGX_VAR)

    ! f(GLMP_VAR) = GLM_ch*u(MAGX_VAR)

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
    
    ru(MAGX_VAR) = u(MAGY_VAR)
    ru(MAGY_VAR) = u(MAGZ_VAR)
    ru(MAGZ_VAR) = u(MAGX_VAR)

    ru(GLMP_VAR) = u(GLMP_VAR)

    rf = xflux(ru)

    !! rotating back
    f(DENS_VAR) = rf(DENS_VAR)

    f(MOMX_VAR) = rf(MOMZ_VAR)
    f(MOMY_VAR) = rf(MOMX_VAR)
    f(MOMZ_VAR) = rf(MOMY_VAR)

    f(ENER_VAR) = rf(ENER_VAR)
    
    f(MAGX_VAR) = rf(MAGZ_VAR)
    f(MAGY_VAR) = rf(MAGX_VAR)
    f(MAGZ_VAR) = rf(MAGY_VAR)

    f(GLMP_VAR) = rf(GLMP_VAR)

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

    ru(MAGX_VAR) = u(MAGZ_VAR)
    ru(MAGY_VAR) = u(MAGX_VAR)
    ru(MAGZ_VAR) = u(MAGY_VAR)
    
    ru(GLMP_VAR) = u(GLMP_VAR)

    rf = xflux(ru)

    !! rotating back
    f(DENS_VAR) = rf(DENS_VAR)

    f(MOMX_VAR) = rf(MOMY_VAR)
    f(MOMY_VAR) = rf(MOMZ_VAR)
    f(MOMZ_VAR) = rf(MOMX_VAR)

    f(ENER_VAR) = rf(ENER_VAR)

    f(MAGX_VAR) = rf(MAGY_VAR)
    f(MAGY_VAR) = rf(MAGZ_VAR)
    f(MAGZ_VAR) = rf(MAGX_VAR)
    
    f(GLMP_VAR) = rf(GLMP_VAR)

end function

!! ========================================================================== !! 

pure elemental function pressure_elemental(dens,momx,momy,momz,ener,magx,magy,magz,glmp) result(p)

    use setup_config_mod, only: kappa,mu0

    real(dp),intent(in)    :: dens, momx, momy, momz, magx,magy,magz, ener, glmp
    real(dp)               :: p

    p = (kappa-1.0_dp) * (ener &
        - 0.5_dp/dens * (momx*momx + momy*momy + momz*momz) &
        - 0.5_dp/mu0  * (magx*magx + magy*magy + magz*magz + glmp*glmp))

end function

pure function pressure_0d(u) result(p)

    use setup_config_mod, only: kappa

    real(dp),intent(in)    :: u(N_VARS)
    real(dp)               :: p

    p = pressure(u(DENS_VAR),u(MOMX_VAR),u(MOMY_VAR),u(MOMZ_VAR),u(ENER_VAR),u(MAGX_VAR),u(MAGY_VAR),u(MAGZ_VAR),u(GLMP_VAR))

end function

pure function pressure_1d(u) result(p)

    real(dp),intent(in)    :: u(:,:)
    real(dp)               :: p(size(u,dim=1))

    p = pressure(u(:,DENS_VAR),u(:,MOMX_VAR),u(:,MOMY_VAR),u(:,MOMZ_VAR),u(:,ENER_VAR),u(:,MAGX_VAR),u(:,MAGY_VAR),u(:,MAGZ_VAR),u(:,GLMP_VAR))

end function

pure function pressure_2d(u) result(p)

    real(dp),intent(in)    :: u(:,:,:)
    real(dp)               :: p(size(u,dim=1),size(u,dim=2))

    p = pressure(u(:,:,DENS_VAR),u(:,:,MOMX_VAR),u(:,:,MOMY_VAR),u(:,:,MOMZ_VAR),u(:,:,ENER_VAR),u(:,:,MAGX_VAR),u(:,:,MAGY_VAR),u(:,:,MAGZ_VAR),u(:,:,GLMP_VAR))

end function

pure function pressure_3d(u) result(p)

    real(dp),intent(in)    :: u(:,:,:,:)
    real(dp)               :: p(size(u,dim=1),size(u,dim=2),size(u,dim=3))

    p = pressure(u(:,:,:,DENS_VAR),u(:,:,:,MOMX_VAR),u(:,:,:,MOMY_VAR),u(:,:,:,MOMZ_VAR),u(:,:,:,ENER_VAR),u(:,:,:,MAGX_VAR),u(:,:,:,MAGY_VAR),u(:,:,:,MAGZ_VAR),u(:,:,:,GLMP_VAR))

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

    p(MAGX_VAR:MAGZ_VAR) = u(MAGX_VAR:MAGZ_VAR)
    p(GLMP_VAR) = u(GLMP_VAR)

end function

pure function prim2cons(p) result(u)
    
    use setup_config_mod, only: kappa,mu0

    real(dp),intent(in)    :: p(N_VARS)
    real(dp)               :: u(N_VARS)

    u(DENS_VAR) = p(DENS_VAR)

    u(MOMX_VAR:MOMZ_VAR) = p(DENS_VAR)*p(VELX_VAR:VELZ_VAR)

    u(ENER_VAR) = p(PRES_VAR)/(kappa-1.0_dp) &
        + 0.5_dp*p(DENS_VAR)*sum(p(VELX_VAR:VELZ_VAR)**2) &
        + 0.5_dp/mu0*(sum(p(MAGX_VAR:MAGZ_VAR)**2) + p(GLMP_VAR)**2)

    u(MAGX_VAR:MAGZ_VAR) = p(MAGX_VAR:MAGZ_VAR)
    u(GLMP_VAR) = p(GLMP_VAR)

end function

!! ========================================================================== !! 

!! pure elemental function cons2entr_elemental(dens,momx,momy,momz,ener) result(s)
!! 
!!     use setup_config_mod, only: kappa
!! 
!!     real(dp),intent(in)    :: dens,momx,momy,momz,ener
!!     real(dp)               :: s
!! 
!!     !! entropy density
!!     s = log(pressure(dens,momx,momy,momz,ener)) - kappa*log(dens)
!!     s = -dens*s/(kappa-1.0_dp)
!! 
!! end function
!! 
!! pure function cons2entr_0d(u) result(s)
!! 
!!     real(dp),intent(in)    :: u(N_VARS)
!!     real(dp)               :: s
!! 
!!     s = cons2entr(u(DENS_VAR),u(MOMX_VAR),u(MOMY_VAR),u(MOMZ_VAR),u(ENER_VAR))
!! 
!! end function
!! 
!! pure function cons2entr_1d(u) result(s)
!! 
!!     real(dp),intent(in)    :: u(N_NODES,N_VARS)
!!     real(dp)               :: s(N_NODES)
!! 
!!     s = cons2entr(u(:,DENS_VAR),u(:,MOMX_VAR),u(:,MOMY_VAR),u(:,MOMZ_VAR),u(:,ENER_VAR))
!! 
!! end function
!! 
!! pure function cons2entr_2d(u) result(s)
!! 
!!     real(dp),intent(in)    :: u(N_NODES,N_NODES,N_VARS)
!!     real(dp)               :: s(N_NODES,N_NODES)
!! 
!!     s = cons2entr(u(:,:,DENS_VAR),u(:,:,MOMX_VAR),u(:,:,MOMY_VAR),u(:,:,MOMZ_VAR),u(:,:,ENER_VAR))
!! 
!! end function
!! 
!! pure function cons2entr_3d(u) result(s)
!! 
!!     real(dp),intent(in)    :: u(N_NODES,N_NODES,N_NODES,N_VARS)
!!     real(dp)               :: s(N_NODES,N_NODES,N_NODES)
!! 
!!     s = cons2entr(u(:,:,:,DENS_VAR),u(:,:,:,MOMX_VAR),u(:,:,:,MOMY_VAR),u(:,:,:,MOMZ_VAR),u(:,:,:,ENER_VAR))
!! 
!! end function
!! 
!! pure function cons2evec(u) result(w)
!! 
!!     use setup_config_mod, only: kappa
!! 
!!     real(dp),intent(in)    :: u(N_VARS)
!!     real(dp)               :: w(N_VARS)
!! 
!!     real(dp) :: p(N_VARS),s,b
!! 
!!     !! primitive variables
!!     p = cons2prim(u)
!! 
!!     !! compressible gas dynamics entropy
!!     s = log(p(PRES_VAR)) - kappa*log(p(DENS_VAR))
!! 
!!     !! beta
!!     b = p(DENS_VAR)/p(PRES_VAR)
!!         
!!     !! entropy variables
!!     w(DENS_VAR) = (kappa - s)/(kappa-1.0_dp) - 0.5_dp*b*(p(VELX_VAR)**2 + p(VELY_VAR)**2 + p(VELZ_VAR)**2)
!!     w(VELX_VAR) = b*p(VELX_VAR)
!!     w(VELY_VAR) = b*p(VELY_VAR)
!!     w(VELZ_VAR) = b*p(VELZ_VAR)
!!     w(PRES_VAR) = -b
!! 
!! end function

pure elemental function soundspeed_elemental(dens,momx,momy,momz,ener,magx,magy,magz,glmp) result(c)

    use setup_config_mod, only: kappa

    real(dp),intent(in)    :: dens,momx,momy,momz,ener,magx,magy,magz,glmp
    real(dp)               :: c

    real(dp) :: u(N_VARS)

    U(DENS_VAR) = dens

    U(MOMX_VAR) = momx
    U(MOMY_VAR) = momy
    U(MOMZ_VAR) = momz

    U(ENER_VAR) = ener

    U(MAGX_VAR) = magx
    U(MAGY_VAR) = magy
    U(MAGZ_VAR) = magz

    U(GLMP_VAR) = glmp

    call Fastest_Signal_Speed(u,c)

    ! c = sqrt(kappa*pressure(dens,momx,momy,momz,ener,magx,magy,magz,glmp)/dens)

end function

pure function soundspeed_0d(u) result(c)

    real(dp),intent(in)    :: u(N_VARS)
    real(dp)               :: c

    c = soundspeed(u(DENS_VAR),&
        u(MOMX_VAR),u(MOMY_VAR),u(MOMZ_VAR),u(ENER_VAR),&
        u(MAGX_VAR),u(MAGY_VAR),u(MAGZ_VAR),u(GLMP_VAR))

end function

pure function soundspeed_1d(u) result(c)

    real(dp),intent(in)    :: u(:,:)
    real(dp)               :: c(size(u,dim=1))

    c = soundspeed(u(:,DENS_VAR),&
        u(:,MOMX_VAR),u(:,MOMY_VAR),u(:,MOMZ_VAR),u(:,ENER_VAR),&
        u(:,MAGX_VAR),u(:,MAGY_VAR),u(:,MAGZ_VAR),u(:,GLMP_VAR))

end function

pure function soundspeed_2d(u) result(c)

    real(dp),intent(in)    :: u(:,:,:)
    real(dp)               :: c(size(u,dim=1),size(u,dim=2))

    c = soundspeed(u(:,:,DENS_VAR),&
        u(:,:,MOMX_VAR),u(:,:,MOMY_VAR),u(:,:,MOMZ_VAR),u(:,:,ENER_VAR),&
        u(:,:,MAGX_VAR),u(:,:,MAGY_VAR),u(:,:,MAGZ_VAR),u(:,:,GLMP_VAR))

end function

pure function soundspeed_3d(u) result(c)

    real(dp),intent(in)    :: u(:,:,:,:)
    real(dp)               :: c(size(u,dim=1),size(u,dim=2),size(u,dim=3))

    c = soundspeed(u(:,:,:,DENS_VAR),&
        u(:,:,:,MOMX_VAR),u(:,:,:,MOMY_VAR),u(:,:,:,MOMZ_VAR),u(:,:,:,ENER_VAR),&
        u(:,:,:,MAGX_VAR),u(:,:,:,MAGY_VAR),u(:,:,:,MAGZ_VAR),u(:,:,:,GLMP_VAR))

end function

!! ========================================================================== !! 

pure elemental function isvalid_elemental(dens,momx,momy,momz,ener,magx,magy,magz,glmp) result(rest)

    real(dp),intent(in)     :: dens,momx,momy,momz,ener,magx,magy,magz,glmp
    logical                 :: rest

    rest = (dens > 0.0_dp) .and. (pressure(dens,momx,momy,momz,ener,magx,magy,magz,glmp) > 0.0_dp)

end function

pure function isvalid_0d(u) result(ok)

    real(dp),intent(in) :: u(N_VARS)
    logical             :: ok

    ok = isvalid(u(DENS_VAR),&
        u(MOMX_VAR),u(MOMY_VAR),u(MOMZ_VAR),u(ENER_VAR),&
        u(MAGX_VAR),u(MAGY_VAR),u(MAGZ_VAR),u(GLMP_VAR))

end function

pure function isvalid_1d(u) result(ok)

    real(dp),intent(in)     :: u(:,:)
    logical                 :: ok

    ok = all(isvalid(u(:,DENS_VAR),&
        u(:,MOMX_VAR),u(:,MOMY_VAR),u(:,MOMZ_VAR),u(:,ENER_VAR),&
        u(:,MAGX_VAR),u(:,MAGY_VAR),u(:,MAGZ_VAR),u(:,GLMP_VAR)))

end function

pure function isvalid_2d(u) result(ok)

    real(dp),intent(in)     :: u(:,:,:)
    logical                 :: ok

    ok = all(isvalid(u(:,:,DENS_VAR),&
        u(:,:,MOMX_VAR),u(:,:,MOMY_VAR),u(:,:,MOMZ_VAR),u(:,:,ENER_VAR),&
        u(:,:,MAGX_VAR),u(:,:,MAGY_VAR),u(:,:,MAGZ_VAR),u(:,:,GLMP_VAR)))

end function

pure function isvalid_3d(u) result(ok)

    real(dp),intent(in)     :: u(:,:,:,:)
    logical                 :: ok

    ok = all(isvalid(u(:,:,:,DENS_VAR),&
        u(:,:,:,MOMX_VAR),u(:,:,:,MOMY_VAR),u(:,:,:,MOMZ_VAR),u(:,:,:,ENER_VAR),&
        u(:,:,:,MAGX_VAR),u(:,:,:,MAGY_VAR),u(:,:,:,MAGZ_VAR),u(:,:,:,GLMP_VAR)))

end function

!! ========================================================================== !! 

pure elemental function xlambdaMax_elemental(dens,momx,momy,momz,ener,magx,magy,magz,glmp) result(res)

    use setup_config_mod, only: kappa

    real(dp),intent(in)    :: dens,momx,momy,momz,ener,magx,magy,magz,glmp
    real(dp)               :: res 

    res = abs(momx/dens) + soundspeed(dens,momx,momy,momz,ener,magx,magy,magz,glmp)

end function

pure function xlambdaMax_0d(u) result(res)

    use setup_config_mod, only: kappa

    real(dp),intent(in)    :: u(N_VARS)
    real(dp)               :: res 

    res = xlambdaMax(u(DENS_VAR),&
        u(MOMX_VAR),u(MOMY_VAR),u(MOMZ_VAR),u(ENER_VAR),&
        u(MAGX_VAR),u(MAGY_VAR),u(MAGZ_VAR),u(GLMP_VAR))

end function

pure function xlambdaMax_1d(u) result(l)

    real(dp),intent(in)    :: u(:,:)
    real(dp)               :: l(size(u,dim=1))

    l = xlambdaMax(u(:,DENS_VAR),&
        u(:,MOMX_VAR),u(:,MOMY_VAR),u(:,MOMZ_VAR),u(:,ENER_VAR),&
        u(:,MAGX_VAR),u(:,MAGY_VAR),u(:,MAGZ_VAR),u(:,GLMP_VAR))

end function

pure function xlambdaMax_2d(u) result(l)

    real(dp),intent(in)    :: u(:,:,:)
    real(dp)               :: l(size(u,dim=1),size(u,dim=2))

    l = xlambdaMax(u(:,:,DENS_VAR),&
        u(:,:,MOMX_VAR),u(:,:,MOMY_VAR),u(:,:,MOMZ_VAR),u(:,:,ENER_VAR),&
        u(:,:,MAGX_VAR),u(:,:,MAGY_VAR),u(:,:,MAGZ_VAR),u(:,:,GLMP_VAR))

end function

pure function xlambdaMax_3d(u) result(l)

    real(dp),intent(in)    :: u(:,:,:,:)
    real(dp)               :: l(size(u,dim=1),size(u,dim=2),size(u,dim=3))

    l = xlambdaMax(u(:,:,:,DENS_VAR),&
        u(:,:,:,MOMX_VAR),u(:,:,:,MOMY_VAR),u(:,:,:,MOMZ_VAR),u(:,:,:,ENER_VAR),&
        u(:,:,:,MAGX_VAR),u(:,:,:,MAGY_VAR),u(:,:,:,MAGZ_VAR),u(:,:,:,GLMP_VAR))

end function

!! ========================================================================== !! 

pure elemental function ylambdaMax_elemental(dens,momx,momy,momz,ener,magx,magy,magz,glmp) result(res)

    real(dp),intent(in)    :: dens,momx,momy,momz,ener,magx,magy,magz,glmp
    real(dp)               :: res

    res = xlambdaMax(dens,momy,momz,momx,ener,magy,magz,magx,glmp)

end function

pure function ylambdaMax_0d(u) result(res)

    use setup_config_mod, only: kappa

    real(dp),intent(in)    :: u(N_VARS)
    real(dp)               :: res 

    res = ylambdaMax(u(DENS_VAR),&
            u(MOMX_VAR),u(MOMY_VAR),u(MOMZ_VAR),u(ENER_VAR),&
            u(MAGX_VAR),u(MAGY_VAR),u(MAGZ_VAR),u(GLMP_VAR))

end function

pure function ylambdaMax_1d(u) result(l)

    real(dp),intent(in)    :: u(:,:)
    real(dp)               :: l(size(u,dim=1))

    l = ylambdaMax(u(:,DENS_VAR),&
        u(:,MOMX_VAR),u(:,MOMY_VAR),u(:,MOMZ_VAR),u(:,ENER_VAR),&
        u(:,MAGX_VAR),u(:,MAGY_VAR),u(:,MAGZ_VAR),u(:,GLMP_VAR))

end function

pure function ylambdaMax_2d(u) result(l)

    real(dp),intent(in)    :: u(:,:,:)
    real(dp)               :: l(size(u,dim=1),size(u,dim=2))

    l = ylambdaMax(u(:,:,DENS_VAR),&
        u(:,:,MOMX_VAR),u(:,:,MOMY_VAR),u(:,:,MOMZ_VAR),u(:,:,ENER_VAR),&
        u(:,:,MAGX_VAR),u(:,:,MAGY_VAR),u(:,:,MAGZ_VAR),u(:,:,GLMP_VAR))

end function

pure function ylambdaMax_3d(u) result(l)

    real(dp),intent(in)    :: u(:,:,:,:)
    real(dp)               :: l(size(u,dim=1),size(u,dim=2),size(u,dim=3))

    l = ylambdaMax(u(:,:,:,DENS_VAR),&
        u(:,:,:,MOMX_VAR),u(:,:,:,MOMY_VAR),u(:,:,:,MOMZ_VAR),u(:,:,:,ENER_VAR),&
        u(:,:,:,MAGX_VAR),u(:,:,:,MAGY_VAR),u(:,:,:,MAGZ_VAR),u(:,:,:,GLMP_VAR))

end function

!! ========================================================================== !! 

pure elemental function zlambdaMax_elemental(dens,momx,momy,momz,ener,magx,magy,magz,glmp) result(res)

    real(dp),intent(in)    :: dens,momx,momy,momz,ener,magx,magy,magz,glmp
    real(dp)               :: res

    res = xlambdaMax(dens,momz,momx,momy,ener,magz,magx,magy,glmp)

end function

pure function zlambdaMax_0d(u) result(res)

    use setup_config_mod, only: kappa

    real(dp),intent(in)    :: u(N_VARS)
    real(dp)               :: res 

    res = zlambdaMax(u(DENS_VAR),&
        u(MOMX_VAR),u(MOMY_VAR),u(MOMZ_VAR),u(ENER_VAR),&
        u(MAGX_VAR),u(MAGY_VAR),u(MAGZ_VAR),u(GLMP_VAR))

end function

pure function zlambdaMax_1d(u) result(l)

    real(dp),intent(in)    :: u(:,:)
    real(dp)               :: l(size(u,dim=1))

    l = zlambdaMax(u(:,DENS_VAR),&
        u(:,MOMX_VAR),u(:,MOMY_VAR),u(:,MOMZ_VAR),u(:,ENER_VAR),&
        u(:,MAGX_VAR),u(:,MAGY_VAR),u(:,MAGZ_VAR),u(:,GLMP_VAR))

end function

pure function zlambdaMax_2d(u) result(l)

    real(dp),intent(in)    :: u(:,:,:)
    real(dp)               :: l(size(u,dim=1),size(u,dim=2))

    l = zlambdaMax(u(:,:,DENS_VAR),&
        u(:,:,MOMX_VAR),u(:,:,MOMY_VAR),u(:,:,MOMZ_VAR),u(:,:,ENER_VAR),&
        u(:,:,MAGX_VAR),u(:,:,MAGY_VAR),u(:,:,MAGZ_VAR),u(:,:,GLMP_VAR))

end function

pure function zlambdaMax_3d(u) result(l)

    real(dp),intent(in)    :: u(:,:,:,:)
    real(dp)               :: l(size(u,dim=1),size(u,dim=2),size(u,dim=3))

    l = zlambdaMax(u(:,:,:,DENS_VAR),&
        u(:,:,:,MOMX_VAR),u(:,:,:,MOMY_VAR),u(:,:,:,MOMZ_VAR),u(:,:,:,ENER_VAR),&
        u(:,:,:,MAGX_VAR),u(:,:,:,MAGY_VAR),u(:,:,:,MAGZ_VAR),u(:,:,:,GLMP_VAR))

end function

pure subroutine Fastest_Signal_Speed(U,c)

    use setup_config_mod, only: mu0,kappa

    real(dp), intent(in)    :: U(N_VARS)
    real(dp), intent(out)   :: c

    real(dp) :: cs2,va2,ca2,pres
    real(dp) :: astar

    pres = pressure(U)

    cs2 = kappa*pres/U(DENS_VAR)            !! sound wave speed
    ca2 = U(MAGX_VAR)**2/U(DENS_VAR)/mu0    !! Alfen wave speed

    va2 = SUM(U(MAGY_VAR:MAGZ_VAR)**2)/U(DENS_VAR)/mu0 + ca2

    astar = SQRT((cs2 + va2)**2 - 4.0_dp*cs2*ca2)

    c = SQRT(0.5_dp*(cs2 + va2 + astar))

end subroutine

pure subroutine Fastest_Signal_Speeds_Roe(U_L,U_R,RoeVel,RoeSnd)

    !! Paper by Cargo & Gallice: "Roe Matrices for Ideal MHD and ...",1997)

    use setup_config_mod, only: mu0,kappa

    real(dp), intent(in)    :: U_L(N_VARS)
    real(dp), intent(in)    :: U_R(N_VARS)

    real(dp), intent(out)   :: RoeVel  !! Roe averaged velocity
    real(dp), intent(out)   :: RoeSnd  !! Roe averaged MHD soundspeed

    real(dp) :: prims_L(N_VARS)
    real(dp) :: prims_R(N_VARS)

    prims_L = cons2prim(U_L)
    prims_R = cons2prim(U_R)

    call FastestWave1D_Roe(U_L,U_R,prims_L,prims_R,RoeVel,RoeSnd)

end subroutine

pure SUBROUTINE FastestWave1D_Roe(ConsL,ConsR,PrimL,PrimR,RoeVelx,cf_Roe)

    use setup_config_mod, only: mu0,kappa

    !----------------------------------------------------------------------------------------------------------------------------------
    ! INPUT VARIABLES
    REAL(dp),INTENT(IN)     :: ConsL(N_VARS) !< vector of left conserved variables, rotated!!!
    REAL(dp),INTENT(IN)     :: ConsR(N_VARS) !< vector of right conserved variables, rotated!!!
    REAL(dp),INTENT(IN)     :: primL(N_VARS) !< vector of left primitive variables, rotated!!!
    REAL(dp),INTENT(IN)     :: primR(N_VARS) !< vector of right primitive variables, rotated!!!
    !----------------------------------------------------------------------------------------------------------------------------------
    ! OUTPUT VARIABLES
    REAL(dp),INTENT(OUT)    :: RoeVelx        !< roe velocity 
    REAL(dp),INTENT(OUT)    :: cf_Roe         !< fastest wave
    !----------------------------------------------------------------------------------------------------------------------------------
    ! LOCAL VARIABLES 
    REAL(dp)             :: SqrtRho_L,SqrtRho_R,sSqrtRoe_LR
    REAL(dp)             :: Roe_L,Roe_R
    REAL(dp)             :: ptot_L,ptot_R
    REAL(dp)             :: sRho_Roe,RoeVel(3),RoeB(3)
    REAL(dp)             :: H_L,H_R,RoeH
    REAL(dp)             :: XX, va2_Roe,c2_Roe,ca2_Roe,astar_Roe

    real(dp), parameter :: smu_0 = 1.0_dp/mu0
    real(dp), parameter :: s2mu_0 = 0.5_dp/mu0

    !==================================================================================================================================

    SqrtRho_L   = SQRT(PrimL(1))
    SqrtRho_R   = SQRT(PrimR(1))
    sSqrtRoe_LR = 1./(SqrtRho_L+SqrtRho_R)
    Roe_L       = SqrtRho_L*sSqrtRoe_LR
    Roe_R       = SqrtRho_R*sSqrtRoe_LR
    ptot_L      = PrimL(5)+s2mu_0*SUM(PrimL(6:8)**2) !Total presssure!
    ptot_R      = PrimR(5)+s2mu_0*SUM(PrimR(6:8)**2) !Total presssure!

    sRho_Roe    = 1./(SqrtRho_L*SqrtRho_R)
    RoeVel(1:3) = Roe_L*PrimL(2:4)+Roe_R*PrimR(2:4)
    RoeB(1:3)   = Roe_L*PrimR(6:8)+Roe_R*PrimL(6:8)

    H_L  = (ConsL(5)+ptot_L)/PrimL(1)
    H_R  = (ConsR(5)+ptot_R)/PrimR(1)
    RoeH = Roe_L*H_L+Roe_R*H_R

    XX   = 0.5*SUM((PrimL(6:8)-PrimR(6:8))**2)*(sSqrtRoe_LR**2)

    va2_Roe   = SUM(RoeB(:)**2)*smu_0*sRho_Roe
    c2_Roe    = (2.-Kappa)*XX+ (Kappa-1)*(RoeH-0.5*SUM(RoeVel(:)**2)-va2_Roe)
    ca2_Roe   = RoeB(1)**2*smu_0*sRho_Roe
    astar_Roe = SQRT((c2_Roe+va2_Roe)**2-4.*c2_Roe*ca2_Roe)
    cf_Roe    = SQRT(0.5*(c2_Roe+va2_Roe+astar_Roe))
    RoeVelx   = RoeVel(1)

END SUBROUTINE

!! pure subroutine Fastest_Signal_Speeds_Roe(U_L,U_R,RoeVel,RoeSnd)
!! 
!!     !! Paper by Cargo & Gallice: "Roe Matrices for Ideal MHD and ...",1997)
!! 
!!     use setup_config_mod, only: mu0,kappa
!! 
!!     real(dp), intent(in)    :: U_L(N_VARS)
!!     real(dp), intent(in)    :: U_R(N_VARS)
!! 
!!     real(dp), intent(out)   :: RoeVel  !! Roe averaged velocity
!!     real(dp), intent(out)   :: RoeSnd  !! Roe averaged MHD soundspeed
!! 
!!     real(dp) :: SqrtRho_L,SqrtRho_R,sSqrtRoe_LR
!!     real(dp) :: Roe_L,Roe_R
!!     real(dp) :: ptot_L,ptot_R
!!     real(dp) :: sRho_Roe,RoeVels(3),RoeB(3)
!!     real(dp) :: H_L,H_R,RoeH
!!     real(dp) :: XX, va2_Roe,c2_Roe,ca2_Roe,astar_Roe
!! 
!!     real(dp) :: p_L
!!     real(dp) :: p_R
!! 
!!     real(dp) :: v_L(3)
!!     real(dp) :: v_R(3)
!! 
!!     real(dp), parameter :: smu_0 = 1.0_dp/mu0
!!     real(dp), parameter :: s2mu_0 = 0.5_dp/mu0
!!     
!!     v_L = U_L(MOMX_VAR:MOMZ_VAR)/U_L(DENS_VAR)
!!     v_R = U_R(MOMX_VAR:MOMZ_VAR)/U_R(DENS_VAR)
!! 
!!     p_L = pressure(U_L)
!!     p_R = pressure(U_R)
!! 
!!     SqrtRho_L   = SQRT(U_L(DENS_VAR))
!!     SqrtRho_R   = SQRT(U_R(DENS_VAR))
!! 
!!     sSqrtRoe_LR = 1.0_dp/(SqrtRho_L+SqrtRho_R)
!!     Roe_L       = SqrtRho_L*sSqrtRoe_LR
!!     Roe_R       = SqrtRho_R*sSqrtRoe_LR
!!     ptot_L      = p_L + s2mu_0*SUM(U_L(MAGX_VAR:MAGZ_VAR)**2) !Total presssure!
!!     ptot_R      = p_R + s2mu_0*SUM(U_R(MAGX_VAR:MAGZ_VAR)**2) !Total presssure!
!! 
!!     sRho_Roe    = 1.0_dp/(SqrtRho_L*SqrtRho_R)
!!     RoeVels     = Roe_L*v_L+Roe_R*v_R
!!     RoeB        = Roe_L*U_R(MAGX_VAR:MAGZ_VAR)+Roe_R*U_L(MAGX_VAR:MAGZ_VAR)
!! 
!!     H_L         = (U_L(ENER_VAR)+ptot_L)/U_L(DENS_VAR)
!!     H_R         = (U_R(ENER_VAR)+ptot_R)/U_R(DENS_VAR)
!!     RoeH        = Roe_L*H_L+Roe_R*H_R
!! 
!!     XX          = 0.5_dp*SUM((U_L(MAGX_VAR:MAGZ_VAR)-U_R(MAGX_VAR:MAGZ_VAR))**2)*(sSqrtRoe_LR**2)
!! 
!!     va2_Roe     = SUM(RoeB**2)*smu_0*sRho_Roe
!!     c2_Roe      = (2.0_dp-Kappa)*XX + (Kappa-1.0_dp)*(RoeH-0.5_dp*SUM(RoeVels**2)-va2_Roe)
!!     ca2_Roe     = RoeB(1)**2*smu_0*sRho_Roe
!!     astar_Roe   = SQRT((c2_Roe+va2_Roe)**2-4.0_dp*c2_Roe*ca2_Roe)
!! 
!!     !! output
!!     RoeSnd      = SQRT(0.5_dp*(c2_Roe+va2_Roe+astar_Roe))
!!     RoeVel      = RoeVels(1)
!! 
!! end subroutine

SUBROUTINE WaveSpeeds1D(Prim,ca,cs,cf)

    use setup_config_mod, only: kappa
    use setup_config_mod, only: mu0

    !----------------------------------------------------------------------------------------------------------------------------------
    ! INPUT VARIABLES
    REAL(dp),INTENT(IN)     :: prim(N_VARS) !< vector of primitive variables, rotated!!!
    !----------------------------------------------------------------------------------------------------------------------------------
    ! OUTPUT VARIABLES
    REAL(dp),INTENT(OUT)    :: ca,cs,cf      !< alfven, sound and fast magnetosonic wave speed
    !----------------------------------------------------------------------------------------------------------------------------------
    ! LOCAL VARIABLES 
    REAL(dp)                :: sRho,c2,va2,ca2
    REAL(dp)                :: astar

    real(dp), parameter :: smu_0 = 1.0_dp/mu0

    !==================================================================================================================================

    sRho=1./prim(1)
    c2=kappa*prim(5)*sRho
    ca2=Prim(6)*Prim(6)*sRho*smu_0 !Alfen wave speed
    va2=SUM(prim(7:8)*prim(7:8))*sRho*smu_0+ca2
    ca=SQRT(ca2)

    astar=SQRT((c2+va2)*(c2+va2)-4.*c2*ca2)
    cs=SQRT(0.5*(c2+va2-astar))
    cf=SQRT(0.5*(c2+va2+astar))

END SUBROUTINE

pure SUBROUTINE FastestWave1D(Prim,cf)

    USE setup_config_mod, ONLY: mu0
    USE setup_config_mod, ONLY: Kappa

    !----------------------------------------------------------------------------------------------------------------------------------
    ! INPUT VARIABLES
    REAL(dp),INTENT(IN)     :: prim(N_VARS) !< vector of primitive variables, rotated!!!
    !----------------------------------------------------------------------------------------------------------------------------------
    ! OUTPUT VARIABLES
    REAL(dp),INTENT(OUT)    :: cf            !< fast magnetosonic wave speed
    !----------------------------------------------------------------------------------------------------------------------------------
    ! LOCAL VARIABLES 
    REAL(dp)                :: sRho,c2,va2,ca2
    REAL(dp)                :: astar
    !==================================================================================================================================

    real(dp), parameter :: smu_0 = 1.0_dp/mu0

    sRho=1./prim(1)
    c2=kappa*prim(5)*sRho
    ca2=Prim(6)*Prim(6)*sRho*smu_0 !Alfen wave speed
    va2=SUM(prim(7:8)*prim(7:8))*sRho*smu_0+ca2

    astar=SQRT((c2+va2)*(c2+va2)-4.*c2*ca2)
    cf=SQRT(0.5*(c2+va2+astar))

END SUBROUTINE FastestWave1D

end module
