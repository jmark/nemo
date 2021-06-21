!# define PP_SPLIT_FORM 1

module two_point_flux_mod

use constants_mod

private

public :: xflux, yflux, zflux

contains

pure function xflux(uL,uR) result(fstar)

    use equations_const_mod
    use setup_config_mod, only: kappa
    use setup_config_mod, only: mu0
    use equations_share_mod, only: GLM_ch

    real(dp), intent(in)    :: uL(N_VARS)
    real(dp), intent(in)    :: uR(N_VARS)
    real(dp)                :: fstar(N_VARS)

    real(dp), parameter :: kappaM1  = kappa - 1.0_dp
    real(dp), parameter :: skappaM1 = 1.0_dp/(kappa - 1.0_dp)
    real(dp), parameter :: smu_0    = 1.0_dp/mu0

    real(dp) :: betaLN,beta_R,beta_L
    real(dp) :: rhoLN,B2_L,B2_R,v2_L,v2_R
    real(dp) :: pTilde,p_L,p_R
    real(dp) :: v_L(3),v_R(3)
    real(dp) :: BAvg(3),vAvg(3)
    real(dp) :: v1_B2Avg
    real(dp) :: vB_Avg
    real(dp) :: psiAvg

    associate(&
        rho_L   => UL(1),   rho_R   => UR(1), &
        rhoV_L  => UL(2:4), rhoV_R  => UR(2:4), &
        E_L     => UL(5)-0.5_dp*smu_0*UL(9)**2, E_R => UR(5)-0.5_dp*smu_0*UR(9)**2, &
        psi_L   => UL(9), psi_R => UR(9), &
        B_L     => UL(6:8), B_R => UR(6:8))

    !! Get the inverse density, velocity, and pressure on left and right.
    v_L = rhoV_L(:)/rho_L
    v_R = rhoV_R(:)/rho_R

    v2_L = SUM(v_L(:)*v_L(:))
    v2_R = SUM(v_R(:)*v_R(:))
    B2_L = SUM(B_L(:)*B_L(:))
    B2_R = SUM(B_R(:)*B_R(:))

    p_L    = kappaM1*(E_L - 0.5_dp*(rho_L*v2_L+smu_0*B2_L))
    p_R    = kappaM1*(E_R - 0.5_dp*(rho_R*v2_R+smu_0*B2_R))
    beta_L = 0.5_dp*rho_L/p_L
    beta_R = 0.5_dp*rho_R/P_R

    !! Get the averages for the numerical flux.
    rhoLN      = ln_mean(rho_L, rho_R)
    betaLN     = ln_mean(beta_L,beta_R)
    vAvg       = 0.5_dp * (v_L + v_R)
    BAvg       = 0.5_dp * (B_L + B_R)
    v1_B2Avg   = 0.5_dp * (v_L(1)*B2_L + v_R(1)*B2_R)
    vB_Avg     = 0.5_dp * (SUM(V_L(:)*B_L(:)) + SUM(V_R(:)*B_R(:)))
    pTilde     = 0.5_dp * ((rho_L+rho_R)/(beta_L+beta_R)+smu_0*0.5_dp*(B2_L+B2_R))
    psiAvg     = 0.5_dp * (psi_L+psi_R)

    !! Entropy conserving and kinetic energy conserving flux.
    Fstar(1) = rhoLN*vAvg(1)
    Fstar(2) = Fstar(1)*vAvg(1) - smu_0*BAvg(1)*BAvg(1) + pTilde
    Fstar(3) = Fstar(1)*vAvg(2) - smu_0*BAvg(1)*BAvg(2)
    Fstar(4) = Fstar(1)*vAvg(3) - smu_0*BAvg(1)*BAvg(3)
    Fstar(7) = vAvg(1)*Bavg(2) - BAvg(1)*vAvg(2)
    Fstar(8) = vAvg(1)*Bavg(3) - BAvg(1)*vAvg(3)
    Fstar(6) = GLM_ch*psiAvg
    Fstar(9) = GLM_ch*BAvg(1)

    !! Total energy flux.
    Fstar(5) = Fstar(1)*0.5_dp*(skappaM1/betaLN - 0.5_dp*(v2_L+v2_R))  &
               + SUM(vAvg(:)*Fstar(2:4)) &
               + smu_0*(SUM(BAvg(:)*Fstar(6:8)) &
                       - 0.5_dp*v1_B2Avg + BAvg(1)*vB_Avg &
                       + Fstar(9)*psiAvg-GLM_ch*0.5_dp*(psi_L*B_L(1)+psi_R*B_R(1)))
    end associate 

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

!==================================================================================================================================
!> Computes the logarithmic mean: (aR-aL)/(LOG(aR)-LOG(aL)) = (aR-aL)/LOG(aR/aL)
!> Problem: if aL~= aR, then 0/0, but should tend to --> 0.5*(aR+aL)
!>
!> introduce xi=aR/aL and f=(aR-aL)/(aR+aL) = (xi-1)/(xi+1) 
!> => xi=(1+f)/(1-f) 
!> => Log(xi) = log(1+f)-log(1-f), and for smaRl f (f^2<1.0E-02) :
!>
!>    Log(xi) ~=     (f - 1/2 f^2 + 1/3 f^3 - 1/4 f^4 + 1/5 f^5 - 1/6 f^6 + 1/7 f^7)
!>                  +(f + 1/2 f^2 + 1/3 f^3 + 1/4 f^4 + 1/5 f^5 + 1/6 f^6 + 1/7 f^7)
!>             = 2*f*(1           + 1/3 f^2           + 1/5 f^4           + 1/7 f^6)
!>  (aR-aL)/Log(xi) = (aR+aL)*f/(2*f*(1 + 1/3 f^2 + 1/5 f^4 + 1/7 f^6)) = (aR+aL)/(2 + 2/3 f^2 + 2/5 f^4 + 2/7 f^6)
!>  (aR-aL)/Log(xi) = 0.5*(aR+aL)*(105/ (105+35 f^2+ 21 f^4 + 15 f^6)
!==================================================================================================================================
pure elemental function LN_MEAN(aL,aR)

    REAL(dp), INTENT(IN) :: aL  !< left value
    REAL(dp), INTENT(IN) :: aR  !< right value

    REAL(dp)            :: LN_MEAN  !< result

    REAL(dp)            :: Xi,u
    REAL(dp), PARAMETER :: eps = 1.0e-4_dp  ! tolerance for f^2, such that switch is smooth in double precision 

    Xi = aR/aL
    u = (Xi*(Xi-2.0_dp)+1.0_dp)/(Xi*(Xi+2.0_dp)+1.0_dp) !u=f^2, f=(aR-aL)/(aR+aL)=(xi-1)/(xi+1)
    LN_MEAN = MERGE((aL+aR)*52.5_dp/(105_dp + u*(35_dp + u*(21_dp +u*15_dp))), & !u <eps (test true)
                  (aR-aL)/LOG(Xi)                                         , & !u>=eps (test false)
                  (u.LT.eps)                                              )   !test

end function

end module
