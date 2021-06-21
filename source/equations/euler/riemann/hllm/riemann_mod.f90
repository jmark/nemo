module riemann_mod

use constants_mod

contains

!! HLLM Riemann solver
pure function xriemann(U_L,U_R) result(F)

    use equations_const_mod
    use equations_mod, only: pressure
    use setup_config_mod, only: kappa

    use equations_mod, only: flux => xflux

    real(dp), intent(in) :: U_L(N_VARS), U_R(N_VARS)
    real(dp)             :: F(N_VARS)
    
    real(dp) :: F_L(N_VARS), F_R(N_VARS)

    real(dp) :: c_L,c_R !! sound speed
    real(dp) :: H_L,H_R !! enthalpy
    real(dp) :: p_L,p_R !! pressure
    real(dp) :: S_L,S_R !! signal wave speed
    real(dp) :: V_L(N_DIMS),V_R(N_DIMS) !! velocities

    real(dp) :: alpha(N_DIMS)
    real(dp) :: roeVel(N_DIMS),roeH,roeC,roeDens
    real(dp) :: rx(N_VARS),ry(N_VARS),rz(N_VARS)

    real(dp) :: sqrtrho_L,sqrtrho_R,ssqrtrho,absvel
    real(dp), parameter :: beta = SQRT(0.5_dp*(kappa-1.0_dp)/kappa)

    real(dp) :: smooth

    !! velocities
    V_L = U_L(MOMX_VAR:MOMZ_VAR)/U_L(DENS_VAR)
    V_R = U_R(MOMX_VAR:MOMZ_VAR)/U_R(DENS_VAR)

    !! pressures
    p_L = pressure(U_L)
    p_R = pressure(U_R)

    !! sound speeds
    c_L = SQRT(kappa*p_L/U_L(DENS_VAR))
    c_R = SQRT(kappa*p_R/U_R(DENS_VAR))

    !! enthalpy
    H_L = (U_L(ENER_VAR) + p_L)/U_L(DENS_VAR)
    H_R = (U_R(ENER_VAR) + p_R)/U_R(DENS_VAR)

    !! auxiliaries
    SqrtRho_L = SQRT(U_L(DENS_VAR))
    SqrtRho_R = SQRT(U_R(DENS_VAR))
    sSqrtRho  = 1.0_dp/(SqrtRho_L+SqrtRho_R)

    !! Roe averages
    RoeVel    = (SqrtRho_R*V_R + SqrtRho_L*V_L) * sSqrtRho
    RoeH      = (SqrtRho_R*H_R + SqrtRho_L*H_L) * sSqrtRho
    absVel    = DOT_PRODUCT(RoeVel,RoeVel)
    RoeC      = SQRT((kappa-1.0_dp)*(RoeH-0.5_dp*absVel))
    RoeDens   = SQRT(U_L(DENS_VAR)*U_R(DENS_VAR))

    !! signal velocities
    S_L = MIN(RoeVel(X_DIR) - RoeC, V_L(X_DIR) - beta*c_L, 0.0_dp)
    S_R = MAX(RoeVel(X_DIR) + RoeC, V_R(X_DIR) + beta*c_R, 0.0_dp)

    if (S_L .ge. 0.0_dp) then
        F = flux(U_L)

    else if (S_R .le. 0.0_dp) then
        F = flux(U_R)

    else
        !! central fluxes
        F_L = flux(U_L)
        F_R = flux(U_R)

        !! wave strenghts
        alpha(X_DIR) = (U_R(DENS_VAR) - U_L(DENS_VAR))  - (p_R-p_L)/(RoeC*RoeC)
        alpha(Y_DIR) = RoeDens*(V_R(Y_DIR) - V_L(Y_DIR))
        alpha(Z_DIR) = RoeDens*(V_R(Z_DIR) - V_L(Z_DIR))

        !! eigenvectors
        rx = (/ 1.0_dp, RoeVel(X_DIR),  RoeVel(Y_DIR),  RoeVel(Z_DIR),  0.5*absVel     /)
        ry = (/ 0.0_dp, 0.0_dp,         1.0_dp,         0.0_dp,         RoeVel(Y_DIR)  /)
        rz = (/ 0.0_dp, 0.0_dp,         0.0_dp,         1.0_dp,         RoeVel(Z_DIR)  /)

        F = (S_R*F_L - S_L*F_R + S_L*S_R * &
                (U_R - U_L - RoeC/(RoeC + ABS(0.5_dp*(S_L+S_R))) * &
                    (rx*Alpha(X_DIR)+ry*Alpha(Y_DIR)+rz*Alpha(Z_DIR))))/(S_R-S_L)

        !! alleviate carbuncles at sharp grid-aligned shock fronts
        smooth = SQRT(abs(p_L - p_R)/(p_L + p_R))
        F = F + smooth * ((S_R*F_L - S_L*F_R + S_L*S_R*(U_R - U_L)) / (S_R - S_L) - F)
    end if

end function

pure function yriemann(uL,uR) result(sf)

    use equations_const_mod

    real(dp), intent(in) :: uL(N_VARS), uR(N_VARS)
    real(dp)             :: sf(N_VARS)

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

    rf = xriemann(ruL,ruR)

    sf(DENS_VAR) = rf(DENS_VAR)
    sf(MOMX_VAR) = rf(MOMZ_VAR)
    sf(MOMY_VAR) = rf(MOMX_VAR)
    sf(MOMZ_VAR) = rf(MOMY_VAR)
    sf(ENER_VAR) = rf(ENER_VAR)
    
end function

pure function zriemann(uL,uR) result(sf)

    use equations_const_mod

    real(dp), intent(in) :: uL(N_VARS), uR(N_VARS)
    real(dp)             :: sf(N_VARS)

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

    rf = xriemann(ruL,ruR)

    sf(DENS_VAR) = rf(DENS_VAR)
    sf(MOMX_VAR) = rf(MOMY_VAR)
    sf(MOMY_VAR) = rf(MOMZ_VAR)
    sf(MOMZ_VAR) = rf(MOMX_VAR)
    sf(ENER_VAR) = rf(ENER_VAR)
    
end function

end module
