module riemann_mod

use constants_mod

contains

!! Euler flux
pure function xflux(u,pres) result(f)

    use equations_const_mod

    real(dp), intent(in)    :: u(N_VARS), pres
    real(dp)                :: f(N_VARS)

    real(dp) :: velx

    velx = u(MOMX_VAR)/u(DENS_VAR)

    f(DENS_VAR) = u(MOMX_VAR)
    f(MOMX_VAR) = velx * u(MOMX_VAR) + pres
    f(MOMY_VAR) = velx * u(MOMY_VAR)
    f(MOMZ_VAR) = velx * u(MOMZ_VAR)
    f(ENER_VAR) = velx * (u(ENER_VAR) + pres)

end function

pure function xriemann(U_L,U_R) result(F)

    use equations_const_mod
    use setup_config_mod, only: kappa
    use equations_mod, only: pressure

    real(dp), intent(in) :: U_L(N_VARS), U_R(N_VARS)
    real(dp)             :: F(N_VARS)
    
    real(dp) :: F_L(N_VARS), F_R(N_VARS)
    real(dp) :: F_HLLE(N_VARS), F_ROEF(N_VARS)

    real(dp) :: p_L,p_R !! pressure
    real(dp) :: c_L,c_R !! pressure
    real(dp) :: H_L,H_R !! enthalpy
    real(dp) :: S_L,S_R !! signal wave speed
    real(dp) :: V_L(N_DIMS),V_R(N_DIMS) !! velocities

    real(dp) :: sqrtrho_L, sqrtrho_R, ssqrtrho
    real(dp) :: roevel(N_DIMS),roeh,roec,absvel,RoeDens

    real(dp) :: dU(N_VARS)      !! jumps
    real(dp) :: alpha(N_VARS)   !! wave strenghts

    real(dp) :: lambda(N_VARS), a_L(N_VARS), a_R(N_VARS), da !! roe eigenvalues
    real(dp) :: r1(N_VARS),r2(N_VARS),r3(N_VARS),r4(N_VARS),r5(N_VARS)  !! roe eigenvectors

    real(dp) :: smooth

    integer :: h

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

    !! central fluxes
    F_L = xflux(U_L,p_L)
    F_R = xflux(U_R,p_R)

    !! ---------------------------------------------------------------------- !!

    !! auxiliaries for Roe averages
    SqrtRho_L = SQRT(U_L(DENS_VAR))
    SqrtRho_R = SQRT(U_R(DENS_VAR))
    sSqrtRho  = 1.0_dp/(SqrtRho_L+SqrtRho_R)

    !! Roe averages
    RoeVel    = (SqrtRho_R*V_R + SqrtRho_L*V_L) * sSqrtRho
    absVel    = DOT_PRODUCT(RoeVel,RoeVel)

    !! Roe averages
    RoeH      = (SqrtRho_R*H_R+SqrtRho_L*H_L) * sSqrtRho
    RoeC      = SQRT((kappa-1)*(RoeH-0.5_dp*absVel))

    RoeDens   = SQRT(U_L(DENS_VAR)*U_R(DENS_VAR))

    !! ---------------------------------------------------------------------- !!
    !! Signal velocity estimates by Davis (1988) and Einfeldt (1988)

    S_L = RoeVel(X_DIR) - RoeC
    S_R = RoeVel(X_DIR) + RoeC

    if (S_L .GE. 0.0_dp) then
        F_HLLE = F_L
    elseif (S_R .LE. 0.0_dp) then
        F_HLLE = F_R
    else
        F_HLLE = (S_R*F_L - S_L*F_R + S_L*S_R*(U_R - U_L)) / (S_R - S_L)
    end if 

    !! ---------------------------------------------------------------------- !!
    !! eigenvalues

    lambda  = (/ RoeVel(X_DIR)-RoeC,RoeVel(X_DIR),RoeVel(X_DIR),RoeVel(X_DIR),RoeVel(X_DIR)+RoeC/)

    !! eigenvectors
    r1 = (/ 1.0_dp,     lambda(1),  RoeVel(2),  RoeVel(3),  RoeH-RoeVel(1)*Roec /)
    r2 = (/ 1.0_dp,     RoeVel(1),  RoeVel(2),  RoeVel(3),  0.5*absVel          /)
    r3 = (/ 0.0_dp,     0.0_dp,     1.0_dp,     0.0_dp,     RoeVel(2)           /)
    r4 = (/ 0.0_dp,     0.0_dp,     0.0_dp,     1.0_dp,     RoeVel(3)           /)
    r5 = (/ 1.0_dp,     lambda(4),  RoeVel(2),  RoeVel(3),  RoeH+RoeVel(1)*Roec /)

    !! jumps
    dU(1)   = U_R(DENS_VAR) - U_L(DENS_VAR)
    dU(2:4) = V_R - V_L
    dU(5)   = p_R - p_L

    !! wave strenghts
    Alpha(1) = 0.5_dp/(Roec*Roec)*(dU(5) - RoeDens*Roec*dU(2))
    Alpha(2) = dU(1) - dU(5)/(Roec*Roec)
    Alpha(3) = RoeDens*dU(3)
    Alpha(4) = RoeDens*dU(4)
    Alpha(5) = 0.5_dp/(Roec*Roec)*(dU(5) + RoeDens*Roec*dU(2))

    !! entropy fix
    a_L(1) = V_L(X_DIR) - c_L
    a_L(2) = V_L(X_DIR)
    a_L(3) = V_L(X_DIR)
    a_L(4) = V_L(X_DIR)
    a_L(5) = V_L(X_DIR) + c_L

    a_R(1) = V_R(X_DIR) - c_R
    a_R(2) = V_R(X_DIR)
    a_R(3) = V_R(X_DIR)
    a_R(4) = V_R(X_DIR)
    a_R(5) = V_R(X_DIR) + c_R

    do h = 1,N_VARS
        da = max(0.0_dp, lambda(h) - a_L(h), a_R(h) - lambda(h))
        if (abs(lambda(h)) .lt. da) then
            lambda(h) = 0.5_dp*(lambda(h)*lambda(h)/da + da)
        else
            lambda(h) = abs(lambda(h))
        end if
    end do 

    !! assemble flux
    F_ROEF = 0.5_dp*((F_L + F_R) - &
           Alpha(1)*lambda(1)*r1 - &
           Alpha(2)*lambda(2)*r2 - &
           Alpha(3)*lambda(3)*r3 - &
           Alpha(4)*lambda(4)*r4 - &
           Alpha(5)*lambda(5)*r5)

    !! ---------------------------------------------------------------------- !!

    smooth = SQRT(abs(p_L - p_R)/(p_L + p_R))
    !smooth = SQRT(abs(p_L/U_L(DENS_VAR) - p_R/U_R(DENS_VAR))/(p_L/U_L(DENS_VAR) + p_R/U_R(DENS_VAR)))

    F = smooth*F_HLLE + (1.0_dp-smooth)*F_ROEF

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
