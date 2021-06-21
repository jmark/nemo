module riemann_mod

use constants_mod

contains

pure function xriemann(U_L,U_R) result(F)

    use equations_const_mod
    use equations_mod, only: pressure
    use setup_config_mod, only: kappa

    use equations_mod, only: flux => xflux

    real(dp), dimension(N_VARS), intent(in) :: U_L, U_R
    real(dp), dimension(N_VARS)             :: F
    
    real(dp) :: F_L(N_VARS), F_R(N_VARS)

    real(dp) :: c_L,c_R !! sound speed
    real(dp) :: H_L,H_R !! enthalpy
    real(dp) :: p_L,p_R !! pressure
    real(dp) :: S_L,S_R !! signal wave speed
    real(dp) :: V_L,V_R !! velocities

    real(dp) :: HTilde, vTilde, aTilde
    real(dp) :: SStar, QStar, UStar(N_VARS)

    real(dp) :: smooth

    real(dp) :: F_HLLE(N_VARS)

    real(dp), parameter :: beta = SQRT(0.5_dp*(kappa-1)/kappa)

    !! velocity
    V_L = U_L(MOMX_VAR)/U_L(DENS_VAR)
    V_R = U_R(MOMX_VAR)/U_R(DENS_VAR)

    !! pressure
    p_L = pressure(U_L)
    p_R = pressure(U_R)

    !! sound speed
    c_L = SQRT(kappa*p_L/U_L(DENS_VAR))
    c_R = SQRT(kappa*p_R/U_R(DENS_VAR))

    !! Signal velocity estimates by Davis (1988) and Einfeldt (1988)

    !! enthalpy
    H_L = (U_L(ENER_VAR) + p_L)/U_L(DENS_VAR)
    H_R = (U_R(ENER_VAR) + p_R)/U_R(DENS_VAR)

    !! Roe averages
    HTilde = (SQRT(U_L(DENS_VAR))*H_L + SQRT(U_R(DENS_VAR))*H_R)/(SQRT(U_L(DENS_VAR)) + SQRT(U_R(DENS_VAR)))
    vTilde = (SQRT(U_L(DENS_VAR))*v_L + SQRT(U_R(DENS_VAR))*v_R)/(SQRT(U_L(DENS_VAR)) + SQRT(U_R(DENS_VAR)))
    aTilde = SQRT((kappa-1.0_dp)*(HTilde - 0.5_dp*vTilde**2))

    !! signal velocities
    S_L = MIN(vTilde - aTilde, V_L - beta*c_L, 0.0_dp)
    S_R = MAX(vTilde + aTilde, V_R + beta*c_R, 0.0_dp)

    if (S_L .GE. 0.0_dp) then
        F = flux(U_L)

    elseif (S_R .LE. 0.0_dp) then
        F = flux(U_R)

    else
        F_L = flux(U_L)
        F_R = flux(U_R)

        smooth = SQRT(abs(p_L - p_R)/(p_L + p_R))
        F_HLLE = (S_R*F_L - S_L*F_R + S_L*S_R*(U_R - U_L)) / (S_R - S_L)

        SStar = (p_R - p_L + U_L(MOMX_VAR)*(S_L - V_L) - U_R(MOMX_VAR)*(S_R - V_R)) / (U_L(DENS_VAR)*(S_L - V_L) - U_R(DENS_VAR)*(S_R - V_R))

        if (S_L .le. 0.0_dp .and. 0.0_dp .le. SStar) then
            QStar = U_L(DENS_VAR)*(S_L - V_L)/(S_L-SStar)

            UStar(DENS_VAR) = QStar
            UStar(MOMX_VAR) = QStar * SStar
            UStar(MOMY_VAR) = QStar * U_L(MOMY_VAR)/U_L(DENS_VAR)
            UStar(MOMZ_VAR) = QStar * U_L(MOMZ_VAR)/U_L(DENS_VAR)
            UStar(ENER_VAR) = QStar * (U_L(ENER_VAR)/U_L(DENS_VAR) + (SStar-V_L)*(SStar + p_L/U_L(DENS_VAR)/(S_L - V_L)))

            F = smooth*F_HLLE + (1.0_dp-smooth)*(F_L + S_L*(UStar - U_L))
            !F = F_L + S_L*(UStar - U_L)

        else
            QStar = U_R(DENS_VAR)*(S_R - V_R)/(S_R-SStar)

            UStar(DENS_VAR) = QStar
            UStar(MOMX_VAR) = QStar * SStar
            UStar(MOMY_VAR) = QStar * U_R(MOMY_VAR)/U_R(DENS_VAR)
            UStar(MOMZ_VAR) = QStar * U_R(MOMZ_VAR)/U_R(DENS_VAR)
            UStar(ENER_VAR) = QStar * (U_R(ENER_VAR)/U_R(DENS_VAR) + (SStar-V_R)*(SStar + p_R/U_R(DENS_VAR)/(S_R - V_R)))

            F = smooth*F_HLLE + (1.0_dp-smooth)*(F_R + S_R*(UStar - U_R))
            !F = F_R + S_R*(UStar - U_R)

        end if 
    end if 

end function

pure function yriemann(inflow, exflow) result(sflux)

    use equations_const_mod

    real(dp), dimension(N_VARS), intent(in) :: inflow, exflow
    real(dp), dimension(N_VARS)             :: sflux

    real(dp) :: rinflow(N_VARS), rexflow(N_VARS), rflux(N_VARS) !! rotated

    rinflow(DENS_VAR) = inflow(DENS_VAR)
    rinflow(MOMX_VAR) = inflow(MOMY_VAR)
    rinflow(MOMY_VAR) = inflow(MOMZ_VAR)
    rinflow(MOMZ_VAR) = inflow(MOMX_VAR)
    rinflow(ENER_VAR) = inflow(ENER_VAR)
    
    rexflow(DENS_VAR) = exflow(DENS_VAR)
    rexflow(MOMX_VAR) = exflow(MOMY_VAR)
    rexflow(MOMY_VAR) = exflow(MOMZ_VAR)
    rexflow(MOMZ_VAR) = exflow(MOMX_VAR)
    rexflow(ENER_VAR) = exflow(ENER_VAR)

    rflux = xriemann(rinflow,rexflow)

    sflux(DENS_VAR) = rflux(DENS_VAR)
    sflux(MOMX_VAR) = rflux(MOMZ_VAR)
    sflux(MOMY_VAR) = rflux(MOMX_VAR)
    sflux(MOMZ_VAR) = rflux(MOMY_VAR)
    sflux(ENER_VAR) = rflux(ENER_VAR)
    
end function

pure function zriemann(inflow, exflow) result(sflux)

    use equations_const_mod

    real(dp), dimension(N_VARS), intent(in) :: inflow, exflow
    real(dp), dimension(N_VARS)             :: sflux

    real(dp) :: rinflow(N_VARS), rexflow(N_VARS), rflux(N_VARS) !! rotated

    rinflow(DENS_VAR) = inflow(DENS_VAR)
    rinflow(MOMX_VAR) = inflow(MOMZ_VAR)
    rinflow(MOMY_VAR) = inflow(MOMX_VAR)
    rinflow(MOMZ_VAR) = inflow(MOMY_VAR)
    rinflow(ENER_VAR) = inflow(ENER_VAR)
    
    rexflow(DENS_VAR) = exflow(DENS_VAR)
    rexflow(MOMX_VAR) = exflow(MOMZ_VAR)
    rexflow(MOMY_VAR) = exflow(MOMX_VAR)
    rexflow(MOMZ_VAR) = exflow(MOMY_VAR)
    rexflow(ENER_VAR) = exflow(ENER_VAR)

    rflux = xriemann(rinflow,rexflow)

    sflux(DENS_VAR) = rflux(DENS_VAR)
    sflux(MOMX_VAR) = rflux(MOMY_VAR)
    sflux(MOMY_VAR) = rflux(MOMZ_VAR)
    sflux(MOMZ_VAR) = rflux(MOMX_VAR)
    sflux(ENER_VAR) = rflux(ENER_VAR)
    
end function

end module
