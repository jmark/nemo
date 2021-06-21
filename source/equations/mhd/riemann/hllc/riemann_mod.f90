# Error The HLLC Riemann solver does not work with GLM divergence cleaning. Do not use this solver.

module riemann_mod

use constants_mod

contains

pure SUBROUTINE RiemannSolverByHLLC(ConsL,ConsR,Flux)

    use equations_mod, only: xflux
    USE equations_mod, ONLY: Cons2Prim

    use equations_mod, only: FastestWave1D
    use equations_mod, only: FastestWave1D_Roe

    USE Equations_share_mod, ONLY: GLM_ch
    USE setup_config_mod, ONLY: Kappa
    USE setup_config_mod, ONLY: mu0

    ! real(dp), parameter :: KappaM1 = kappa-1.0_dp
    ! real(dp), parameter :: KappaM2 = kappa-2.0_dp
    ! real(dp), parameter :: sKappaM1 = 1.0_dp/(kappa-1.0_dp)

    real(dp), parameter :: smu_0 = 1.0_dp/mu0
    real(dp), parameter :: s2mu_0 = 0.5_dp/mu0

    !----------------------------------------------------------------------------------------------------------------------------------
    ! INPUT VARIABLES
    !----------------------------------------------------------------------------------------------------------------------------------
    ! INPUT VARIABLES
    REAL(dp),INTENT(IN)  :: ConsL(1:N_VARS) !<  left conservative state  
    REAL(dp),INTENT(IN)  :: ConsR(1:N_VARS) !< right conservative state
    !----------------------------------------------------------------------------------------------------------------------------------
    ! OUTPUT VARIABLES
    REAL(dp),INTENT(OUT) :: Flux(1:N_VARS) !<numerical flux
    !----------------------------------------------------------------------------------------------------------------------------------
    ! LOCAL VARIABLES
    !----------------------------------------------------------------------------------------------------------------------------------
    REAL(dp)             :: FluxL(N_VARS), FluxR(N_VARS)
    REAL(dp)             :: PrimL(1:N_VARS), PrimR(1:N_VARS) !< left/right primitive state
    REAL(dp)             :: U_HLL(N_VARS)
    REAL(dp)             :: U_star_L(N_VARS), U_star_R(N_VARS)
    REAL(dp)             :: SL, SR, SM, s_SL_SM, s_SR_SM
    REAL(dp)             :: ptot_L, p_star
    REAL(dp)             :: rho_L, rho_R
    REAL(dp)             :: srho_HLL
    REAL(dp)             :: Bx_star,By_star,Bz_star,E_star
    REAL(dp)             :: vx_L, vx_R
    REAL(dp)             :: Bx_L, Bx_R
    REAL(dp)             :: cf_L,cf_R,cf_Roe,RoeVelx 

    !==================================================================================================================================

    ! CALL ConsToPrim(PrimL(:),ConsL(:))
    ! CALL ConsToPrim(PrimR(:),ConsR(:))

    PrimL(:) = cons2prim(ConsL(:))
    PrimR(:) = cons2prim(ConsR(:))

    !--------------------------------------------------------------------------------
    !use Roe mean wavespeeds from  Roe meanvalues 
    !   (paper by Cargo & Gallice: "Roe Matrices for Ideal MHD and ...",1997)

    CALL FastestWave1D(PrimL,cf_L)
    CALL FastestWave1D_Roe(ConsL,ConsR,PrimL,PrimR,RoeVelx,cf_Roe)

    SL = MIN(PrimL(2)-cf_L, RoeVelx-cf_Roe   )


    IF (SL .GT. 0.0)  THEN
      ! CALL EvalAdvectionFlux1D(ConsL,Flux)
        flux = xflux(consL)
      RETURN
    END IF

    CALL FastestWave1D(PrimR,cf_R)
    SR = MAX(RoeVelx+cf_Roe, PrimR(2)+cf_R   )

    IF (SR .LT. 0.0) THEN
      !CALL EvalAdvectionFlux1D(ConsR,Flux)
        flux = xflux(consR)
      RETURN
    END IF

    ! NOW SL<0 and SR>0

    !CALL EvalAdvectionFlux1D(ConsL,FluxL)
    !CALL EvalAdvectionFlux1D(ConsR,FluxR)

    fluxL = xflux(consL)
    fluxR = xflux(consR)

    ! CALL EvalHLLState(ConsL,ConsR,SL,SR,FluxL,FluxR,U_HLL)

    U_HLL = (SR*ConsR-SL*ConsL-FluxR+FluxL)/(SR-SL)

    rho_L = PrimL(1)
    vx_L  = PrimL(2)
    Bx_L  = PrimL(6) 

    rho_R = PrimR(1)
    vx_R  = PrimR(2)
    Bx_R  = PrimR(6)

    ptot_L   = PrimL(5)+s2mu_0*SUM(PrimL(6:8)**2) !Total presssure!

    !SM = ((rho_L*vx_L*(vx_L-SL)+ptot_L-smu_0*Bx_L*Bx_L)-(rho_R*vx_R*(vx_R-SR)+ptot_R-smu_0*Bx_R*Bx_R)) &
    !     /(rho_L*(vx_L-SL)-rho_R*(vx_R-SR))
    sRho_HLL = 1./U_HLL(1)
    SM       = U_HLL(2)*sRho_HLL

    Bx_star  = U_HLL(6)
    By_star  = U_HLL(7)
    Bz_star  = U_HLL(8)
    p_star   = ptot_L + rho_L*(vx_L-SL)*(vx_L-SM)-s2mu_0*(Bx_L*Bx_L-Bx_star*Bx_star)
    E_star   = p_star*SM + smu_0*(-Bx_star*SUM(U_HLL(6:8)*U_HLL(2:4))*sRho_HLL  &
                                  +U_HLL(9)*(GLM_ch*Bx_star-0.5*U_HLL(9)*SM  )&
                                 )


    !-------------------------------------------------


    IF ((SL .LE. 0.0) .AND. (0.0 .LT. SM)) THEN
      s_SL_SM=1./(SL-SM)
      U_star_L(1) = rho_L*(SL-vx_L)*s_SL_SM
      U_star_L(2) = U_star_L(1)*SM
      U_star_L(3) = (ConsL(3)*SL-FluxL(3) -smu_0*Bx_star*By_star)*s_SL_SM
      U_star_L(4) = (ConsL(4)*SL-FluxL(4) -smu_0*Bx_star*Bz_star)*s_SL_SM
      U_star_L(5) = (E_star+ ConsL(5)*SL-FluxL(5) )*s_SL_SM
      U_star_L(6) = Bx_star 
      U_star_L(7) = By_star
      U_star_L(8) = Bz_star
      U_star_L(9) = ConsL(9)
      Flux = FluxL+SL*(U_star_L-ConsL)
    ELSE
      s_SR_SM=1./(SR-SM)
      U_star_R(1)  = rho_R*(SR-vx_R)*s_SR_SM
      U_star_R(2) = U_star_R(1)*SM
      U_star_R(3) = (ConsR(3)*SR-FluxR(3) -smu_0*Bx_star*By_star)*s_SR_SM
      U_star_R(4) = (ConsR(4)*SR-FluxR(4) -smu_0*Bx_star*Bz_star)*s_SR_SM
      U_star_R(5) = (E_star+ ConsR(5)*SR-FluxR(5) )*s_SR_SM
      U_star_R(6) = Bx_star 
      U_star_R(7) = By_star
      U_star_R(8) = Bz_star
      U_star_R(9) = ConsR(9)!0.0
      Flux = FluxR+SR*(U_star_R-ConsR)
    END IF

END SUBROUTINE RiemannSolverByHLLC

pure function xriemann(uL, uR) result(sf)

    real(dp), intent(in) :: uL(N_VARS), uR(N_VARS)
    real(dp)             :: sf(N_VARS)
    
    call RiemannSolverByHLLC(uL,uR,sf)

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
    
    ruL(MAGX_VAR) = uL(MAGY_VAR)
    ruL(MAGY_VAR) = uL(MAGZ_VAR)
    ruL(MAGZ_VAR) = uL(MAGX_VAR)

    ruL(GLMP_VAR) = uL(GLMP_VAR)


    ruR(DENS_VAR) = uR(DENS_VAR)

    ruR(MOMX_VAR) = uR(MOMY_VAR)
    ruR(MOMY_VAR) = uR(MOMZ_VAR)
    ruR(MOMZ_VAR) = uR(MOMX_VAR)

    ruR(ENER_VAR) = uR(ENER_VAR)

    ruR(MAGX_VAR) = uR(MAGY_VAR)
    ruR(MAGY_VAR) = uR(MAGZ_VAR)
    ruR(MAGZ_VAR) = uR(MAGX_VAR)

    ruR(GLMP_VAR) = uR(GLMP_VAR)

    rf = xriemann(ruL,ruR)

    sf(DENS_VAR) = rf(DENS_VAR)

    sf(MOMX_VAR) = rf(MOMZ_VAR)
    sf(MOMY_VAR) = rf(MOMX_VAR)
    sf(MOMZ_VAR) = rf(MOMY_VAR)

    sf(ENER_VAR) = rf(ENER_VAR)
    
    sf(MAGX_VAR) = rf(MAGZ_VAR)
    sf(MAGY_VAR) = rf(MAGX_VAR)
    sf(MAGZ_VAR) = rf(MAGY_VAR)

    sf(GLMP_VAR) = rf(GLMP_VAR)

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
    
    ruL(MAGX_VAR) = uL(MAGZ_VAR)
    ruL(MAGY_VAR) = uL(MAGX_VAR)
    ruL(MAGZ_VAR) = uL(MAGY_VAR)

    ruL(GLMP_VAR) = uL(GLMP_VAR)

 
    ruR(DENS_VAR) = uR(DENS_VAR)

    ruR(MOMX_VAR) = uR(MOMZ_VAR)
    ruR(MOMY_VAR) = uR(MOMX_VAR)
    ruR(MOMZ_VAR) = uR(MOMY_VAR)

    ruR(ENER_VAR) = uR(ENER_VAR)

    ruR(MAGX_VAR) = uR(MAGZ_VAR)
    ruR(MAGY_VAR) = uR(MAGX_VAR)
    ruR(MAGZ_VAR) = uR(MAGY_VAR)

    ruR(GLMP_VAR) = uR(GLMP_VAR)

    rf = xriemann(ruL,ruR)

    sf(DENS_VAR) = rf(DENS_VAR)

    sf(MOMX_VAR) = rf(MOMY_VAR)
    sf(MOMY_VAR) = rf(MOMZ_VAR)
    sf(MOMZ_VAR) = rf(MOMX_VAR)

    sf(ENER_VAR) = rf(ENER_VAR)
    
    sf(MAGX_VAR) = rf(MAGY_VAR)
    sf(MAGY_VAR) = rf(MAGZ_VAR)
    sf(MAGZ_VAR) = rf(MAGX_VAR)

    sf(GLMP_VAR) = rf(GLMP_VAR)

end function

end module
