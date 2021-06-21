# define PP_GLM 1

module riemann_mod

use constants_mod
use setup_config_mod, only: kappa
use setup_config_mod, only: mu0
use equations_share_mod, only: glm_ch

integer, parameter :: PP_nVar = N_VARS

real(dp), parameter :: KappaM1 = kappa-1.0_dp
real(dp), parameter :: KappaM2 = kappa-2.0_dp
real(dp), parameter :: sKappaM1 = 1.0_dp/(kappa-1.0_dp)

real(dp), parameter :: smu_0 = 1.0_dp/mu0
real(dp), parameter :: s2mu_0 = 0.5_dp/mu0

contains

pure function xriemann(U_L,U_R) result(F)

    real(dp), intent(in) :: U_L(N_VARS), U_R(N_VARS)
    real(dp)             :: F(N_VARS)
    
    call EntropyStable9WaveFlux(U_L,U_R,F)

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

!! Following code has been copied from the 'FLUXO' code project.

!==================================================================================================================================
!> The following three routines are for the entropy stable 9wave solver. 
!>    See: Derigs et al. (2018). "Ideal GLM-MHD: About the entropy consistent nine-wave magnetic field divergence diminishing ideal
!>         magnetohydrodynamics equations. Journal of Computational Physics, 364, 420-467."
!> ATTENTION: 1) Central flux is entropy conservaing for MHD, kinetic Energy conservation only in the Euler case
!>            2) mu_0 added, total energy contribution is 1/(2mu_0)(|B|^2+psi^2), in energy flux: 1/mu_0*(B.B_t + psi*psi_t) 
!>            3) Dissipation for each characteristic wave seperately
!==================================================================================================================================

pure SUBROUTINE EntropyStable9WaveFlux(UL,UR,Fstar)

! ! MODULES
! USE MOD_PreProc
! USE MOD_Flux_Average, ONLY:LN_MEAN
! USE MOD_Equation_Vars,ONLY:kappa,kappaM1,skappaM1,smu_0,s2mu_0,consToEntropy

IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL(dp),DIMENSION(PP_nVar),INTENT(IN)  :: UL      !< left state
REAL(dp),DIMENSION(PP_nVar),INTENT(IN)  :: UR      !< right state
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL(dp),DIMENSION(PP_nVar),INTENT(OUT) :: Fstar   !<  flux in x
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL(dp)            :: beta_R,beta_L
REAL(dp)            :: srho_L,srho_R
REAL(dp)            :: u2_L,u2_R
REAL(dp)            :: p_L,p_R
REAL(dp)            :: u_L(3),u_R(3)
REAL(dp)            :: V_jump(PP_nVar)
REAL(dp)            :: Dmatrix(PP_nVar),Rmatrix(PP_nVar,PP_nVar),RT(PP_nVar,PP_nVar)

!==================================================================================================================================
ASSOCIATE(  rho_L =>   UL(1),  rho_R =>   UR(1), &
!
           rhoU_L => UL(2:4), rhoU_R => UR(2:4), &
#ifdef PP_GLM
              E_L =>UL(5)-s2mu_0*UL(9)**2, E_R =>UR(5)-s2mu_0*UR(9)**2, &
#else
              E_L =>UL(5)   ,    E_R =>UR(5), &
#endif
              B_L => UL(6:8),    B_R => UR(6:8)  )

! Get auxiliar variables
srho_L=1./rho_L
srho_R=1./rho_R
u_L = rhoU_L(:)*srho_L
u_R = rhoU_R(:)*srho_R
u2_L = SUM(u_L(:)*u_L(:))
u2_R = SUM(u_R(:)*u_R(:))
p_L    = kappaM1*(E_L - 0.5*(rho_L*u2_L+smu_0*SUM(B_L(:)*B_L(:))))
p_R    = kappaM1*(E_R - 0.5*(rho_R*u2_R+smu_0*SUM(B_R(:)*B_R(:))))
beta_L = 0.5*rho_L/p_L !beta=rho/(2*p)
beta_R = 0.5*rho_R/P_R

! Jump in entropy vars
V_jump(1)         = (kappa*(LOG(rho_R*srho_L))-LOG(p_R/p_L))*skappaM1 &
                       - (beta_R*u2_R         -beta_L*u2_L  )  
V_jump(2:4)       =  2.0*(beta_R*u_R(:)       -beta_L*u_L(:))        ! 2*beta*v
V_jump(5)         = -2.0*(beta_R              -beta_L       )        !-2*beta
V_jump(6:PP_nVar) =  2.0*(beta_R*UR(6:PP_nVar)-beta_L*UL(6:PP_nVar)) ! 2*beta*B
    
call EntropyStable9WaveFlux_VolFluxAndDissipMatrices(UL,UR,Fstar,Dmatrix,Rmatrix)
RT = TRANSPOSE(Rmatrix)

! Compute entropy-stable fluxes
Fstar = Fstar - 0.5*MATMUL(Rmatrix,Dmatrix*MATMUL(RT,V_jump))

END ASSOCIATE 
END SUBROUTINE EntropyStable9WaveFlux

!==================================================================================================================================
!> Computation of the entropy conserving flux and the dissipation matrices for the entropy stable 9wave solver
!==================================================================================================================================
pure subroutine EntropyStable9WaveFlux_VolFluxAndDissipMatrices(UL,UR,Fstar,Dmatrix,Rmatrix)

!!   USE MOD_Equation_Vars,ONLY:kappaM1, smu_0, s2mu_0, sKappaM1, Kappa
!!   USE MOD_Flux_Average, ONLY:LN_MEAN
!! #ifdef PP_GLM
!!   USE MOD_Equation_Vars,ONLY:GLM_ch
!! #endif

  real(dp), intent(in)  :: UL(PP_nVar)      !< left state
  real(dp), intent(in)  :: UR(PP_nVar)      !< right state
  real(dp), intent(out) :: Fstar  (PP_nVar)
  real(dp), intent(out) :: Dmatrix(PP_nVar)
  real(dp), intent(out) :: Rmatrix(PP_nVar,PP_nVar)

  !-local-variables-----------------------------------------
  ! Mean states
  real(dp)    :: sbetaLN,betaAvg,rhoLN,pAvg,psiAvg,BAvg(3),uAvg(3)
  real(dp)    :: p_L,p_R,u_L(3), u_R(3), pTilde
  ! Additional

  REAL(dp) :: B2_L,B2_R,B2Avg
  REAL(dp) :: u2_L,u2_R,u2Avg
  REAL(dp) :: beta_R,beta_L
  REAL(dp) :: srho_L,srho_R
  REAL(dp) :: u1_B2Avg, uB_Avg
  real(dp) :: Tmatrix(PP_nVar)
  REAL(dp) :: aA,aLN,abeta,bb1A,bb2A,bb3A,bbA,aabbA,ca,xx,xxx,cf,cs,bperpA,beta2A,beta3A,alphaf,alphas,sgnb1
  REAL(dp) :: psiSplus,psiSminus,psiFplus,psiFminus, srhoLN, rhoAvg, pLN
  real(dp) :: LambdaMax, phi

  !---------------------------------------------------------
  
  ASSOCIATE(  rho_L =>   UL(1),  rho_R =>   UR(1), &
             rhoU_L => UL(2:4), rhoU_R => UR(2:4), &
#ifdef PP_GLM
                E_L =>UL(5)-0.5*smu_0*UL(9)**2, E_R =>UR(5)-0.5*smu_0*UR(9)**2, &
              psi_L =>UL(9)   ,  psi_R =>UR(9), &
#else
                E_L =>UL(5)   ,    E_R =>UR(5), &
#endif
                B_L => UL(6:8),    B_R => UR(6:8)  )
  ! Get auxiliar variables
  srho_L=1./rho_L
  srho_R=1./rho_R
  u_L = rhoU_L(:)*srho_L
  u_R = rhoU_R(:)*srho_R
  u2_L = SUM(u_L(:)*u_L(:))
  u2_R = SUM(u_R(:)*u_R(:))
  B2_L = SUM(B_L(:)*B_L(:))
  B2_R = SUM(B_R(:)*B_R(:))
  p_L    = kappaM1*(E_L - 0.5*(rho_L*u2_L+smu_0*B2_L))
  p_R    = kappaM1*(E_R - 0.5*(rho_R*u2_R+smu_0*B2_R))
  beta_L = 0.5*rho_L/p_L !beta=rho/(2*p)
  beta_R = 0.5*rho_R/P_R

  ! Get the averages for the numerical flux
  rhoLN      = LN_MEAN( rho_L, rho_R)
  sbetaLN    = 1./LN_MEAN(beta_L,beta_R)
  uAvg       = 0.5 * ( u_L +  u_R)
  u2Avg      = 0.5 * (u2_L + u2_R)
  BAvg       = 0.5 * ( B_L +  B_R)
  B2Avg      = 0.5 * (B2_L + B2_R)
  u1_B2Avg   = 0.5 * (u_L(1)*B2_L       + u_R(1)*B2_R)
  uB_Avg     = 0.5 * (SUM(u_L(:)*B_L(:))+ SUM(u_R(:)*B_R(:)))
  betaAvg    = 0.5 * (beta_L + beta_R)                                                               
  pAvg       = 0.5*(rho_L+rho_R)/(beta_L+beta_R) !rhoMEAN/(2*betaMEAN)
#ifdef PP_GLM
  psiAvg     = 0.5*(psi_L+psi_R)
#endif
  pTilde     = pAvg+ s2mu_0*B2Avg !+1/(2mu_0)({{|B|^2}}...)
  
  ! Entropy conserving and kinetic energy conserving flux
Fstar(1) = rhoLN*uAvg(1)
Fstar(2) = Fstar(1)*uAvg(1) - smu_0*Bavg(1)*BAvg(1) + pTilde
Fstar(3) = Fstar(1)*uAvg(2) - smu_0*Bavg(1)*BAvg(2)
Fstar(4) = Fstar(1)*uAvg(3) - smu_0*Bavg(1)*BAvg(3)
Fstar(7) = uAvg(1)*Bavg(2) - BAvg(1)*uAvg(2)
Fstar(8) = uAvg(1)*Bavg(3) - BAvg(1)*uAvg(3)
#ifdef PP_GLM
Fstar(6) = GLM_ch*psiAvg
Fstar(9) = GLM_ch*BAvg(1)
#endif

Fstar(5) = Fstar(1)*0.5*(skappaM1*sbetaLN - u2Avg)  &
           + SUM(uAvg(:)*Fstar(2:4)) &
           +smu_0*( SUM(BAvg(:)*Fstar(6:8)) &
                   -0.5*u1_B2Avg +BAvg(1)*uB_Avg &
#ifdef PP_GLM
                   +Fstar(9)*psiAvg-GLM_ch*0.5*(B_L(1)*psi_L+B_R(1)*psi_R)     &
#endif
                   )
  ! MHD waves selective dissipation

    ! Compute additional averages
    rhoAvg = 0.5*(rho_L + rho_R)
    pLN = 0.5*rhoLN*sbetaLN
    u2Avg = u_L(1)*u_R(1) + u_L(2)*u_R(2) + u_L(3)*u_R(3)  ! Redefinition of u2Avg -> \bar{\norm{u}^2}
    srhoLN = 1./rhoLN
    aA = SQRT(kappa*pAvg/rhoLN)
    aLN = SQRT(kappa*pLN/rhoLN)
    abeta = SQRT(0.5*kappa/betaAvg)
    bb1A = BAvg(1)/SQRT(rhoLN)
    bb2A = BAvg(2)/SQRT(rhoLN)
    bb3A = BAvg(3)/SQRT(rhoLN)
    bbA = SQRT(bb1A*bb1A + bb2A*bb2A + bb3A*bb3A)
    aabbA = aA*aA + bbA*bbA

    ! Alfven speed
    ! Derigs et al. (2018), (4.63)
    ca = ABS(bb1A)

    ! Control round-off errors
    xx = aabbA*aabbA - 4.0*aA*aA*bb1A*bb1A
    if(xx .lt. 0.0) xx = 0.0
    xxx = aabbA + SQRT(xx)

    ! Fast magnetoacoustic speed
    ! Derigs et al. (2018), (4.63)
    cf = SQRT(0.5 * xxx)

    ! Control round-off errors
    xxx = aabbA - SQRT(xx)
    if(xxx .lt. 0.0) xxx = 0.0

    ! Slow magnetoacoustic speed
    ! Derigs et al. (2018), (4.63)
    cs = SQRT(0.5 * xxx)

    bperpA = SQRT(bb2A*bb2A + bb3A*bb3A)
    ! In case of very small bperpA, the values of betaA_{2,3}
    ! are indeterminable so we make them pairwise orthogonal
    ! Derigs et al. (2018), (4.63)
    if(bperpA .gt. 1.0E-14) then
      beta2A = bb2A/bperpA
      beta3A = bb3A/bperpA
    else
      bperpA = 0.0
      beta2A = 1.0/sqrt(2.0)
      beta3A = 1.0/sqrt(2.0)
    endif

    ! Avoid negative round-off errors when computing alphaf (auxiliary variable)
    xx = 0.0
    if((cf*cf-cs*cs) .gt. 0.0) xx = (aA*aA-cs*cs)/(cf*cf-cs*cs)
    if(xx .gt. 0.0) then
      alphaf = SQRT(xx)
    else
      alphaf = 0.0
    end if

    ! Avoid negative round-off errors when computing alphas (auxiliary variable)
    xx = 0.0
    if((cf*cf-cs*cs) .gt. 0.0) xx = (cf*cf-aA*aA)/(cf*cf-cs*cs)
    if(xx .gt. 0.0) then
      alphas = SQRT(xx)
    else
      alphas = 0.0
    end if

    sgnb1 = sign(1.,bb1A)

    ! Derigs et al. (2018), (4.63)
    psiSplus =  0.5_dp*alphas*rhoLN*u2Avg - abeta*alphaf*rhoLN*bperpA + &
                alphas*rhoLN*aLN*aLN*skappaM1 + alphas*cs*rhoLN*uAvg(1) + &
                alphaf*cf*rhoLN*sgnb1*(uAvg(2)*beta2A + uAvg(3)*beta3A)
    psiSminus = 0.5_dp*alphas*rhoLN*u2Avg - abeta*alphaf*rhoLN*bperpA + &
                alphas*rhoLN*aLN*aLN*skappaM1 - alphas*cs*rhoLN*uAvg(1) - &
                alphaf*cf*rhoLN*sgnb1*(uAvg(2)*beta2A + uAvg(3)*beta3A)

    psiFplus =  0.5_dp*alphaf*rhoLN*u2avg + abeta*alphas*rhoLN*bperpA + &
                alphaf*rhoLN*aLN*aLN*skappaM1 + alphaf*cf*rhoLN*uAvg(1) - &
                alphas*cs*rhoLN*sgnb1*(uAvg(2)*beta2A + uAvg(3)*beta3A)
    psiFminus = 0.5_dp*alphaf*rhoLN*u2avg + abeta*alphas*rhoLN*bperpA + &
                alphaf*rhoLN*aLN*aLN*skappaM1 - alphaf*cf*rhoLN*uAvg(1) + &
                alphas*cs*rhoLN*sgnb1*(uAvg(2)*beta2A + uAvg(3)*beta3A)

    ! + fast magnetoacoustic wave
    ! Derigs et al. (2018), (4.68)
    Rmatrix(:,1) = (/ alphaf*rhoLN, &
                      alphaf*rhoLN*(uAvg(1) + cf), &
                      rhoLN*(alphaf*uAvg(2) - alphas*cs*beta2A*sgnb1), &
                      rhoLN*(alphaf*uAvg(3) - alphas*cs*beta3A*sgnb1), &
                      psiFplus, &
                      0.0, &
                      alphas*abeta*beta2A*SQRT(rhoLN), &
                      alphas*abeta*beta3A*SQRT(rhoLN), &
                      0.0 /)

    ! + Alfven wave
    ! Derigs et al. (2018), (4.67)
    Rmatrix(:,2) = (/ 0.0, &
                      0.0,  &
                      rhoLN*SQRT(rhoAvg)*beta3A, &
                      -rhoLN*SQRT(rhoAvg)*beta2A, &
                      -rhoLN*SQRT(rhoAvg)*(beta2A*uAvg(3) - beta3A*uAvg(2)), &
                      0.0, &
                      -rhoLN*beta3A, &
                      rhoLN*beta2A, &
                      0.0 /)

    ! + slow magnetoacoustic wave
    ! Derigs et al. (2018), (4.69)
    Rmatrix(:,3) = (/ alphas*rhoLN, &
                      alphas*rhoLN*(uAvg(1) + cs), &
                      rhoLN*(alphas*uAvg(2) + alphaf*cf*beta2A*sgnb1), &
                      rhoLN*(alphas*uAvg(3) + alphaf*cf*beta3A*sgnb1), &
                      psiSplus, &
                      0.0, &
                      -alphaf*abeta*beta2A*SQRT(rhoLN), &
                      -alphaf*abeta*beta3A*SQRT(rhoLN), &
                      0.0 /)

    ! + GLM wave
    ! Dergs et al. (2018), eq. (4.65)
    Rmatrix(:,4) = (/ 0.0, &
                      0.0, &
                      0.0, &
                      0.0, &
                      BAvg(1) + psiAvg, &
                      1.0, &
                      0.0, &
                      0.0, &
                      1.0 /)

    ! Entropy wave
    ! Derigs et al. (2018), (4.66)
    Rmatrix(:,5) = (/ 1.0, &
                      uAvg(1), &
                      uAvg(2), &
                      uAvg(3), &
                      0.5*u2Avg, &
                      0.0, &
                      0.0, &
                      0.0, &
                      0.0 /)

    ! - GLM wave
    ! Dergs et al. (2018), eq. (4.65)
    Rmatrix(:,6) = (/ 0.0, &
                      0.0, &
                      0.0, &
                      0.0, &
                      BAvg(1) - psiAvg, &
                      1.0, &
                      0.0, &
                      0.0, &
                      -1.0 /)

    ! - slow magnetoacoustic wave
    ! Derigs et al. (2018), (4.69)
    Rmatrix(:,7) = (/ alphas*rhoLN, &
                      alphas*rhoLN*(uAvg(1) - cs), &
                      rhoLN*(alphas*uAvg(2) - alphaf*cf*beta2A*sgnb1), &
                      rhoLN*(alphas*uAvg(3) - alphaf*cf*beta3A*sgnb1), &
                      psiSminus, &
                      0.0, &
                      -alphaf*abeta*beta2A*SQRT(rhoLN), &
                      -alphaf*abeta*beta3A*SQRT(rhoLN), &
                      0.0 /)

    ! - Alfven wave
    ! Derigs et al. (2018), (4.67)
    Rmatrix(:,8) = (/ 0.0, &
                      0.0, &
                      -rhoLN*SQRT(rhoAvg)*beta3A, &
                      rhoLN*SQRT(rhoAvg)*beta2A, &
                      rhoLN*SQRT(rhoAvg)*(beta2A*uAvg(3) - beta3A*uAvg(2)), &
                      0.0, &
                      -rhoLN*beta3A, &
                      rhoLN*beta2A, &
                      0.0 /)

    ! - fast magnetoacoustic wave
    ! Derigs et al. (2018), (4.68)
    Rmatrix(:,9) = (/ alphaf*rhoLN, &
                      alphaf*rhoLN*(uAvg(1) - cf), &
                      rhoLN*(alphaf*uAvg(2) + alphas*cs*beta2A*sgnb1), &
                      rhoLN*(alphaf*uAvg(3) + alphas*cs*beta3A*sgnb1), &
                      psiFminus, &
                      0.0, &
                      alphas*abeta*beta2A*SQRT(rhoLN), &
                      alphas*abeta*beta3A*SQRT(rhoLN), &
                      0.0 /)
    
    ! A blend of the 9waves solver and LLF
      phi = sqrt(ABS(1.0-(p_R*rho_R/(p_L*rho_L)))/(1.0+(p_R*rho_R/(p_L*rho_L))))
      LambdaMax = MAX(ABS( uAvg(1) + cf ),ABS( uAvg(1) - cf ))
      
      Dmatrix(1) = (1.0_dp-phi)*ABS( uAvg(1) + cf ) + phi*LambdaMax
      Dmatrix(2) = (1.0_dp-phi)*ABS( uAvg(1) + ca ) + phi*LambdaMax
      Dmatrix(3) = (1.0_dp-phi)*ABS( uAvg(1) + cs ) + phi*LambdaMax
      Dmatrix(4) = (1.0_dp-phi)*ABS( uAvg(1) + GLM_ch ) + phi*LambdaMax
      Dmatrix(5) = (1.0_dp-phi)*ABS( uAvg(1)      ) + phi*LambdaMax
      Dmatrix(6) = (1.0_dp-phi)*ABS( uAvg(1) - GLM_ch ) + phi*LambdaMax
      Dmatrix(7) = (1.0_dp-phi)*ABS( uAvg(1) - cs ) + phi*LambdaMax
      Dmatrix(8) = (1.0_dp-phi)*ABS( uAvg(1) - ca ) + phi*LambdaMax
      Dmatrix(9) = (1.0_dp-phi)*ABS( uAvg(1) - cf ) + phi*LambdaMax
      
    ! Pure 9waves solver (no LLF)
!#      
!#      Dmatrix(1) = ABS( uAvg(1) + cf ) ! + fast magnetoacoustic wave
!#      Dmatrix(2) = ABS( uAvg(1) + ca ) ! + Alfven wave
!#      Dmatrix(3) = ABS( uAvg(1) + cs ) ! + slow magnetoacoustic wave
!#      Dmatrix(4) = ABS( uAvg(1) + GLM_ch ) ! + GLM wave
!#      Dmatrix(5) = ABS( uAvg(1)      ) ! / Entropy wave
!#      Dmatrix(6) = ABS( uAvg(1) - GLM_ch ) ! - GLM wave
!#      Dmatrix(7) = ABS( uAvg(1) - cs ) ! - slow magnetoacoustic wave
!#      Dmatrix(8) = ABS( uAvg(1) - ca ) ! - Alfven wave
!#      Dmatrix(9) = ABS( uAvg(1) - cf ) ! - fast magnetoacoustic wave
    
    ! LLF solver
!#      LambdaMax = MAX(ABS( uAvg(1) + cf ),ABS( uAvg(1) - cf ))
!#      
!#      Dmatrix(1) = LambdaMax
!#      Dmatrix(2) = LambdaMax
!#      Dmatrix(3) = LambdaMax
!#      Dmatrix(4) = LambdaMax
!#      Dmatrix(5) = LambdaMax
!#      Dmatrix(6) = LambdaMax
!#      Dmatrix(7) = LambdaMax
!#      Dmatrix(8) = LambdaMax
!#      Dmatrix(9) = LambdaMax
    
    ! Diagonal scaling matrix as described in Winters et al., eq. (4.15)
    Tmatrix = 0.0
    Tmatrix(1) = 0.5_dp/kappa/rhoLN ! + f
    Tmatrix(2) = 0.25_dp/betaAvg/rhoLN/rhoLN ! + a
    Tmatrix(3) = Tmatrix(1) ! + s
    Tmatrix(4) = 0.25_dp / betaAvg ! + GLM
    Tmatrix(5) = rhoLN*kappaM1/kappa ! E
    Tmatrix(6) = Tmatrix(4) ! - GLM
    Tmatrix(7) = Tmatrix(1) ! - s
    Tmatrix(8) = Tmatrix(2) ! - a
    Tmatrix(9) = Tmatrix(1) ! - f

    ! Scale D matrix
    Dmatrix = Dmatrix*Tmatrix

end ASSOCIATE
end subroutine EntropyStable9WaveFlux_VolFluxAndDissipMatrices  

!==================================================================================================================================
!> Transformation from conservative variables U to entropy vector, dS/dU, S = -rho*s/(kappa-1), s=ln(p)-kappa*ln(rho)
!==================================================================================================================================
PURE FUNCTION ConsToEntropy(cons) RESULT(Entropy)

!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL(dp),DIMENSION(PP_nVar),INTENT(IN)  :: cons    !< vector of conservative variables
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL(dp),DIMENSION(PP_nVar)             :: entropy !< vector of entropy variables
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL(dp)                                :: srho,u,v,w,v2s2,rho_sp,s
!==================================================================================================================================
srho   = 1.0_dp/cons(1)
u      = cons(2)*srho
v      = cons(3)*srho
w      = cons(4)*srho
v2s2   = 0.5_dp*(u*u+v*v+w*w)
rho_sp = cons(1)/(KappaM1*(cons(5)-cons(1)*v2s2-s2mu_0*SUM(cons(6:PP_nVar)*cons(6:PP_nVar)))) ! pressure includes psi^2 if PP_nVar=9
!s      = LOG(p) - kappa*LOG(cons(1))
s      = - LOG(rho_sp*(cons(1)**kappaM1))

! Convert to entropy variables
entropy(1)         =  (kappa-s)*skappaM1 - rho_sp*v2s2  !(kappa-s)/(kappa-1)-beta*|v|^2
entropy(2)         =  rho_sp*u                  ! 2*beta*u
entropy(3)         =  rho_sp*v                  ! 2*beta*v
entropy(4)         =  rho_sp*w                  ! 2*beta*w
entropy(5)         = -rho_sp                    !-2*beta
entropy(6:PP_nVar) =  rho_sp*cons(6:PP_nVar)    ! 2*beta*B +2*beta*psi

END FUNCTION ConsToEntropy

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
PURE FUNCTION LN_MEAN(aL,aR)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL(dp),INTENT(IN) :: aL  !< left value
REAL(dp),INTENT(IN) :: aR  !< right value
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL(dp)            :: LN_MEAN  !< result
!----------------------------------------------------------------------------------------------------------------------------------
! LOCaR VaLIABLES
REAL(dp)           :: Xi,u
REAL(dp),PARAMETER :: eps=1.0E-4  ! tolerance for f^2, such that switch is smooth in double precision 
!==================================================================================================================================
Xi = aR/aL
u=(Xi*(Xi-2.)+1.)/(Xi*(Xi+2.)+1.) !u=f^2, f=(aR-aL)/(aR+aL)=(xi-1)/(xi+1)
LN_MEAN=MERGE((aL+aR)*52.5_dp/(105_dp + u*(35_dp + u*(21_dp +u*15_dp))), & !u <eps (test true)
              (aR-aL)/LOG(Xi)                                         , & !u>=eps (test false)
              (u.LT.eps)                                              )   !test
END FUNCTION LN_MEAN

end module
