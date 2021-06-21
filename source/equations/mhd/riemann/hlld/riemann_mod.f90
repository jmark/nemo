# Error The HLLD Riemann solver does not work with GLM divergence cleaning. Do not use this solver.

module riemann_mod

use constants_mod

contains

!==================================================================================================================================
!> HLLD solver following "A multi-state HLL approximate Riemann solver for ideal magnetohydrodynamics"
!> Takahiro Miyoshi  Kanya Kusano, JCP 208 (2005). Follow implementation of Elise Estibals
!> Input state already rotated to normal system, and 
!> ONLY WORKS FOR mu_0=1!!!
!==================================================================================================================================

pure SUBROUTINE RiemannSolverByHLLD(ConsL_in,ConsR_in,Flux)

    use equations_mod, only: xflux
    USE equations_mod, ONLY: Cons2Prim
    USE equations_mod, ONLY: prim2cons

    USE equations_mod, ONLY: Fastest_Signal_Speed
    USE setup_config_mod, ONLY: Kappa
    USE setup_config_mod, ONLY: mu0

    real(dp), parameter :: KappaM1 = 1.0_dp/kappa
    real(dp), parameter :: KappaM2 = kappa-2.0_dp
    real(dp), parameter :: sKappaM1 = 1.0_dp/(kappa-1.0_dp)

    real(dp), parameter :: smu_0 = 1.0_dp/mu0
    real(dp), parameter :: s2mu_0 = 0.5_dp/mu0

    REAL(dp),INTENT(IN)  :: ConsL_in(1:N_VARS) !<  left conservative state  
    REAL(dp),INTENT(IN)  :: ConsR_in(1:N_VARS) !< right conservative state

    REAL(dp),INTENT(OUT) :: Flux(1:N_VARS) !<numerical flux

    REAL(dp)                 :: ConsL(1:N_VARS), ConsR(1:N_VARS)
    REAL(dp)                 :: PrimL(1:N_VARS), PrimR(1:N_VARS)
    REAL(dp)                 :: FluxL(N_VARS), FluxR(N_VARS)
    REAL(dp)                 :: SL, SR, SM
    REAL(dp)                 :: SLS, SRS
    REAL(dp)                 :: FastestWave 
    REAL(dp)                 :: sSL_SM,sSR_SM
    REAL(dp)                 :: cf_L,cf_R
    REAL(dp)                 :: temp,stemp
    REAL(dp)                 :: signBn,Bn,ptotS 
    REAL(dp)                 :: rhoLS,sqrtRhoLS,rhoLSS,uB,uBS
    REAL(dp)                 :: rhoRS,sqrtRhoRS,rhoRSS
    REAL(dp)                 :: pmagL,ptotL,eLS,eLSS !,rhoL,pL,unL,eL
    REAL(dp)                 :: pmagR,ptotR,eRS,eRSS !,rhoR,pR,unR,eR
    REAL(dp),DIMENSION(1:3)  :: uLS,utL,utLS  !,uL
    REAL(dp),DIMENSION(1:3)  :: uRS,utR,utRS  !,uR
    REAL(dp),DIMENSION(1:3)  :: BLS,BtL,BtLS  !,BL
    REAL(dp),DIMENSION(1:3)  :: BRS,BtR,BtRS  !,BR
    REAL(dp),DIMENSION(1:3)  :: uSS,BSS,utSS,BtSS
    REAL(dp),DIMENSION(1:N_VARS)  :: U_LS, U_LSS
    REAL(dp),DIMENSION(1:N_VARS)  :: U_RS, U_RSS
    REAL(dp), PARAMETER      :: hlld_eps = 1.0e-8_dp
    REAL(dp), PARAMETER      :: hlld_small_eps = 1.0e-12_dp

    !==================================================================================================================================

    ! CALL ConsToPrim(PrimL(:),ConsL_in(:))
    ! CALL ConsToPrim(PrimR(:),ConsR_in(:))

    primL = Cons2Prim(consL_in)
    primR = Cons2Prim(consR_in)

    ! make B_n continuous 
    Bn=0.5*(PrimL(6)+PrimR(6))
    PrimR(6)=Bn
    PrimL(6)=Bn

    ! CALL PrimToCons(PrimL(:),ConsL(:))
    ! CALL PrimToCons(PrimR(:),ConsR(:))

    consL = prim2cons(primL)
    consR = prim2cons(primR)

    !--------------------------------------------------------------------------------

    CALL Fastest_Signal_Speed(PrimL,cf_L)
    CALL Fastest_Signal_Speed(PrimR,cf_R)

    FastestWave=MAX(cf_L,cf_R)

    SL = MIN(PrimL(2),PrimR(2))-FastestWave


    IF (SL .GT. 0.0)  THEN
      !CALL EvalAdvectionFlux1D(ConsL,Flux)
        flux = xflux(consL)
      RETURN
    END IF

    SR = MAX(PrimL(2),PrimR(2))+FastestWave

    IF (SR .LT. 0.0) THEN
      ! CALL EvalAdvectionFlux1D(ConsR,Flux)
        flux = xflux(consR)
      RETURN
    END IF

    ! NOW SL<0 and SR>0
    ! CALL EvalAdvectionFlux1D(ConsL,FluxL)
    ! CALL EvalAdvectionFlux1D(ConsR,FluxR)

    fluxL = xflux(consL)
    fluxR = xflux(consR)

    !-------------------------------------------------

    !Left and Right variables
    ASSOCIATE( rhoL => PrimL(1)          , rhoR => PrimR(1)          ,  &  
               unL  => PrimL(2)          , unR  => PrimR(2)          ,  & 
               pL   => PrimL(5)          , pR   => PrimR(5)          ,  & 
               eL   => ConsL(5)          , eR   => ConsR(5)          ,  & 
               uL   => PrimL(2:4)        , uR   => PrimR(2:4)        ,  & 
               BL   => PrimL(6:8)        , BR   => PrimR(6:8)           )

    utL(1:3)  = (/0.0_dp,uL(2),uL(3)/) ; utR(1:3)  = (/0.0_dp,uR(2),uR(3)/) 
    BtL(1:3)  = (/0.0_dp,BL(2),BL(3)/) ; BtR(1:3)  = (/0.0_dp,BR(2),BR(3)/) 

    !magnetic and total pressures
    pmagL = 0.5_dp*(BL(1)**2 + BL(2)**2 + BL(3)**2) ; ptotL = pL + pmagL
    pmagR = 0.5_dp*(BR(1)**2 + BR(2)**2 + BR(3)**2) ; ptotR = pR + pmagR

    !Compute middle slope
    SM  = (rhoR*unR*(SR-unR) - rhoL * unL*(SL-unL) - ptotR + ptotL) / &
         & (rhoR*(SR-unR) - rhoL*(SL-unL))

    sSL_SM=1.0_dp/(SL-SM)
    sSR_SM=1.0_dp/(SR-SM)
    rhoLS = rhoL * (SL-unL) * sSL_SM
    rhoRS = rhoR * (SR-unR) * sSR_SM 

    sqrtRhoLS=SQRT(rhoLS)
    sqrtRhoRS=SQRT(rhoRS)

    SLS = SM - ABS(Bn) / SqrtRhoLS
    SRS = SM + ABS(Bn) / SqrtRhoRS

    !Compute intermediate states
    ptotS = ptotL + rhoL * (SL - unL) * (SM - unL)

    !Left * state
    temp = rhoL * (SL-unL) * (SL - SM) - Bn*Bn
    IF (ABS(temp)<hlld_eps) THEN
      utLS = utL
      BtLS = BtL
    ELSE
      stemp=1./temp
      utLS = utL - (Bn*(SM-unL)*stemp)*BtL
      BtLS = ((rhoL*(SL-unL)*(SL-unL)-Bn*Bn)*stemp)* BtL
    END IF
    BLS = (/Bn,0.0_dp,0.0_dp/) + BtLS
    uLS = (/SM,0.0_dp,0.0_dp/) + utLS

    uBS = uLS(1)*BLS(1) + uLS(2)*BLS(2) + uLS(3)*BLS(3)
    uB  = uL(1) *BL(1)  + uL(2) *BL(2)  + uL(3) *BL(3)

    eLS = ( eL*(SL-unL) +  ptotS*SM - ptotL*unL + Bn*(uB-uBS) ) *sSL_SM

    U_LS(1)   = rhoLS       
    U_LS(2:4) = rhoLS * uLS 
    U_LS(5)   = eLS         
    U_LS(6:8) = BLS         
    U_LS(9) = ConsL(9)

    !Right * state
    temp = rhoR * (SR-unR) * (SR - SM) - Bn*Bn
    IF (ABS(temp)<hlld_eps) THEN
       utRS = utR
       BtRS = BtR
    ELSE
       stemp=1./temp
       utRS = utR - (Bn*(SM-unR)*stemp)*BtR
       BtRS = ((rhoR*(SR-unR)*(SR-unR)-Bn*Bn)*stemp) * BtR
    END IF
    BRS = (/Bn,0.0_dp,0.0_dp/) + BtRS
    uRS = (/SM,0.0_dp,0.0_dp/) + utRS

    uBS = uRS(1)*BRS(1) + uRS(2)*BRS(2) + uRS(3)*BRS(3)
    uB  = uR(1) *BR(1)  + uR(2) *BR(2)  + uR(3) *BR(3)

    eRS = (eR*(SR-unR) - ptotR*unR + ptotS*SM  + Bn*(uB-uBS) ) *sSR_SM

    U_RS(1)   = rhoRS
    U_RS(2:4) = rhoRS * uRS
    U_RS(5)   = eRS
    U_RS(6:8) = BRS
    U_RS(9) = ConsR(9)

    IF (ABS(Bn)<hlld_small_eps) THEN
       !If Bn=0 four states reduce to two
       U_LSS = U_LS
       U_RSS = U_RS
    ELSE
       rhoLSS = rhoLS
       rhoRSS = rhoRS

       temp = SqrtRhoLS + SqrtRhoRS

       IF (Bn<0.0) THEN
          signBn = -1.0_dp
       ELSE
          signBn =  1.0_dp
       END IF
       stemp=1./temp 
       utSS = (SqrtRhoLS*utLS + SqrtRhoRS*utRS + signBn * (BtRS - BtLS)) * stemp
       BtSS = (SqrtRhoLS*BtRS + SqrtRhoRS*BtLS + signBn * (utRS - utLS) * sqrtRhoLS*sqrtRhoRS) *stemp
       
       uSS  = (/SM,0.0_dp,0.0_dp/) + utSS
       BSS  = (/Bn,0.0_dp,0.0_dp/) + BtSS
       uBS  = utSS(1)*BtSS(1) + utSS(2)*BtSS(2) + utSS(3)*BtSS(3)
       
       uB   = utLS(1)*BtLS(1) + utLS(2)*BtLS(2) + utLS(3)*BLS(3)
       eLSS = eLS - signBn * sqrtRhoLS * (uB - uBS)
      
       uB   = uRS(1)*BtRS(1) + utRS(2)*BtRS(2) + utRS(3)*BtRS(3)
       eRSS = eRS + signBn * sqrtRhoRS * (uB - uBS)

       U_LSS(1)   = rhoLSS          ; U_RSS(1)   = rhoRSS
       U_LSS(2:4) = rhoLSS*uSS      ; U_RSS(2:4) = rhoRSS*uSS
       U_LSS(5)   = eLSS            ; U_RSS(5)   = eRSS
       U_LSS(6:8) = BSS             ; U_RSS(6:8) = BSS
       U_LSS(9) = ConsL(9)          ; U_RSS(9)   = consR(9)
    END IF !|Bn|<0

    END ASSOCIATE !rhoL,rhoR,...

    !Fluxes
    IF (0.0_dp<=SLS) THEN
      !U=ULS
      Flux(:) = FluxL(:) + SL*(U_LS(:) - ConsL(:))
    ELSE IF (0.0_dp<=SM) THEN
      !U=ULSS
      Flux(:) = FluxL(:) + SLS*U_LSS(:) - (SLS-SL)*U_LS(:) - SL*ConsL(:)
    ELSE IF (0.0_dp<=SRS) THEN
      !U=URSS
      Flux(:) = FluxR(:) + SRS*U_RSS(:) - (SRS-SR)*U_RS(:) - SR*ConsR(:)
    ELSE
      !U=URS
      Flux(:) = FluxR(:) + SR*(U_RS(:) - ConsR(:))
    END IF

END SUBROUTINE RiemannSolverByHLLD

pure function xriemann(uL, uR) result(sf)

    real(dp), intent(in) :: uL(N_VARS), uR(N_VARS)
    real(dp)             :: sf(N_VARS)
    
    call RiemannSolverByHLLD(uL,uR,sf)

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
