module riemann_mod

use constants_mod

contains

pure SUBROUTINE RiemannSolverByRoe(ConsL,ConsR,Flux)

    use equations_mod, only: xflux
    USE equations_mod, ONLY: Cons2Prim
    USE equations_mod, ONLY: Cons2Prim

    USE equations_mod, ONLY: WaveSpeeds1D
    USE setup_config_mod, ONLY: Kappa
    USE setup_config_mod, ONLY: mu0

    real(dp), parameter :: KappaM1 = kappa-1.0_dp
    real(dp), parameter :: KappaM2 = kappa-2.0_dp
    real(dp), parameter :: sKappaM1 = 1.0_dp/(kappa-1.0_dp)

    real(dp), parameter :: smu_0 = 1.0_dp/mu0
    real(dp), parameter :: s2mu_0 = 0.5_dp/mu0

    ! INPUT VARIABLES
    REAL(dp),INTENT(IN)  :: ConsL(1:N_VARS) !<  left conservative state  
    REAL(dp),INTENT(IN)  :: ConsR(1:N_VARS) !< right conservative state

    !----------------------------------------------------------------------------------------------------------------------------------
    ! OUTPUT VARIABLES
    REAL(dp),INTENT(OUT) :: Flux(1:N_VARS) !<numerical flux
    !----------------------------------------------------------------------------------------------------------------------------------
    ! LOCAL VARIABLES
    !----------------------------------------------------------------------------------------------------------------------------------
    INTEGER          :: i

    REAL(dp)             :: FluxL(N_VARS), FluxR(N_VARS)
    REAL(dp)             :: PrimL(N_VARS), PrimR(N_VARS)

    REAL(dp)             :: ptot_L, ptot_R
    REAL(dp)             :: SqrtRho_L,SqrtRho_R,sSqrtRoe_LR
    REAL(dp)             :: Roe_L,Roe_R
    REAL(dp)             :: Rho_Roe,sRho_Roe,SqrtRho_Roe,sSqrtRho_Roe,RoeVel(3),RoeB(3),RoeVel2,RoeB2
    REAL(dp)             :: H_L,H_R,RoeH
    REAL(dp)             :: XX, va2_Roe,c2_Roe,sc2_Roe,ca2_Roe,cs2_Roe,cf2_Roe,astar_Roe
    REAL(dp)             :: c_Roe,ca_Roe,cs_Roe,cf_Roe
    REAL(dp)             :: lambda(1:7), eta(1:7), EVr(1:8,1:7)
    REAL(dp)             :: delta_rho,delta_v(3),delta_B(3),delta_p
    REAL(dp)             :: beta_y,beta_z,beta_v,beta_dv,beta_dB,BtMag,sBx 
    REAL(dp)             :: alpha_f,alpha_s,scf2_cs2,as_cs_s,as_c_ssqRho,af_cf_s,af_c_ssqRho
    REAL(dp)             :: efix_delta,sefix_delta

    efix_delta=1.0E-06_dp
    sefix_delta=1.0E+06_dp

    PrimL(:) = cons2prim(ConsL(:))
    PrimR(:) = cons2prim(ConsR(:))

    !--------------------------------------------------------------------------------
    !use Roe mean wavespeeds from  Roe meanvalues 
    !   (paper by Cargo & Gallice: "Roe Matrices for Ideal MHD and ...",1997)

    SqrtRho_L   = SQRT(PrimL(1))
    SqrtRho_R   = SQRT(PrimR(1))
    sSqrtRoe_LR = 1./(SqrtRho_L+SqrtRho_R)
    Roe_L       = SqrtRho_L*sSqrtRoe_LR
    Roe_R       = SqrtRho_R*sSqrtRoe_LR
    ptot_L      = PrimL(5)+s2mu_0*SUM(PrimL(6:8)**2) !Total presssure!
    ptot_R      = PrimR(5)+s2mu_0*SUM(PrimR(6:8)**2) !Total presssure!

    Rho_Roe     = (SqrtRho_L*SqrtRho_R)
    sRho_Roe    = 1./Rho_Roe
    SqrtRho_Roe = SQRT(Rho_Roe)
    sSqrtRho_Roe= 1./SqrtRho_Roe
    RoeVel(1:3) = Roe_L*PrimL(2:4)+Roe_R*PrimR(2:4)
    RoeVel2     = SUM(RoeVel(:)**2)
    RoeB(1:3)   = Roe_L*PrimR(6:8)+Roe_R*PrimL(6:8)
    RoeB2       = SUM(RoeB(:)**2)

    H_L  = (ConsL(5)+ptot_L)/PrimL(1)
    H_R  = (ConsR(5)+ptot_R)/PrimR(1)
    RoeH = Roe_L*H_L+Roe_R*H_R

    XX   = 0.5*SUM((PrimL(6:8)-PrimR(6:8))**2)*(sSqrtRoe_LR**2)

    va2_Roe   = RoeB2*smu_0*sRho_Roe
    c2_Roe    = (2.-Kappa)*XX+ KappaM1*(RoeH-0.5*RoeVel2-va2_Roe)
    sc2_Roe   = 1./c2_Roe
    c_Roe     = SQRT(c2_Roe) 
    ca2_Roe   = RoeB(1)**2*smu_0*sRho_Roe
    ca_Roe    = SQRT(ca2_Roe) 
    astar_Roe = SQRT((c2_Roe+va2_Roe)**2-4.*c2_Roe*ca2_Roe) 
    cf2_Roe   = 0.5*(c2_Roe+va2_Roe+astar_Roe)
    cf_Roe    = SQRT(cf2_Roe)
    cs2_Roe   = 0.5*(c2_Roe+va2_Roe-astar_Roe)
    cs_Roe    = SQRT(cs2_Roe)
    !--------------------------------------------------------------------------------

    FluxL(1:N_VARS) = xFlux(ConsL(1:N_VARS))
    FluxR(1:N_VARS) = xFlux(ConsR(1:N_VARS))

    ! the following implementation is based on the astrophysics code Pluto (roe.c), which implements the Cargo & Gallice Roe solver
    sBx  = SIGN(1.0_dp,RoeB(1))
    Btmag=SQRT(RoeB(2)*RoeB(2)+RoeB(3)*RoeB(3))
    IF(Btmag.GT. 1.0E-9)THEN
      beta_y=RoeB(2)/BtMag
      beta_z=RoeB(3)/BtMag
    ELSE
      beta_y=SQRT(0.5)
      beta_z=beta_y
    END IF

    IF(cf2_Roe .EQ. cs2_Roe) THEN
      alpha_f  = 1.
      alpha_s  = 0.
    ELSEIF(c2_Roe .LE. cs2_Roe)THEN
      alpha_f  = 0.
      alpha_s  = 1.
    ELSEIF(cf2_Roe .LE. c2_Roe)THEN
      alpha_f  = 1.
      alpha_s  = 0.
    ELSE
      scf2_cs2 = 1./(cf2_Roe - cs2_Roe)
      alpha_f  = (c2_Roe  - cs2_Roe)*scf2_cs2
      alpha_s  = (cf2_Roe -  c2_Roe)*scf2_cs2
      alpha_f  = SQRT(MAX(0.0_dp, alpha_f))
      alpha_s  = SQRT(MAX(0.0_dp, alpha_s))
    END IF

    ! calculate differences in primitive variables
    Delta_rho    = PrimR(1) - PrimL(1)
    Delta_v(1:3) = PrimR(2:4) - PrimL(2:4)
    Delta_B(1:3) = PrimR(6:8) - PrimL(6:8)
    Delta_p   = KappaM1*( (0.5*RoeVel2-XX)*delta_rho -SUM(RoeVel*(ConsR(2:4)-ConsL(2:4)))+(ConsR(5)-ConsL(5))-smu_0*SUM(RoeB*delta_B)) 

    ! mean eigenvalues
    lambda(1)    = RoeVel(1)-cf_Roe
    lambda(2)    = RoeVel(1)-ca_Roe
    lambda(3)    = RoeVel(1)-cs_Roe
    lambda(4)    = RoeVel(1)
    lambda(5)    = RoeVel(1)+cs_Roe
    lambda(6)    = RoeVel(1)+ca_Roe
    lambda(7)    = RoeVel(1)+cf_Roe

    ! mean eigenvectors
    ! u-c_f
    beta_v  = beta_y* RoeVel(2) + beta_z* RoeVel(3)
    beta_dv = beta_y*delta_v(2) + beta_z*delta_v(3)
    beta_dB = beta_y*delta_B(2) + beta_z*delta_B(3)

    as_cs_s    = alpha_s*cs_Roe*sBx
    as_c_ssqRho= alpha_s*c_Roe*sSqrtRho_Roe

    EVr(1,1) = alpha_f
    EVr(2,1) = alpha_f*lambda(1)
    EVr(3,1) = alpha_f*RoeVel(2)+as_cs_s*beta_y
    EVr(4,1) = alpha_f*RoeVel(3)+as_cs_s*beta_z
    EVr(5,1) = alpha_f*(RoeH-va2_Roe-RoeVel(1)*cf_Roe)+as_cs_s*beta_v+as_c_ssqRho*Btmag
    EVr(6,1) = 0.
    EVr(7,1) = as_c_ssqRho*beta_y 
    EVr(8,1) = as_c_ssqRho*beta_z


    eta(1) = 0.5*sc2_Roe*(  alpha_f*(XX*delta_rho+delta_p)  & 
                          + Rho_Roe*(as_cs_s*beta_dv-alpha_f*cf_Roe*delta_v(1) + as_c_ssqRho*beta_dB))

    ! u+c_f

    EVr(1,7) = alpha_f
    EVr(2,7) = alpha_f*lambda(7)
    EVr(3,7) = alpha_f*RoeVel(2)-as_cs_s*beta_y
    EVr(4,7) = alpha_f*RoeVel(3)-as_cs_s*beta_z
    EVr(5,7) = alpha_f*(RoeH-va2_Roe+RoeVel(1)*cf_Roe)-as_cs_s*beta_v+as_c_ssqRho*Btmag
    EVr(6,7) = 0.
    EVr(7,7) = EVr(7,1)
    EVr(8,7) = EVr(8,1)


    eta(7) = 0.5*sc2_Roe*(  alpha_f*(XX*delta_rho+delta_p)  &
                          - Rho_Roe*(as_cs_s*beta_dv-alpha_f*cf_Roe*delta_v(1) - as_c_ssqRho*beta_dB))

    ! u -c_s
    af_cf_s=alpha_f*cf_Roe*sBx
    af_c_ssqRho=alpha_f*c_Roe*sSqrtRho_Roe

    EVr(1,3) = alpha_s
    EVr(2,3) = alpha_s*lambda(3)
    EVr(3,3) = alpha_s*RoeVel(2)-af_cf_s*beta_y
    EVr(4,3) = alpha_s*RoeVel(3)-af_cf_s*beta_z
    EVr(5,3) = alpha_s*(RoeH-va2_Roe-RoeVel(1)*cs_Roe)-af_cf_s*beta_v-af_c_ssqRho*Btmag
    EVr(6,3) = 0.
    EVr(7,3) = -af_c_ssqRho*beta_y 
    EVr(8,3) = -af_c_ssqRho*beta_z

    eta(3) = 0.5*sc2_Roe*(  alpha_s*(XX*delta_rho+delta_p) &
                          - Rho_Roe*(af_cf_s*beta_dv+alpha_s*cs_Roe*delta_v(1) + af_c_ssqRho*beta_dB))

    ! u +c_s

    EVr(1,5) = alpha_s
    EVr(2,5) = alpha_s*lambda(5)
    EVr(3,5) = alpha_s*RoeVel(2)+af_cf_s*beta_y
    EVr(4,5) = alpha_s*RoeVel(3)+af_cf_s*beta_z
    EVr(5,5) = alpha_s*(RoeH-va2_Roe+RoeVel(1)*cs_Roe)+af_cf_s*beta_v-af_c_ssqRho*Btmag
    EVr(6,5) = 0.
    EVr(7,5) = EVr(7,3)
    EVr(8,5) = EVr(8,3) 

    eta(5) = 0.5*sc2_Roe*(  alpha_s*(XX*delta_rho+delta_p) &
                          + Rho_Roe*(af_cf_s*beta_dv+alpha_s*cs_Roe*delta_v(1) - af_c_ssqRho*beta_dB))

    ! u-c_a
    beta_v  = beta_z* RoeVel(2) - beta_y* RoeVel(3) !overwrite!
    beta_dv = beta_z*delta_v(2) - beta_y*delta_v(3) !overwrite!
    beta_dB = beta_z*delta_B(2) - beta_y*delta_B(3) !overwrite!

    EVr(1,2)=0.
    EVr(2,2)=0.
    EVr(3,2)=-Rho_Roe*beta_z
    EVr(4,2)= Rho_Roe*beta_y
    EVr(5,2)=-Rho_Roe*beta_v
    EVr(6,2)=0. 
    EVr(7,2)=-sBx*SqrtRho_roe*beta_z 
    EVr(8,2)= sBx*SqrtRho_roe*beta_y 

    eta(2) = 0.5*( -beta_dv  - sBx*sSqrtRho_Roe*beta_dB)

    ! u+c_a

    EVr(1,6)=0.
    EVr(2,6)=0.
    EVr(3,6)=-EVr(3,2)
    EVr(4,6)=-EVr(4,2)
    EVr(5,6)=-EVr(5,2)
    EVr(6,6)=0. 
    EVr(7,6)=EVr(7,2)
    EVr(8,6)=EVr(8,2)

    eta(6) = 0.5*( beta_dv  - sBx*sSqrtRho_Roe*beta_dB)

    !u 

    Evr(1,4)   = 1.
    Evr(2,4)   = RoeVel(1)
    Evr(3,4)   = RoeVel(2)
    Evr(4,4)   = RoeVel(3)
    Evr(5,4)   = 0.5*RoeVel2+ KappaM2*sKappaM1*XX
    Evr(6:8,4) = 0.

    eta(4) =  sc2_Roe*((c2_Roe - XX)*delta_rho - delta_p)

    ! entropy fix

    IF(ABS(lambda(1)) .LT. 0.5*efix_delta) lambda(1)= lambda(1)*lambda(1)*sefix_delta+0.25*efix_delta
    IF(ABS(lambda(7)) .LT. 0.5*efix_delta) lambda(7)= lambda(7)*lambda(7)*sefix_delta+0.25*efix_delta
    IF(ABS(lambda(3)) .LT. 0.5*efix_delta) lambda(3)= lambda(3)*lambda(3)*sefix_delta+0.25*efix_delta
    IF(ABS(lambda(6)) .LT. 0.5*efix_delta) lambda(6)= lambda(6)*lambda(6)*sefix_delta+0.25*efix_delta

    !CHECK Roe Matrix condition FR-FL = A*(UR-UL) = R*lambda*eta
    !Flux(:)=FluxR-FluxL
    !DO i=1,7
    !  Flux(1:8)=Flux(1:8) - eta(i)*lambda(i)*EVr(1:8,i)
    !END DO
    !DO i=1,8
    !  IF(ABS(Flux(i))>1.0E-08)THEN
    !    WRITE(*,*)'Roe matrix not satisfied, var=',i,'Flux(i)=',Flux(i)
    !  END IF 
    !END DO
    !IF(MAXVAL(ABS(Flux))>1.0E-08) STOP


    ! assemble Roe flux
    Flux(:)=FluxL(:)+FluxR(:)
    DO i=1,7
      Flux(1:8)=Flux(1:8) - eta(i)*ABS(lambda(i))*EVr(1:8,i)
    END DO
    Flux(:)=0.5*Flux(:)

end subroutine

pure function xriemann(uL, uR) result(sf)

    real(dp), intent(in) :: uL(N_VARS), uR(N_VARS)
    real(dp)             :: sf(N_VARS)
    
    call RiemannSolverByRoe(uL,uR,sf)

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
