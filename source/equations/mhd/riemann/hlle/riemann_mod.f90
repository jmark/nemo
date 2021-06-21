# Error The HLLE Riemann solver does not work with GLM divergence cleaning. Do not use this solver.

module riemann_mod

use constants_mod
use setup_config_mod, only: kappa
use setup_config_mod, only: mu0
use equations_share_mod, only: glm_ch

real(dp), parameter :: KappaM1 = kappa-1.0_dp
real(dp), parameter :: KappaM2 = kappa-2.0_dp

real(dp), parameter :: smu_0 = 1.0_dp/mu0
real(dp), parameter :: s2mu_0 = 0.5_dp/mu0

contains

pure SUBROUTINE EvalAdvectionFlux1D(U_Face,F_Face)

    ! MODULES
    ! USE MOD_Equation_Vars,ONLY:KappaM1,KappaM2,smu_0,s2mu_0
    ! #ifdef PP_GLM
    ! USE MOD_Equation_vars,ONLY:GLM_ch
    ! #endif 

    !----------------------------------------------------------------------------------------------------------------------------------
    ! INPUT VARIABLES
    REAL(dp),INTENT(IN)     :: U_Face(N_VARS)    !< input conservative state
    !----------------------------------------------------------------------------------------------------------------------------------
    ! OUTPUT VARIABLES
    REAL(dp),INTENT(OUT)    :: F_Face(N_VARS)    !< F_Face(iVar), Cartesian flux in "x" direction
    !----------------------------------------------------------------------------------------------------------------------------------
    ! LOCAL VARIABLES 
    REAL(dp)                :: srho               ! reciprocal values for density and kappa/Pr
    REAL(dp)                :: v1,v2,v3,ptilde    ! velocity and pressure(including magnetic pressure
    REAL(dp)                :: vb                 ! magnetic field, bb2=|bvec|^2
    !==================================================================================================================================
    ! auxiliary variables
    srho = 1. / U_Face(1) ! 1/rho
    v1   = U_Face(2)*srho ! u
    v2   = U_Face(3)*srho ! v
    v3   = U_Face(4)*srho ! w
    ASSOCIATE( b1 => U_Face(6), &
               b2 => U_Face(7), &
               b3 => U_Face(8), &
               Etotal => U_Face(5)-s2mu_0*U_Face(9)**2  &
    )

    vb   = (b1*v1+b2*v2+b3*v3)
    ! ptilde includes magnetic pressure
    ptilde = kappaM1*(Etotal-0.5*U_Face(1)*(v1*v1+v2*v2+v3*v3))-kappaM2*s2mu_0*(b1*b1+b2*b2+b3*b3)
    ! Advection fluxes x-direction
    F_Face(1)= U_Face(2)          ! rho*u
    F_Face(2)= U_Face(2)*v1+ptilde    -smu_0*b1*b1  ! rho*u²+p     -1/mu_0*b1*b1
    F_Face(3)= U_Face(2)*v2           -smu_0*b1*b2  ! rho*u*v      -1/mu_0*b1*b2
    F_Face(4)= U_Face(2)*v3           -smu_0*b1*b3  ! rho*u*w      -1/mu_0*b1*b3
    F_Face(5)=(Etotal + ptilde)*v1    -smu_0*b1*vb  ! (rho*e+p)*u  -1/mu_0*b1*(v dot B)
    F_Face(6)=0.
    F_Face(7)=v1*b2-b1*v2
    F_Face(8)=v1*b3-b1*v3

    F_Face(5)=F_Face(5) + smu_0 *GLM_ch*b1*U_Face(9)
    F_Face(6)=F_Face(6) +        GLM_ch   *U_Face(9)
    F_Face(9)=                   GLM_ch*b1
    END ASSOCIATE ! b1 => UFace(6) ...

END SUBROUTINE EvalAdvectionFlux1D

pure SUBROUTINE EvalAdvectionFlux1D_noglm(U_Face,F_Face)

    ! MODULES
    ! USE MOD_Equation_Vars,ONLY:KappaM1,KappaM2,smu_0,s2mu_0
    ! #ifdef PP_GLM
    ! USE MOD_Equation_vars,ONLY:GLM_ch
    ! #endif 

    !----------------------------------------------------------------------------------------------------------------------------------
    ! INPUT VARIABLES
    REAL(dp),INTENT(IN)     :: U_Face(N_VARS)    !< input conservative state
    !----------------------------------------------------------------------------------------------------------------------------------
    ! OUTPUT VARIABLES
    REAL(dp),INTENT(OUT)    :: F_Face(N_VARS)    !< F_Face(iVar), Cartesian flux in "x" direction
    !----------------------------------------------------------------------------------------------------------------------------------
    ! LOCAL VARIABLES 
    REAL(dp)                :: srho               ! reciprocal values for density and kappa/Pr
    REAL(dp)                :: v1,v2,v3,ptilde    ! velocity and pressure(including magnetic pressure
    REAL(dp)                :: vb                 ! magnetic field, bb2=|bvec|^2
    !==================================================================================================================================
    ! auxiliary variables

    F_Face = 0.0_dp

    srho = 1. / U_Face(1) ! 1/rho
    v1   = U_Face(2)*srho ! u
    v2   = U_Face(3)*srho ! v
    v3   = U_Face(4)*srho ! w
    ASSOCIATE( b1 => U_Face(6), &
               b2 => U_Face(7), &
               b3 => U_Face(8), &
               Etotal => U_Face(5)-s2mu_0*U_Face(9)**2  &
    )

    vb   = (b1*v1+b2*v2+b3*v3)
    ! ptilde includes magnetic pressure
    ptilde = kappaM1*(Etotal-0.5*U_Face(1)*(v1*v1+v2*v2+v3*v3))-kappaM2*s2mu_0*(b1*b1+b2*b2+b3*b3)
    ! Advection fluxes x-direction
    F_Face(1)= U_Face(2)          ! rho*u
    F_Face(2)= U_Face(2)*v1+ptilde    -smu_0*b1*b1  ! rho*u²+p     -1/mu_0*b1*b1
    F_Face(3)= U_Face(2)*v2           -smu_0*b1*b2  ! rho*u*v      -1/mu_0*b1*b2
    F_Face(4)= U_Face(2)*v3           -smu_0*b1*b3  ! rho*u*w      -1/mu_0*b1*b3
    F_Face(5)=(Etotal + ptilde)*v1    -smu_0*b1*vb  ! (rho*e+p)*u  -1/mu_0*b1*(v dot B)
    F_Face(6)=0.
    F_Face(7)=v1*b2-b1*v2
    F_Face(8)=v1*b3-b1*v3

    F_Face(5)=F_Face(5) + smu_0 *GLM_ch*b1*U_Face(9)
    F_Face(6)=F_Face(6) +        GLM_ch   *U_Face(9)
    ! F_Face(9)=                   GLM_ch*b1
    END ASSOCIATE ! b1 => UFace(6) ...

END SUBROUTINE EvalAdvectionFlux1D_noglm

pure SUBROUTINE EvalAdvectionFlux1D_glm(U_Face,F_Face)

    ! MODULES
    ! USE MOD_Equation_Vars,ONLY:KappaM1,KappaM2,smu_0,s2mu_0
    ! #ifdef PP_GLM
    ! USE MOD_Equation_vars,ONLY:GLM_ch
    ! #endif 

    !----------------------------------------------------------------------------------------------------------------------------------
    ! INPUT VARIABLES
    REAL(dp),INTENT(IN)     :: U_Face(N_VARS)    !< input conservative state
    !----------------------------------------------------------------------------------------------------------------------------------
    ! OUTPUT VARIABLES
    REAL(dp),INTENT(OUT)    :: F_Face(N_VARS)    !< F_Face(iVar), Cartesian flux in "x" direction
    !----------------------------------------------------------------------------------------------------------------------------------
    ! LOCAL VARIABLES 
    REAL(dp)                :: srho               ! reciprocal values for density and kappa/Pr
    REAL(dp)                :: v1,v2,v3,ptilde    ! velocity and pressure(including magnetic pressure
    REAL(dp)                :: vb                 ! magnetic field, bb2=|bvec|^2
    !==================================================================================================================================
    ! auxiliary variables
    ! srho = 1. / U_Face(1) ! 1/rho
    ! v1   = U_Face(2)*srho ! u
    ! v2   = U_Face(3)*srho ! v
    ! v3   = U_Face(4)*srho ! w
    ! ASSOCIATE( b1 => U_Face(6), &
    !            b2 => U_Face(7), &
    !            b3 => U_Face(8), &
    !            Etotal => U_Face(5)-s2mu_0*U_Face(9)**2  &
    ! )

    ASSOCIATE( b1 => U_Face(6))

    ! vb   = (b1*v1+b2*v2+b3*v3)
    ! ! ptilde includes magnetic pressure
    ! ptilde = kappaM1*(Etotal-0.5*U_Face(1)*(v1*v1+v2*v2+v3*v3))-kappaM2*s2mu_0*(b1*b1+b2*b2+b3*b3)
    ! ! Advection fluxes x-direction
    ! F_Face(1)= U_Face(2)          ! rho*u
    ! F_Face(2)= U_Face(2)*v1+ptilde    -smu_0*b1*b1  ! rho*u²+p     -1/mu_0*b1*b1
    ! F_Face(3)= U_Face(2)*v2           -smu_0*b1*b2  ! rho*u*v      -1/mu_0*b1*b2
    ! F_Face(4)= U_Face(2)*v3           -smu_0*b1*b3  ! rho*u*w      -1/mu_0*b1*b3
    ! F_Face(5)=(Etotal + ptilde)*v1    -smu_0*b1*vb  ! (rho*e+p)*u  -1/mu_0*b1*(v dot B)
    ! F_Face(6)=0.
    ! F_Face(7)=v1*b2-b1*v2
    ! F_Face(8)=v1*b3-b1*v3
        
    F_Face = 0.0_dp

    ! F_Face(5) = smu_0*GLM_ch*b1*U_Face(9)
    ! F_Face(6) =       GLM_ch   *U_Face(9)
    F_Face(9) =                 GLM_ch*b1

    END ASSOCIATE ! b1 => UFace(6) ...

END SUBROUTINE EvalAdvectionFlux1D_glm



PURE SUBROUTINE ConsToPrim(prim,cons)

    !----------------------------------------------------------------------------------------------------------------------------------
    ! INPUT VARIABLES
    REAL(dp),INTENT(IN)     :: cons(N_VARS) !< vector of conservative variables
    !----------------------------------------------------------------------------------------------------------------------------------
    ! OUTPUT VARIABLES
    REAL(dp),INTENT(OUT)    :: prim(N_VARS) !< vector of primitive variables
    !----------------------------------------------------------------------------------------------------------------------------------
    ! LOCAL VARIABLES 
    REAL(dp)                :: sRho    ! 1/Rho
    !==================================================================================================================================
    sRho=1./cons(1)
    ! rho
    prim(1)=cons(1)
    ! velocity
    prim(2:4)=cons(2:4)*sRho
    ! GAS PRESSURE (WITHOUT magnetic pressure)
    prim(5)=KappaM1*(cons(5)-0.5*SUM(cons(2:4)*prim(2:4))-s2mu_0*SUM(cons(6:N_VARS)*cons(6:N_VARS)))
    ! B,psi
    prim(6:N_VARS)=cons(6:N_VARS)

END SUBROUTINE ConsToPrim

!==================================================================================================================================
!> calculate fastest wave speed 
!==================================================================================================================================
pure SUBROUTINE FastestWave1D(Prim,cf)

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
    sRho=1./prim(1)
    c2=kappa*prim(5)*sRho
    ca2=Prim(6)*Prim(6)*sRho*smu_0 !Alfen wave speed
    va2=SUM(prim(7:8)*prim(7:8))*sRho*smu_0+ca2

    astar=SQRT((c2+va2)*(c2+va2)-4.*c2*ca2)
    cf=SQRT(0.5*(c2+va2+astar))

END SUBROUTINE FastestWave1D

!==================================================================================================================================
!> calculate fastest roe wave speed 
!> use Roe mean wavespeeds from  Roe meanvalues 
!>   (paper by Cargo & Gallice: "Roe Matrices for Ideal MHD and ...",1997)
!==================================================================================================================================
pure SUBROUTINE FastestWave1D_Roe(ConsL,ConsR,PrimL,PrimR,RoeVelx,cf_Roe)

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

END SUBROUTINE FastestWave1D_Roe

pure SUBROUTINE RiemannSolverByHLL(ConsL,ConsR,Flux)

    ! USE MOD_PreProc
    ! USE MOD_Flux,          ONLY: EvalAdvectionFlux1D
    ! USE MOD_Equation_vars, ONLY: FastestWave1D
    ! USE MOD_Equation_vars, ONLY: FastestWave1D_Roe
    ! USE MOD_Equation_vars, ONLY: ConsToPrim

    REAL(dp),INTENT(IN)  :: ConsL(1:N_VARS) !<  left conservative state  
    REAL(dp),INTENT(IN)  :: ConsR(1:N_VARS) !< right conservative state
    !----------------------------------------------------------------------------------------------------------------------------------
    ! OUTPUT VARIABLES
    REAL(dp),INTENT(OUT) :: Flux(1:N_VARS) !<numerical flux
    !----------------------------------------------------------------------------------------------------------------------------------
    ! LOCAL VARIABLES                                                                                                               !
    !----------------------------------------------------------------------------------------------------------------------------------
    REAL(dp)             :: FluxL(1:N_VARS), FluxR(1:N_VARS)
    REAL(dp)             :: PrimL(1:N_VARS), PrimR(1:N_VARS) !< left/right primitive state
    REAL(dp)             :: SL, SR
    REAL(dp)             :: cf_L,cf_R,cf_Roe,RoeVelx 

    REAL(dp) :: lmax

    !==================================================================================================================================

    CALL ConsToPrim(PrimL(:),ConsL(:))
    CALL ConsToPrim(PrimR(:),ConsR(:))

    CALL FastestWave1D(PrimL,cf_L)
    CALL FastestWave1D_Roe(ConsL,ConsR,PrimL,PrimR,RoeVelx,cf_Roe)
    SL = MIN(PrimL(2)-cf_L, RoeVelx-cf_Roe,0.0_dp)

    ! IF (SL .GT. 0.0)  THEN
    !   CALL EvalAdvectionFlux1D_noglm(ConsL,Flux)
    !   RETURN
    ! END IF

    CALL FastestWave1D(PrimR,cf_R)
    SR = MAX(RoeVelx+cf_Roe, PrimR(2)+cf_R,0.0_dp)

    ! IF (SR .LT. 0.0) THEN
    !   CALL EvalAdvectionFlux1D_noglm(ConsR,Flux)
    !   RETURN
    ! END IF

    ! NOW SL<0 and SR>0
    ! CALL EvalAdvectionFlux1D_noglm(ConsL,FluxL)
    ! CALL EvalAdvectionFlux1D_noglm(ConsR,FluxR)

    CALL EvalAdvectionFlux1D(ConsL,FluxL)
    CALL EvalAdvectionFlux1D(ConsR,FluxR)
    Flux = (SR*FluxL - SL*FluxR + SL*SR*(ConsR-ConsL))/(SR-SL)

    ! -- !

    ! CALL EvalAdvectionFlux1D_glm(ConsL,FluxL)
    ! CALL EvalAdvectionFlux1D_glm(ConsR,FluxR)

    ! lmax = max(abs(PrimL(2))+cf_L,abs(PrimR(2))+cf_R)

    ! ! Flux(5) = Flux(5) + 0.5_dp*(FluxL(5) + FluxR(5) + lmax*(ConsL(5)-ConsR(5)))
    ! ! Flux(6) = Flux(6) + 0.5_dp*(FluxL(6) + FluxR(6) + lmax*(ConsL(6)-ConsR(6)))
    ! Flux(9) = Flux(9) + 0.5_dp*(FluxL(9) + FluxR(9) + lmax*(ConsL(9)-ConsR(9)))

    ! -- !

    ! CALL EvalAdvectionFlux1D(ConsL,FluxL)
    ! CALL EvalAdvectionFlux1D(ConsR,FluxR)

    ! lmax = max(abs(PrimL(2))+cf_L,abs(PrimR(2))+cf_R)

    ! Flux = 0.5_dp*(FluxL + FluxR + lmax*(ConsL-ConsR))

END SUBROUTINE RiemannSolverByHLL

!! HLL Riemann solver
pure function xriemann(U_L,U_R) result(F)

    ! use equations_const_mod

    ! use equations_mod, only: flux => xflux
    ! use equations_mod, only: Fastest_Signal_Speed
    ! use equations_mod, only: Fastest_Signal_Speeds_Roe

    real(dp), intent(in) :: U_L(N_VARS), U_R(N_VARS)
    real(dp)             :: F(N_VARS)
    
    ! real(dp) :: F_L(N_VARS), F_R(N_VARS)

    ! real(dp) :: V_L,V_R !! velocities
    ! real(dp) :: vTilde, aTilde
    ! real(dp) :: c_L,c_R !! sound speed
    ! real(dp) :: s_L,s_R !! signal speeds

    call RiemannSolverByHLL(U_L,U_R,F)

    !! V_L = U_L(MOMX_VAR)/U_L(DENS_VAR)
    !! V_R = U_R(MOMX_VAR)/U_R(DENS_VAR)

    !! CALL Fastest_Signal_Speeds_Roe(U_L,U_R,vTilde,aTilde)

    !! CALL Fastest_Signal_Speed(U_L,c_L)
    !! CALL Fastest_Signal_Speed(U_R,c_R)

    !! S_L = MIN(vTilde-aTilde,v_L-c_L,0.0_dp)
    !! S_R = MAX(vTilde+aTilde,v_R+c_R,0.0_dp)

    !! !S_L = MIN(v_L-c_L,v_R-c_R,0.0_dp)
    !! !S_R = MAX(v_L+c_L,v_R+c_R,0.0_dp)

    !! !if (S_L .GE. 0.0_dp .and. S_R .GT. 0.0_dp) then
    !! ! if (S_L .GT. 0.0_dp) then

    !! !     F = flux(U_L)

    !! ! else if (S_R .LT. 0.0_dp) then

    !! !     F = flux(U_R)

    !! ! else

    !! !     F_L = flux(U_L)
    !! !     F_R = flux(U_R)

    !! !     F = (S_R*F_L - S_L*F_R + S_L*S_R*(U_R - U_L)) / (S_R - S_L)

    !! ! end if 

    !! F_L = flux(U_L)
    !! F_R = flux(U_R)

    !! F = (S_R*F_L - S_L*F_R + S_L*S_R*(U_R - U_L)) / (S_R - S_L)

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
