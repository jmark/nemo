!# define PP_NON_CONSERVATIVE_FLUX 1

module non_conservative_flux_mod

use constants_mod
use equations_const_mod

contains

pure function Non_Conservative_two_point_xFlux(UI,UL,UR) result(fl)

    use equations_const_mod
    use setup_config_mod, only: mu0

    real(dp), intent(in)    :: uI(N_VARS) !! inner state
    real(dp), intent(in)    :: uL(N_VARS) !! left state
    real(dp), intent(in)    :: uR(N_VARS) !! right state

    real(dp)                :: FL(N_VARS)

    real(dp) :: v(3),Bavg,Gavg

    real(dp)                :: mhdF(N_VARS)
    real(dp)                :: glmF(N_VARS)

    mhdF = 0.0_dp
    glmF = 0.0_dp

    Bavg = 0.5_dp*(UL(MAGX_VAR) + UR(MAGX_VAR))
    Gavg = 0.5_dp*(UL(GLMP_VAR) + UR(GLMP_VAR))

    !! Velocity.
    v = UI(MOMX_VAR:MOMZ_VAR)/UI(DENS_VAR)

    !! Powell term.
    mhdF(MOMX_VAR:MOMZ_VAR)    = Bavg * UI(MAGX_VAR:MAGZ_VAR)/mu0
    mhdF(ENER_VAR)             = Bavg * SUM(UI(MAGX_VAR:MAGZ_VAR)*v)/mu0
    mhdF(MAGX_VAR:MAGZ_VAR)    = Bavg * v

    !! GLM term.
    glmF((/ENER_VAR,GLMP_VAR/)) = Gavg*v(X_DIR) * (/UI(GLMP_VAR)/mu0,1.0_dp/)

    FL = mhdF + glmF

end function

pure function Non_Conservative_two_point_yFlux(UI,UL,UR) result(F)

    use equations_const_mod

    real(dp), intent(in)    :: UI(N_VARS)
    real(dp), intent(in)    :: UL(N_VARS)
    real(dp), intent(in)    :: UR(N_VARS)
    real(dp)                :: F(N_VARS)

    real(dp) :: ruI(N_VARS), ruL(N_VARS), ruR(N_VARS), rf(N_VARS) !! rotated

    ruI(DENS_VAR) = uI(DENS_VAR)

    ruI(MOMX_VAR) = uI(MOMY_VAR)
    ruI(MOMY_VAR) = uI(MOMZ_VAR)
    ruI(MOMZ_VAR) = uI(MOMX_VAR)

    ruI(ENER_VAR) = uI(ENER_VAR)
    
    ruI(MAGX_VAR) = uI(MAGY_VAR)
    ruI(MAGY_VAR) = uI(MAGZ_VAR)
    ruI(MAGZ_VAR) = uI(MAGX_VAR)

    ruI(GLMP_VAR) = uI(GLMP_VAR)

    !! -- !!

    ruL(DENS_VAR) = uL(DENS_VAR)

    ruL(MOMX_VAR) = uL(MOMY_VAR)
    ruL(MOMY_VAR) = uL(MOMZ_VAR)
    ruL(MOMZ_VAR) = uL(MOMX_VAR)

    ruL(ENER_VAR) = uL(ENER_VAR)
    
    ruL(MAGX_VAR) = uL(MAGY_VAR)
    ruL(MAGY_VAR) = uL(MAGZ_VAR)
    ruL(MAGZ_VAR) = uL(MAGX_VAR)

    ruL(GLMP_VAR) = uL(GLMP_VAR)

    !! -- !!

    ruR(DENS_VAR) = uR(DENS_VAR)

    ruR(MOMX_VAR) = uR(MOMY_VAR)
    ruR(MOMY_VAR) = uR(MOMZ_VAR)
    ruR(MOMZ_VAR) = uR(MOMX_VAR)

    ruR(ENER_VAR) = uR(ENER_VAR)

    ruR(MAGX_VAR) = uR(MAGY_VAR)
    ruR(MAGY_VAR) = uR(MAGZ_VAR)
    ruR(MAGZ_VAR) = uR(MAGX_VAR)

    ruR(GLMP_VAR) = uR(GLMP_VAR)

    !! -- !!

    rf = Non_Conservative_two_point_xFlux(rUI,rUL,rUR)

    !! -- !!

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

pure function Non_Conservative_two_point_zFlux(UI,UL,UR) result(F)

    use equations_const_mod

    real(dp), intent(in)    :: UI(N_VARS)
    real(dp), intent(in)    :: UL(N_VARS)
    real(dp), intent(in)    :: UR(N_VARS)
    real(dp)                :: F(N_VARS)

    real(dp) :: ruI(N_VARS), ruL(N_VARS), ruR(N_VARS), rf(N_VARS) !! rotated

    ruI(DENS_VAR) = uI(DENS_VAR)
                               
    ruI(MOMX_VAR) = uI(MOMZ_VAR)
    ruI(MOMY_VAR) = uI(MOMX_VAR)
    ruI(MOMZ_VAR) = uI(MOMY_VAR)
                               
    ruI(ENER_VAR) = uI(ENER_VAR)
                               
    ruI(MAGX_VAR) = uI(MAGZ_VAR)
    ruI(MAGY_VAR) = uI(MAGX_VAR)
    ruI(MAGZ_VAR) = uI(MAGY_VAR)
                               
    ruI(GLMP_VAR) = uI(GLMP_VAR)

    !! -- !!

    ruL(DENS_VAR) = uL(DENS_VAR)
                               
    ruL(MOMX_VAR) = uL(MOMZ_VAR)
    ruL(MOMY_VAR) = uL(MOMX_VAR)
    ruL(MOMZ_VAR) = uL(MOMY_VAR)
                               
    ruL(ENER_VAR) = uL(ENER_VAR)
                               
    ruL(MAGX_VAR) = uL(MAGZ_VAR)
    ruL(MAGY_VAR) = uL(MAGX_VAR)
    ruL(MAGZ_VAR) = uL(MAGY_VAR)
                               
    ruL(GLMP_VAR) = uL(GLMP_VAR)

    !! -- !!

    ruR(DENS_VAR) = uR(DENS_VAR)
                               
    ruR(MOMX_VAR) = uR(MOMZ_VAR)
    ruR(MOMY_VAR) = uR(MOMX_VAR)
    ruR(MOMZ_VAR) = uR(MOMY_VAR)
                               
    ruR(ENER_VAR) = uR(ENER_VAR)
                               
    ruR(MAGX_VAR) = uR(MAGZ_VAR)
    ruR(MAGY_VAR) = uR(MAGX_VAR)
    ruR(MAGZ_VAR) = uR(MAGY_VAR)
                               
    ruR(GLMP_VAR) = uR(GLMP_VAR)

    !! -- !!

    rf = Non_Conservative_two_point_xFlux(rUI,rUL,rUR)

    !! -- !!

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

end module
