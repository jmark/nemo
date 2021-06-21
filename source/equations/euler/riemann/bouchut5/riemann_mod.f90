module riemann_mod

use constants_mod

enum, bind(c)

    !! "NUNK" ... FLASH speak
    enumerator :: B5_NUNK_VARS_BEGIN = 0

    enumerator :: DENS_VAR
    enumerator :: MOMX_VAR
    enumerator :: MOMY_VAR
    enumerator :: MOMZ_VAR
    enumerator :: ENER_VAR

    enumerator :: MAGX_VAR
    enumerator :: MAGY_VAR
    enumerator :: MAGZ_VAR

    enumerator :: VELX_VAR
    enumerator :: VELY_VAR
    enumerator :: VELZ_VAR

    enumerator :: PRES_VAR
    enumerator :: EINT_VAR

    enumerator :: GAMC_VAR

    enumerator :: B5_NUNK_VARS_END

endenum

integer, parameter :: B5_NUNK_VARS = B5_NUNK_VARS_END - 1

enum, bind(c)

    enumerator :: B5_FLUX_VARS_BEGIN = 0

    enumerator :: DENS_FLUX
    enumerator :: XMOM_FLUX
    enumerator :: YMOM_FLUX
    enumerator :: ZMOM_FLUX
    enumerator :: ENER_FLUX

    enumerator :: MAGX_FLUX
    enumerator :: MAGY_FLUX
    enumerator :: MAGZ_FLUX

    enumerator :: VELX_FLUX
    enumerator :: VELY_FLUX
    enumerator :: VELZ_FLUX

    enumerator :: PRES_FLUX
    enumerator :: EINT_FLUX

    enumerator :: GAMC_FLUX

    enumerator :: B5_FLUX_VARS_END

endenum

integer, parameter :: B5_FLUX_VARS = B5_FLUX_VARS_END - 1

contains

!! 5-wave Riemann solver (Bouchut, 2008)
pure function xriemann(U_L,U_R) result(F)

    use equations_mod, only: pressure
    use setup_config_mod, only: kappa

    real(dp), intent(in) :: U_L(N_VARS),U_R(N_VARS)
    real(dp)             :: F(N_VARS)
    
    real(dp) :: speed, vint

    real(dp) :: lm,lp,lam,lap,lbm,lbp

    real(dp) :: cal, car, cbl, cbr, cbovrhol, cbovrhor, caovrhol, caovrhor, aql, aqr, aql2, aqr2
    real(dp) :: pixl, pixr, Bx2l, Bx2r, asl2, asr2, asl, asr, BxBtl,BxBtr, Xl, Xr
    real(dp) :: xsmr, xsml, B2l, B2r, alpha, Etotl, Etotr, ptotl, ptotr

    real(dp), DIMENSION(N_DIMS) :: pil, pir, ustar,  Bstarl, Bstarr, pinstar, ctabl, ctabr, ustarmid
    real(dp) :: pi2l, pi2r,pistar2, rhostarl, rhostarr, eintstarl, eintstarr
    real(dp) :: B2starl, B2starr, u2star, pistar,Etotstarl, Etotstarr
    real(dp) :: dBn

    real(dp) :: Um(B5_NUNK_VARS), Up(B5_NUNK_VARS), Flux(B5_FLUX_VARS)

    integer, parameter :: ivn  = VELX_VAR
    integer, parameter :: ibn  = MAGX_VAR
    integer, parameter :: idir = 1

    Up(DENS_VAR:ENER_VAR) = U_L
    Um(DENS_VAR:ENER_VAR) = U_R

    Up(MAGX_VAR:MAGZ_VAR) = 0.0_dp
    Um(MAGX_VAR:MAGZ_VAR) = 0.0_dp

    Up(VELX_VAR:VELZ_VAR) = Up(MOMX_VAR:MOMZ_VAR)/Up(DENS_VAR)
    Um(VELX_VAR:VELZ_VAR) = Um(MOMX_VAR:MOMZ_VAR)/Um(DENS_VAR)

    Up(PRES_VAR) = pressure(Up(DENS_VAR:ENER_VAR))
    Um(PRES_VAR) = pressure(Um(DENS_VAR:ENER_VAR))

    Up(EINT_VAR) = Up(PRES_VAR)/(kappa-1)/Up(DENS_VAR)
    Um(EINT_VAR) = Um(PRES_VAR)/(kappa-1)/Um(DENS_VAR)

    Up(GAMC_VAR) = kappa
    Um(GAMC_VAR) = kappa

    Flux  = 0.0_dp
    speed = 0.0_dp
    vint  = 0.0_dp

    alpha = 0.5_dp*(max(Up(GAMC_VAR),Um(GAMC_VAR))+1.0_dp) + 1e-4_dp

    Bx2l = Up(ibn)**2 
    Bx2r = Um(ibn)**2


    !! Total pressure and energy
    B2l   = Up(MAGX_VAR)**2 + Up(MAGY_VAR)**2 + Up(MAGZ_VAR)**2
    ptotl = Up(PRES_VAR) + 0.5_dp*B2l
    pixl  = ptotl - Bx2l    
    Etotl = Up(DENS_VAR)*(Up(EINT_VAR) + 0.5_dp*(Up(VELX_VAR)**2 + Up(VELY_VAR)**2 + Up(VELZ_VAR)**2)) + 0.5_dp*B2l

    B2r   = Um(MAGX_VAR)**2 + Um(MAGY_VAR)**2 + Um(MAGZ_VAR)**2
    ptotr = Um(PRES_VAR) + 0.5_dp*B2r
    pixr  = ptotr - Bx2r    
    Etotr = Um(DENS_VAR)*(Um(EINT_VAR) + 0.5_dp*(Um(VELX_VAR)**2 + Um(VELY_VAR)**2 + Um(VELZ_VAR)**2)) + 0.5_dp*B2r


    pil         = 0.0_dp
    pil(idir)   = ptotl
    pil         = pil - Up(ibn)*Up(MAGX_VAR:MAGZ_VAR)
    pi2l        = dot_product(pil,pil)

    pir         = 0.0_dp
    pir(idir)   = ptotr
    pir         = pir - Um(ibn)*Um(MAGX_VAR:MAGZ_VAR)
    pi2r        = dot_product(pir,pir)


    !! Sound soeed
    asl2        = Up(GAMC_VAR)*Up(PRES_VAR)/Up(DENS_VAR)
    asr2        = Um(GAMC_VAR)*Um(PRES_VAR)/Um(DENS_VAR)
    asl2        = max(0.0_dp,asl2)
    asr2        = max(0.0_dp,asr2)
    asl         = sqrt(asl2)
    asr         = sqrt(asr2)

    BxBtl       = abs(Up(ibn))*sqrt(max(0.0_dp,B2l-Bx2l))
    BxBtr       = abs(Um(ibn))*sqrt(max(0.0_dp,B2r-Bx2r))
    aql2        = abs(asl2 + (B2l + BxBtl)/Up(DENS_VAR))
    aqr2        = abs(asr2 + (B2r + BxBtr)/Um(DENS_VAR))
    aql         = sqrt(aql2)
    aqr         = sqrt(aqr2)

    Xl          = max(Up(ivn)-Um(ivn),0.0_dp) + max(pixr-pixl,0.0_dp)/(Up(DENS_VAR)*aql + Um(DENS_VAR)*aqr)
    Xr          = max(Up(ivn)-Um(ivn),0.0_dp) + max(pixl-pixr,0.0_dp)/(Um(DENS_VAR)*aqr + Up(DENS_VAR)*aql)

    xsml        = 1.0_dp - Xl/(aql+alpha*Xl)
    xsmr        = 1.0_dp - Xr/(aqr+alpha*Xr)

    caovrhol    = sqrt((Bx2l + BxBtl)/Up(DENS_VAR)/xsml)
    caovrhor    = sqrt((Bx2r + BxBtr)/Um(DENS_VAR)/xsmr)
    cbovrhol    = aql/sqrt(xsml) + alpha*Xl
    cbovrhor    = aqr/sqrt(xsmr) + alpha*Xr

    !! Avoid ill-posedness. 
    caovrhol    = max(caovrhol, 1.0e-12_dp*cbovrhol)
    caovrhor    = max(caovrhor, 1.0e-12_dp*cbovrhor)

    cal         = Up(DENS_VAR)*caovrhol
    car         = Um(DENS_VAR)*caovrhor
    cbl         = Up(DENS_VAR)*cbovrhol
    cbr         = Um(DENS_VAR)*cbovrhor

    !! Velocities
    ustarmid        = (cal*Up(VELX_VAR:VELZ_VAR) + car*Um(VELX_VAR:VELZ_VAR) + pil-pir) / (cal+car)
    ustar(idir)     = (cbl*Up(ivn) + cbr*Um(ivn) + pixl - pixr) / (cbl+cbr)
    ustarmid(idir)  = ustar(idir)

    lbm         = Up(ivn) - cbovrhol
    lbp         = Um(ivn) + cbovrhor

    lam         = Up(ivn) - cal/Up(DENS_VAR)
    lap         = Um(ivn) + car/Um(DENS_VAR)

    lm          = min(lam,lbm)
    lp          = max(lap,lbp)

    !! "Output"
    speed       = max(lp,-lm)
    vint        = ustar(idir)

    !! Computation of states and fluxes:

    if (lm .ge. 0.0_dp) then

        Flux(DENS_FLUX)           = Up(DENS_VAR)*Up(ivn)
        Flux(XMOM_FLUX:ZMOM_FLUX) = Flux(DENS_FLUX)*Up(VELX_VAR:VELZ_VAR) + pil
        Flux(MAGX_FLUX:MAGZ_FLUX) = Up(ivn) * Up(MAGX_VAR:MAGZ_VAR) - Up(ibn)*Up(VELX_VAR:VELZ_VAR)
        Flux(ENER_FLUX)           = Etotl*Up(ivn) + dot_product(pil,Up(VELX_VAR:VELZ_VAR))

    else if (ustar(idir) .gt. 0.0_dp) then

        if (lbm .gt. 0.0_dp) then

            ustar           = (cal*Up(VELX_VAR:VELZ_VAR) + car*Um(VELX_VAR:VELZ_VAR) +pil - pir) / (cal+car)
            ustar(idir)     = Up(ivn)
            pinstar         = (car*pil + cal*pir - cal*car*(Um(VELX_VAR:VELZ_VAR) - Up(VELX_VAR:VELZ_VAR))) / (cal+car)
            pinstar(idir)   = pixl

            rhostarl        = Up(DENS_VAR)

            Bstarl          = Up(MAGX_VAR:MAGZ_VAR) + rhostarl*Up(ibn)*(pil-pinstar)/(cal**2)
            Bstarl(idir)    = Up(ibn)
            u2star          = dot_product(ustar,ustar)
            B2starl         = dot_product(Bstarl,Bstarl)
            pistar2         = dot_product(pinstar,pinstar)
            eintstarl       = Up(EINT_VAR) + 0.5_dp*(B2l/Up(DENS_VAR) - B2starl/rhostarl) - (pi2l - pistar2)/(2.0_dp*cal**2)
            Etotstarl       = rhostarl*(0.5_dp*u2star + eintstarl) + 0.5_dp*B2starl

            Flux(DENS_FLUX)             = rhostarl*ustar(idir)
            Flux(XMOM_FLUX:ZMOM_FLUX)   = rhostarl*ustar(idir)*ustar + pinstar
            Flux(MAGX_FLUX:MAGZ_FLUX)   = ustar(idir)*Bstarl - Up(ibn)*ustar
            Flux(ENER_FLUX)             = Etotstarl*ustar(idir) + dot_product(pinstar(1:3),ustar(1:3))

        else if (lam .lt. 0.0_dp) then

            ctabl       = (/cal, cal, cal/)
            ctabl(idir) = cbl
            ctabr       = (/car, car, car/)
            ctabr(idir) = cbr

            ustar           = (ctabl*Up(VELX_VAR:VELZ_VAR) + ctabr*Um(VELX_VAR:VELZ_VAR) + pil - pir) / (ctabl+ctabr)
            !ustar(idir)    = (cbl*Up(ivn,i-1)+cbr*Um(ivn,i)+pixl-pixr) / (cbl+cbr)
            pinstar         = (ctabr*pil+ctabl*pir-ctabl*ctabr*(Um(VELX_VAR:VELZ_VAR) - Up(VELX_VAR:VELZ_VAR))) / (ctabl+ctabr)
            !pinstar(idir)  = (cbr*pixl+cbl*pixr-cbl*cbr*(Um(ivn,i)-Up(ivn,i-1))) / (cnl+cnr)

            rhostarl        = 1.0_dp/(1.0_dp/Up(DENS_VAR) + (pil(idir) - pinstar(idir))/(cbl**2))
            Bstarl          = rhostarl*(Up(MAGX_VAR:MAGZ_VAR)/Up(DENS_VAR) + Up(ibn)*(pil-pinstar)/(ctabl**2))
            u2star          = dot_product(ustar,ustar)
            B2starl         = dot_product(Bstarl,Bstarl)
            pi2l            = dot_product(pil/ctabl,pil/ctabl)
            pistar2         = dot_product(pinstar/ctabl,pinstar/ctabl)
            eintstarl       = Up(EINT_VAR) + 0.5_dp*(B2l/Up(DENS_VAR) - B2starl/rhostarl) - 0.5_dp*(pi2l-pistar2)
            Etotstarl       = rhostarl*(0.5_dp*u2star + eintstarl) + 0.5_dp*B2starl

            Flux(DENS_FLUX)             = rhostarl*ustar(idir)
            Flux(XMOM_FLUX:ZMOM_FLUX)   = rhostarl*ustar(idir)*ustar + pinstar
            Flux(MAGX_FLUX:MAGZ_FLUX)   = ustar(idir)*Bstarl - Up(ibn)*ustar
            Flux(ENER_FLUX)             = Etotstarl*ustar(idir)+dot_product(pinstar(1:3),ustar(1:3))

        else !lam > 0, lbm < 0

            ustar           = Up(VELX_VAR:VELZ_VAR)
            ustar(idir)     = vint
            pinstar         = pil
            pinstar(idir)   = (cbr*pixl+cbl*pixr-cbl*cbr*(Um(ivn)-Up(ivn))) / (cbl+cbr)

            rhostarl        = 1.0_dp/(1.0_dp/Up(DENS_VAR) + (pil(idir)-pinstar(idir))/(cbl**2))
            Bstarl          = rhostarl*(Up(MAGX_VAR:MAGZ_VAR)/Up(DENS_VAR))
            Bstarl(idir)    = Up(ibn)
            u2star          = dot_product(ustar,ustar)
            B2starl         = dot_product(Bstarl,Bstarl)
            !pi2l           = dot_product(pil,pil)
            !pistar2        = dot_product(pinstar,pinstar)
            eintstarl       = Up(EINT_VAR) + 0.5_dp*(B2l/Up(DENS_VAR) - B2starl/rhostarl)-(pixl**2-pinstar(idir)**2)/(2.0_dp*cbl**2)
            Etotstarl       = rhostarl*(u2star/2. + eintstarl) + B2starl/2

            Flux(DENS_FLUX)             = rhostarl*ustar(idir)
            Flux(XMOM_FLUX:ZMOM_FLUX)   = rhostarl*ustar(idir)*ustar + pinstar
            Flux(MAGX_FLUX:MAGZ_FLUX)   = ustar(idir)*Bstarl - Up(ibn)*ustar
            Flux(ENER_FLUX)             = Etotstarl*ustar(idir) + dot_product(pinstar(1:3),ustar(1:3))

        endif ! Done with left wave fan

    else if (lp.gt.0.0_dp) then

        if (lbp .lt. 0.0_dp ) then !lap>0

            ustar           = (cal*Up(VELX_VAR:VELZ_VAR) + car*Um(VELX_VAR:VELZ_VAR) + pil-pir) / (cal+car)
            ustar(idir)     = Um(ivn)
            pinstar         = (car*pil + cal*pir - cal*car*(Um(VELX_VAR:VELZ_VAR) - Up(VELX_VAR:VELZ_VAR))) / (cal+car)
            pinstar(idir)   = pixr

            rhostarr        = Um(DENS_VAR)

            Bstarr          = Um(MAGX_VAR:MAGZ_VAR) + rhostarr*Um(ibn)*(pir-pinstar)/(car**2)
            Bstarr(idir)    = Um(ibn)
            u2star          = dot_product(ustar,ustar)
            B2starr         = dot_product(Bstarr,Bstarr)
            pistar2         = dot_product(pinstar,pinstar)
            eintstarr       = Um(EINT_VAR) + (B2r/Um(DENS_VAR) - B2starr/rhostarr)/2 - (pi2r-pistar2)/(2*car**2)
            Etotstarr       = rhostarr*(u2star/2 + eintstarr) + B2starr/2

            Flux(DENS_FLUX)             = rhostarr*ustar(idir)
            Flux(XMOM_FLUX:ZMOM_FLUX)   = rhostarr*ustar(idir)*ustar + pinstar
            Flux(MAGX_FLUX:MAGZ_FLUX)   = ustar(idir)*Bstarr - Um(ibn)*ustar
            Flux(ENER_FLUX)             = Etotstarr*ustar(idir) + dot_product(pinstar(1:3),ustar(1:3))

        else if (lap .gt. 0) then !lbp>0, lap>0

            ctabl           = (/cal, cal, cal/)
            ctabl(idir)     = cbl
            ctabr           = (/car, car,car/)
            ctabr(idir)     = cbr

            ustar           = (ctabl*Up(VELX_VAR:VELZ_VAR) + ctabr*Um(VELX_VAR:VELZ_VAR) + pil-pir) / (ctabl+ctabr)
            !ustar(idir)    = (cbl*Up(ivn,i-1)+cbr*Um(ivn,i)+pixl-pixr) / (cbl+cbr)
            pinstar         = (ctabr*pil+ctabl*pir-ctabl*ctabr*(Um(VELX_VAR:VELZ_VAR) - Up(VELX_VAR:VELZ_VAR))) / (ctabl+ctabr)
            !pinstar(idir)  = (cbr*pixl+cbl*pixr-cbl*cbr*(Um(ivn,i)-Up(ivn,i-1))) / (cnl+cnr)

            rhostarr        = 1.0_dp/(1.0_dp/Um(DENS_VAR) + (pir(idir) - pinstar(idir))/(cbr**2))
            Bstarr          = rhostarr*(Um(MAGX_VAR:MAGZ_VAR)/Um(DENS_VAR) + Um(ibn)*(pir-pinstar)/(ctabr**2))
            u2star          = dot_product(ustar,ustar)
            B2starr         = dot_product(Bstarr,Bstarr)
            pi2r            = dot_product(pir/ctabr,pir/ctabr)
            pistar2         = dot_product(pinstar/ctabr,pinstar/ctabr)
            eintstarr       = Um(EINT_VAR) + (B2r/Um(DENS_VAR) - B2starr/rhostarr)/2 - (pi2r-pistar2)/2
            Etotstarr       = rhostarr*(u2star/2 + eintstarr) + B2starr/2

            Flux(DENS_FLUX)             = rhostarr*ustar(idir)
            Flux(XMOM_FLUX:ZMOM_FLUX)   = rhostarr*ustar(idir)*ustar + pinstar
            Flux(MAGX_FLUX:MAGZ_FLUX)   = ustar(idir)*Bstarr - Um(ibn)*ustar
            Flux(ENER_FLUX)             = Etotstarr*ustar(idir) + dot_product(pinstar(1:3),ustar(1:3))

        else !lbp>0, lap<=0

            ustar           = Um(VELX_VAR:VELZ_VAR)
            ustar(idir)     = vint
            pinstar         = pir
            pinstar(idir)   = (cbr*pixl + cbl*pixr - cbl*cbr*(Um(ivn)-Up(ivn))) / (cbl+cbr)

            rhostarr        = 1.0_dp/(1.0_dp/Um(DENS_VAR) + (pir(idir) - pinstar(idir))/(cbr**2))
            Bstarr          = rhostarr*(Um(MAGX_VAR:MAGZ_VAR)/Um(DENS_VAR))
            Bstarr(idir)    = Um(ibn)
            u2star          = dot_product(ustar,ustar)
            B2starr         = dot_product(Bstarr,Bstarr)
            !pi2l=dot_product(pil,pil)
            !pistar2=dot_product(pinstar,pinstar)
            eintstarr       = Um(EINT_VAR) + (B2r/Um(DENS_VAR) - B2starr/rhostarr)/2 - (pixr**2-pinstar(idir)**2)/(2*cbr**2)
            Etotstarr       = rhostarr*(u2star/2 + eintstarr) + B2starr/2

            Flux(DENS_FLUX)             = rhostarr*ustar(idir)
            Flux(XMOM_FLUX:ZMOM_FLUX)   = rhostarr*ustar(idir)*ustar + pinstar
            Flux(MAGX_FLUX:MAGZ_FLUX)   = ustar(idir)*Bstarr - Um(ibn)*ustar
            Flux(ENER_FLUX)             = Etotstarr*ustar(idir) + dot_product(pinstar(1:3),ustar(1:3))

        endif !! End of right inner states 

    else ! Leftwards supermagnetosonic

        Flux(DENS_FLUX)             = Um(DENS_VAR)*Um(ivn)
        Flux(XMOM_FLUX:ZMOM_FLUX)   = Flux(DENS_FLUX)*Um(VELX_VAR:VELZ_VAR) + pir
        Flux(MAGX_FLUX:MAGZ_FLUX)   = Um(ivn)*Um(MAGX_VAR:MAGZ_VAR) - Um(ibn)*Um(VELX_VAR:VELZ_VAR)
        Flux(ENER_FLUX)             = Etotr*Um(ivn) + dot_product(pir,Um(VELX_VAR:VELZ_VAR))

    endif

    !! Nonconservative bit added to source
    !! dBn = Um(ibn)-Up(ibn)
    !! if (vint .le. 0.0_dp) then
    !!     Src(MAGX_VAR:MAGZ_VAR) = Src(MAGX_VAR:MAGZ_VAR) - ustarmid*dBn/dx
    !! else
    !!     Src(MAGX_VAR:MAGZ_VAR) = Src(MAGX_VAR:MAGZ_VAR) - ustarmid*dBn/dx
    !! endif

    F(DENS_VAR)             = Flux(DENS_FLUX)
    F(MOMX_VAR:MOMZ_VAR)    = Flux(XMOM_FLUX:ZMOM_FLUX)
    F(ENER_VAR)             = Flux(ENER_FLUX)

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
