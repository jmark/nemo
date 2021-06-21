module setup_manufacsol_mod

    use constants_mod
    use setup_config_mod
    use equations_const_mod

contains

pure function init_manufacsol(x,y,z) result(u)

    real(dp), intent(in)    :: x,y,z
    real(dp)                :: u(N_VARS)

    u(DENS_VAR) = d0 + d1*Sin((2*Pi*df1*x)/L) + d2*Cos((2*Pi*df2*y)/L) + d3*Sin((2*Pi*df3*z)/L) + d4*Cos((2*Pi*df4*x*y*z)/L**3);
    u(MOMX_VAR) = u0 + u1*Cos((2*Pi*uf1*x)/L) + u2*Sin((2*Pi*uf2*y)/L) + u3*Cos((2*Pi*uf3*z)/L) + u4*Sin((2*Pi*uf4*x*y*z)/L**3);
    u(MOMY_VAR) = v0 + v1*Sin((2*Pi*vf1*x)/L) + v2*Cos((2*Pi*vf2*y)/L) + v3*Sin((2*Pi*vf3*z)/L) + v4*Cos((2*Pi*vf4*x*y*z)/L**3);
    ! u(MOMZ_VAR) = w0 + w1*Cos((2*Pi*wf1*x)/L) + w2*Sin((2*Pi*wf2*y)/L) + w3*Cos((2*Pi*wf3*z)/L) + w4*Sin((2*Pi*wf4*x*y*z)/L**3);
    u(ENER_VAR) = p0 + p1*Sin((2*Pi*pf1*x)/L) + p2*Cos((2*Pi*pf2*y)/L) + p3*Sin((2*Pi*pf3*z)/L) + p4*Sin((2*Pi*pf4*x*y*z)/L**3);

end function

pure function calc_manufacsol(x,y,z) result(r)

    real(dp), intent(in)    :: x,y,z
    real(dp)                :: r(N_VARS)

    r = 0.0

    !! F function

    r(DENS_VAR) = r(DENS_VAR) + &
        ((2*Pi*u4*uf4*y*z*Cos((2*Pi*uf4*x*y*z)/L**3))/L**3 - (2*Pi*u1*uf1*Sin((2*Pi*uf1*x)/L))/L)*(d0 + d2*Cos((2*df2*Pi*y)/L) + d4*Cos((2*df4*Pi*x*y*z)/L**3) + d1*Sin((2*df1*Pi*x)/L) + d3*Sin((2*df3*Pi*z)/L)) + &
        &  ((2*d1*df1*Pi*Cos((2*df1*Pi*x)/L))/L - (2*d4*df4*Pi*y*z*Sin((2*df4*Pi*x*y*z)/L**3))/L**3)*(u0 + u1*Cos((2*Pi*uf1*x)/L) + u3*Cos((2*Pi*uf3*z)/L) + u2*Sin((2*Pi*uf2*y)/L) + u4*Sin((2*Pi*uf4*x*y*z)/L**3));

    r(MOMX_VAR) = r(MOMX_VAR) + &
        & (2*p1*pf1*Pi*Cos((2*pf1*Pi*x)/L))/L + (2*p4*pf4*Pi*y*z*Cos((2*pf4*Pi*x*y*z)/L**3))/L**3 + &
        & 2*((2*Pi*u4*uf4*y*z*Cos((2*Pi*uf4*x*y*z)/L**3))/L**3 - (2*Pi*u1*uf1*Sin((2*Pi*uf1*x)/L))/L)*(d0 + d2*Cos((2*df2*Pi*y)/L) + d4*Cos((2*df4*Pi*x*y*z)/L**3) + d1*Sin((2*df1*Pi*x)/L) + d3*Sin((2*df3*Pi*z)/L))* &
        &  (u0 + u1*Cos((2*Pi*uf1*x)/L) + u3*Cos((2*Pi*uf3*z)/L) + u2*Sin((2*Pi*uf2*y)/L) + u4*Sin((2*Pi*uf4*x*y*z)/L**3)) + &
        & ((2*d1*df1*Pi*Cos((2*df1*Pi*x)/L))/L - (2*d4*df4*Pi*y*z*Sin((2*df4*Pi*x*y*z)/L**3))/L**3)*(u0 + u1*Cos((2*Pi*uf1*x)/L) + u3*Cos((2*Pi*uf3*z)/L) + u2*Sin((2*Pi*uf2*y)/L) + u4*Sin((2*Pi*uf4*x*y*z)/L**3))**2;

    r(MOMY_VAR) = r(MOMY_VAR) + &
        & ((2*Pi*u4*uf4*y*z*Cos((2*Pi*uf4*x*y*z)/L**3))/L**3 - (2*Pi*u1*uf1*Sin((2*Pi*uf1*x)/L))/L)*(d0 + d2*Cos((2*df2*Pi*y)/L) + d4*Cos((2*df4*Pi*x*y*z)/L**3) + d1*Sin((2*df1*Pi*x)/L) + d3*Sin((2*df3*Pi*z)/L))* &
        &  (v0 + v2*Cos((2*Pi*vf2*y)/L) + v4*Cos((2*Pi*vf4*x*y*z)/L**3) + v1*Sin((2*Pi*vf1*x)/L) + v3*Sin((2*Pi*vf3*z)/L)) + &
        & (v0 + v2*Cos((2*Pi*vf2*y)/L) + v4*Cos((2*Pi*vf4*x*y*z)/L**3) + v1*Sin((2*Pi*vf1*x)/L) + v3*Sin((2*Pi*vf3*z)/L))*((2*d1*df1*Pi*Cos((2*df1*Pi*x)/L))/L - (2*d4*df4*Pi*y*z*Sin((2*df4*Pi*x*y*z)/L**3))/L**3)* &
        &  (u0 + u1*Cos((2*Pi*uf1*x)/L) + u3*Cos((2*Pi*uf3*z)/L) + u2*Sin((2*Pi*uf2*y)/L) + u4*Sin((2*Pi*uf4*x*y*z)/L**3)) + &
        & (d0 + d2*Cos((2*df2*Pi*y)/L) + d4*Cos((2*df4*Pi*x*y*z)/L**3) + d1*Sin((2*df1*Pi*x)/L) + d3*Sin((2*df3*Pi*z)/L))* &
        &  (u0 + u1*Cos((2*Pi*uf1*x)/L) + u3*Cos((2*Pi*uf3*z)/L) + u2*Sin((2*Pi*uf2*y)/L) + u4*Sin((2*Pi*uf4*x*y*z)/L**3))*((2*Pi*v1*vf1*Cos((2*Pi*vf1*x)/L))/L - (2*Pi*v4*vf4*y*z*Sin((2*Pi*vf4*x*y*z)/L**3))/L**3);

    ! r(MOMZ_VAR) = r(MOMZ_VAR) + &
    !     & ((2*Pi*w4*wf4*y*z*Cos((2*Pi*wf4*x*y*z)/L**3))/L**3 - (2*Pi*w1*wf1*Sin((2*Pi*wf1*x)/L))/L)*(d0 + d2*Cos((2*df2*Pi*y)/L) + d4*Cos((2*df4*Pi*x*y*z)/L**3) + d1*Sin((2*df1*Pi*x)/L) + d3*Sin((2*df3*Pi*z)/L))* &
    !     &  (u0 + u1*Cos((2*Pi*uf1*x)/L) + u3*Cos((2*Pi*uf3*z)/L) + u2*Sin((2*Pi*uf2*y)/L) + u4*Sin((2*Pi*uf4*x*y*z)/L**3)) + &
    !     & ((2*Pi*u4*uf4*y*z*Cos((2*Pi*uf4*x*y*z)/L**3))/L**3 - (2*Pi*u1*uf1*Sin((2*Pi*uf1*x)/L))/L)*(d0 + d2*Cos((2*df2*Pi*y)/L) + d4*Cos((2*df4*Pi*x*y*z)/L**3) + d1*Sin((2*df1*Pi*x)/L) + d3*Sin((2*df3*Pi*z)/L))* &
    !     &  (w0 + w1*Cos((2*Pi*wf1*x)/L) + w3*Cos((2*Pi*wf3*z)/L) + w2*Sin((2*Pi*wf2*y)/L) + w4*Sin((2*Pi*wf4*x*y*z)/L**3)) + &
    !     & ((2*d1*df1*Pi*Cos((2*df1*Pi*x)/L))/L - (2*d4*df4*Pi*y*z*Sin((2*df4*Pi*x*y*z)/L**3))/L**3)*(u0 + u1*Cos((2*Pi*uf1*x)/L) + u3*Cos((2*Pi*uf3*z)/L) + u2*Sin((2*Pi*uf2*y)/L) + u4*Sin((2*Pi*uf4*x*y*z)/L**3))* &
    !     &  (w0 + w1*Cos((2*Pi*wf1*x)/L) + w3*Cos((2*Pi*wf3*z)/L) + w2*Sin((2*Pi*wf2*y)/L) + w4*Sin((2*Pi*wf4*x*y*z)/L**3));
         
    r(ENER_VAR) = r(ENER_VAR) + &
        & ((2*Pi*u4*uf4*y*z*Cos((2*Pi*uf4*x*y*z)/L**3))/L**3 - (2*Pi*u1*uf1*Sin((2*Pi*uf1*x)/L))/L)* &
        &    (p0 + p2*Cos((2*pf2*Pi*y)/L) + p1*Sin((2*pf1*Pi*x)/L) + p3*Sin((2*pf3*Pi*z)/L) + p4*Sin((2*pf4*Pi*x*y*z)/L**3) +  &
        &      (p0 + p2*Cos((2*pf2*Pi*y)/L) + p1*Sin((2*pf1*Pi*x)/L) + p3*Sin((2*pf3*Pi*z)/L) + p4*Sin((2*pf4*Pi*x*y*z)/L**3))/(-1 + kappa) +  &
        &      ((d0 + d2*Cos((2*df2*Pi*y)/L) + d4*Cos((2*df4*Pi*x*y*z)/L**3) + d1*Sin((2*df1*Pi*x)/L) + d3*Sin((2*df3*Pi*z)/L))*&
        &         ((v0 + v2*Cos((2*Pi*vf2*y)/L) + v4*Cos((2*Pi*vf4*x*y*z)/L**3) + v1*Sin((2*Pi*vf1*x)/L) + v3*Sin((2*Pi*vf3*z)/L))**2 +  &
        &           (u0 + u1*Cos((2*Pi*uf1*x)/L) + u3*Cos((2*Pi*uf3*z)/L) + u2*Sin((2*Pi*uf2*y)/L) + u4*Sin((2*Pi*uf4*x*y*z)/L**3))**2 +  &
        &           (w0 + w1*Cos((2*Pi*wf1*x)/L) + w3*Cos((2*Pi*wf3*z)/L) + w2*Sin((2*Pi*wf2*y)/L) + w4*Sin((2*Pi*wf4*x*y*z)/L**3))**2))/2.) +  &
        &   (u0 + u1*Cos((2*Pi*uf1*x)/L) + u3*Cos((2*Pi*uf3*z)/L) + u2*Sin((2*Pi*uf2*y)/L) + u4*Sin((2*Pi*uf4*x*y*z)/L**3))*&
        &    ((2*p1*pf1*Pi*Cos((2*pf1*Pi*x)/L))/L + (2*p4*pf4*Pi*y*z*Cos((2*pf4*Pi*x*y*z)/L**3))/L**3 + ((2*p1*pf1*Pi*Cos((2*pf1*Pi*x)/L))/L + (2*p4*pf4*Pi*y*z*Cos((2*pf4*Pi*x*y*z)/L**3))/L**3)/(-1 + kappa) +  &
        &      ((d0 + d2*Cos((2*df2*Pi*y)/L) + d4*Cos((2*df4*Pi*x*y*z)/L**3) + d1*Sin((2*df1*Pi*x)/L) + d3*Sin((2*df3*Pi*z)/L))*&
        &         (2*((2*Pi*u4*uf4*y*z*Cos((2*Pi*uf4*x*y*z)/L**3))/L**3 - (2*Pi*u1*uf1*Sin((2*Pi*uf1*x)/L))/L)*&
        &            (u0 + u1*Cos((2*Pi*uf1*x)/L) + u3*Cos((2*Pi*uf3*z)/L) + u2*Sin((2*Pi*uf2*y)/L) + u4*Sin((2*Pi*uf4*x*y*z)/L**3)) +  &
        &           2*(v0 + v2*Cos((2*Pi*vf2*y)/L) + v4*Cos((2*Pi*vf4*x*y*z)/L**3) + v1*Sin((2*Pi*vf1*x)/L) + v3*Sin((2*Pi*vf3*z)/L))*&
        &            ((2*Pi*v1*vf1*Cos((2*Pi*vf1*x)/L))/L - (2*Pi*v4*vf4*y*z*Sin((2*Pi*vf4*x*y*z)/L**3))/L**3) +  &
        &           2*((2*Pi*w4*wf4*y*z*Cos((2*Pi*wf4*x*y*z)/L**3))/L**3 - (2*Pi*w1*wf1*Sin((2*Pi*wf1*x)/L))/L)*&
        &            (w0 + w1*Cos((2*Pi*wf1*x)/L) + w3*Cos((2*Pi*wf3*z)/L) + w2*Sin((2*Pi*wf2*y)/L) + w4*Sin((2*Pi*wf4*x*y*z)/L**3))))/2. +  &
        &      (((2*d1*df1*Pi*Cos((2*df1*Pi*x)/L))/L - (2*d4*df4*Pi*y*z*Sin((2*df4*Pi*x*y*z)/L**3))/L**3)*&
        &         ((v0 + v2*Cos((2*Pi*vf2*y)/L) + v4*Cos((2*Pi*vf4*x*y*z)/L**3) + v1*Sin((2*Pi*vf1*x)/L) + v3*Sin((2*Pi*vf3*z)/L))**2 +  &
        &           (u0 + u1*Cos((2*Pi*uf1*x)/L) + u3*Cos((2*Pi*uf3*z)/L) + u2*Sin((2*Pi*uf2*y)/L) + u4*Sin((2*Pi*uf4*x*y*z)/L**3))**2 +  &
        &           (w0 + w1*Cos((2*Pi*wf1*x)/L) + w3*Cos((2*Pi*wf3*z)/L) + w2*Sin((2*Pi*wf2*y)/L) + w4*Sin((2*Pi*wf4*x*y*z)/L**3))**2))/2.);
         
    !! G function

    r(DENS_VAR) = r(DENS_VAR) + &
         & (v0 + v2*Cos((2*Pi*vf2*y)/L) + v4*Cos((2*Pi*vf4*x*y*z)/L**3) + v1*Sin((2*Pi*vf1*x)/L) + v3*Sin((2*Pi*vf3*z)/L))*((-2*d2*df2*Pi*Sin((2*df2*Pi*y)/L))/L - (2*d4*df4*Pi*x*z*Sin((2*df4*Pi*x*y*z)/L**3))/L**3) + &
         & (d0 + d2*Cos((2*df2*Pi*y)/L) + d4*Cos((2*df4*Pi*x*y*z)/L**3) + d1*Sin((2*df1*Pi*x)/L) + d3*Sin((2*df3*Pi*z)/L))*((-2*Pi*v2*vf2*Sin((2*Pi*vf2*y)/L))/L - (2*Pi*v4*vf4*x*z*Sin((2*Pi*vf4*x*y*z)/L**3))/L**3);

    r(MOMX_VAR) = r(MOMX_VAR) + &
        & ((2*Pi*u2*uf2*Cos((2*Pi*uf2*y)/L))/L + (2*Pi*u4*uf4*x*z*Cos((2*Pi*uf4*x*y*z)/L**3))/L**3)*(d0 + d2*Cos((2*df2*Pi*y)/L) + d4*Cos((2*df4*Pi*x*y*z)/L**3) + d1*Sin((2*df1*Pi*x)/L) + d3*Sin((2*df3*Pi*z)/L))*&
        &   (v0 + v2*Cos((2*Pi*vf2*y)/L) + v4*Cos((2*Pi*vf4*x*y*z)/L**3) + v1*Sin((2*Pi*vf1*x)/L) + v3*Sin((2*Pi*vf3*z)/L)) + &
        &  (v0 + v2*Cos((2*Pi*vf2*y)/L) + v4*Cos((2*Pi*vf4*x*y*z)/L**3) + v1*Sin((2*Pi*vf1*x)/L) + v3*Sin((2*Pi*vf3*z)/L))*((-2*d2*df2*Pi*Sin((2*df2*Pi*y)/L))/L - (2*d4*df4*Pi*x*z*Sin((2*df4*Pi*x*y*z)/L**3))/L**3)*&
        &   (u0 + u1*Cos((2*Pi*uf1*x)/L) + u3*Cos((2*Pi*uf3*z)/L) + u2*Sin((2*Pi*uf2*y)/L) + u4*Sin((2*Pi*uf4*x*y*z)/L**3)) + &
        &  (d0 + d2*Cos((2*df2*Pi*y)/L) + d4*Cos((2*df4*Pi*x*y*z)/L**3) + d1*Sin((2*df1*Pi*x)/L) + d3*Sin((2*df3*Pi*z)/L))*&
        &   (u0 + u1*Cos((2*Pi*uf1*x)/L) + u3*Cos((2*Pi*uf3*z)/L) + u2*Sin((2*Pi*uf2*y)/L) + u4*Sin((2*Pi*uf4*x*y*z)/L**3))*((-2*Pi*v2*vf2*Sin((2*Pi*vf2*y)/L))/L - (2*Pi*v4*vf4*x*z*Sin((2*Pi*vf4*x*y*z)/L**3))/L**3);
         
    r(MOMY_VAR) = r(MOMY_VAR) + &
        & (2*p4*pf4*Pi*x*z*Cos((2*pf4*Pi*x*y*z)/L**3))/L**3 - (2*p2*pf2*Pi*Sin((2*pf2*Pi*y)/L))/L + &
        &  (v0 + v2*Cos((2*Pi*vf2*y)/L) + v4*Cos((2*Pi*vf4*x*y*z)/L**3) + v1*Sin((2*Pi*vf1*x)/L) + v3*Sin((2*Pi*vf3*z)/L))**2*&
        &   ((-2*d2*df2*Pi*Sin((2*df2*Pi*y)/L))/L - (2*d4*df4*Pi*x*z*Sin((2*df4*Pi*x*y*z)/L**3))/L**3) + &
        &  2*(d0 + d2*Cos((2*df2*Pi*y)/L) + d4*Cos((2*df4*Pi*x*y*z)/L**3) + d1*Sin((2*df1*Pi*x)/L) + d3*Sin((2*df3*Pi*z)/L))*&
        &   (v0 + v2*Cos((2*Pi*vf2*y)/L) + v4*Cos((2*Pi*vf4*x*y*z)/L**3) + v1*Sin((2*Pi*vf1*x)/L) + v3*Sin((2*Pi*vf3*z)/L))*((-2*Pi*v2*vf2*Sin((2*Pi*vf2*y)/L))/L - (2*Pi*v4*vf4*x*z*Sin((2*Pi*vf4*x*y*z)/L**3))/L**3);
         
    ! r(MOMZ_VAR) = r(MOMZ_VAR) + &
    !     & ((2*Pi*w2*wf2*Cos((2*Pi*wf2*y)/L))/L + (2*Pi*w4*wf4*x*z*Cos((2*Pi*wf4*x*y*z)/L**3))/L**3)*(d0 + d2*Cos((2*df2*Pi*y)/L) + d4*Cos((2*df4*Pi*x*y*z)/L**3) + d1*Sin((2*df1*Pi*x)/L) + d3*Sin((2*df3*Pi*z)/L))*&
    !     &  (v0 + v2*Cos((2*Pi*vf2*y)/L) + v4*Cos((2*Pi*vf4*x*y*z)/L**3) + v1*Sin((2*Pi*vf1*x)/L) + v3*Sin((2*Pi*vf3*z)/L)) + &
    !     & (v0 + v2*Cos((2*Pi*vf2*y)/L) + v4*Cos((2*Pi*vf4*x*y*z)/L**3) + v1*Sin((2*Pi*vf1*x)/L) + v3*Sin((2*Pi*vf3*z)/L))*((-2*d2*df2*Pi*Sin((2*df2*Pi*y)/L))/L - (2*d4*df4*Pi*x*z*Sin((2*df4*Pi*x*y*z)/L**3))/L**3)*&
    !     &  (w0 + w1*Cos((2*Pi*wf1*x)/L) + w3*Cos((2*Pi*wf3*z)/L) + w2*Sin((2*Pi*wf2*y)/L) + w4*Sin((2*Pi*wf4*x*y*z)/L**3)) + &
    !     & (d0 + d2*Cos((2*df2*Pi*y)/L) + d4*Cos((2*df4*Pi*x*y*z)/L**3) + d1*Sin((2*df1*Pi*x)/L) + d3*Sin((2*df3*Pi*z)/L))*((-2*Pi*v2*vf2*Sin((2*Pi*vf2*y)/L))/L - (2*Pi*v4*vf4*x*z*Sin((2*Pi*vf4*x*y*z)/L**3))/L**3)*&
    !     &  (w0 + w1*Cos((2*Pi*wf1*x)/L) + w3*Cos((2*Pi*wf3*z)/L) + w2*Sin((2*Pi*wf2*y)/L) + w4*Sin((2*Pi*wf4*x*y*z)/L**3));

    r(ENER_VAR) = r(ENER_VAR) + &
        & ((-2*Pi*v2*vf2*Sin((2*Pi*vf2*y)/L))/L - (2*Pi*v4*vf4*x*z*Sin((2*Pi*vf4*x*y*z)/L**3))/L**3)*&
        &  (p0 + p2*Cos((2*pf2*Pi*y)/L) + p1*Sin((2*pf1*Pi*x)/L) + p3*Sin((2*pf3*Pi*z)/L) + p4*Sin((2*pf4*Pi*x*y*z)/L**3) + &
        &    (p0 + p2*Cos((2*pf2*Pi*y)/L) + p1*Sin((2*pf1*Pi*x)/L) + p3*Sin((2*pf3*Pi*z)/L) + p4*Sin((2*pf4*Pi*x*y*z)/L**3))/(-1 + kappa) + &
        &    ((d0 + d2*Cos((2*df2*Pi*y)/L) + d4*Cos((2*df4*Pi*x*y*z)/L**3) + d1*Sin((2*df1*Pi*x)/L) + d3*Sin((2*df3*Pi*z)/L))*&
        &       ((v0 + v2*Cos((2*Pi*vf2*y)/L) + v4*Cos((2*Pi*vf4*x*y*z)/L**3) + v1*Sin((2*Pi*vf1*x)/L) + v3*Sin((2*Pi*vf3*z)/L))**2 + &
        &         (u0 + u1*Cos((2*Pi*uf1*x)/L) + u3*Cos((2*Pi*uf3*z)/L) + u2*Sin((2*Pi*uf2*y)/L) + u4*Sin((2*Pi*uf4*x*y*z)/L**3))**2 + &
        &         (w0 + w1*Cos((2*Pi*wf1*x)/L) + w3*Cos((2*Pi*wf3*z)/L) + w2*Sin((2*Pi*wf2*y)/L) + w4*Sin((2*Pi*wf4*x*y*z)/L**3))**2))/2.) + &
        & (v0 + v2*Cos((2*Pi*vf2*y)/L) + v4*Cos((2*Pi*vf4*x*y*z)/L**3) + v1*Sin((2*Pi*vf1*x)/L) + v3*Sin((2*Pi*vf3*z)/L))*&
        &  ((2*p4*pf4*Pi*x*z*Cos((2*pf4*Pi*x*y*z)/L**3))/L**3 - (2*p2*pf2*Pi*Sin((2*pf2*Pi*y)/L))/L + ((2*p4*pf4*Pi*x*z*Cos((2*pf4*Pi*x*y*z)/L**3))/L**3 - (2*p2*pf2*Pi*Sin((2*pf2*Pi*y)/L))/L)/(-1 + kappa) + &
        &    ((d0 + d2*Cos((2*df2*Pi*y)/L) + d4*Cos((2*df4*Pi*x*y*z)/L**3) + d1*Sin((2*df1*Pi*x)/L) + d3*Sin((2*df3*Pi*z)/L))*&
        &       (2*((2*Pi*u2*uf2*Cos((2*Pi*uf2*y)/L))/L + (2*Pi*u4*uf4*x*z*Cos((2*Pi*uf4*x*y*z)/L**3))/L**3)*&
        &          (u0 + u1*Cos((2*Pi*uf1*x)/L) + u3*Cos((2*Pi*uf3*z)/L) + u2*Sin((2*Pi*uf2*y)/L) + u4*Sin((2*Pi*uf4*x*y*z)/L**3)) + &
        &         2*(v0 + v2*Cos((2*Pi*vf2*y)/L) + v4*Cos((2*Pi*vf4*x*y*z)/L**3) + v1*Sin((2*Pi*vf1*x)/L) + v3*Sin((2*Pi*vf3*z)/L))*&
        &          ((-2*Pi*v2*vf2*Sin((2*Pi*vf2*y)/L))/L - (2*Pi*v4*vf4*x*z*Sin((2*Pi*vf4*x*y*z)/L**3))/L**3) + &
        &         2*((2*Pi*w2*wf2*Cos((2*Pi*wf2*y)/L))/L + (2*Pi*w4*wf4*x*z*Cos((2*Pi*wf4*x*y*z)/L**3))/L**3)*&
        &          (w0 + w1*Cos((2*Pi*wf1*x)/L) + w3*Cos((2*Pi*wf3*z)/L) + w2*Sin((2*Pi*wf2*y)/L) + w4*Sin((2*Pi*wf4*x*y*z)/L**3))))/2. + &
        &    (((-2*d2*df2*Pi*Sin((2*df2*Pi*y)/L))/L - (2*d4*df4*Pi*x*z*Sin((2*df4*Pi*x*y*z)/L**3))/L**3)*&
        &       ((v0 + v2*Cos((2*Pi*vf2*y)/L) + v4*Cos((2*Pi*vf4*x*y*z)/L**3) + v1*Sin((2*Pi*vf1*x)/L) + v3*Sin((2*Pi*vf3*z)/L))**2 + &
        &         (u0 + u1*Cos((2*Pi*uf1*x)/L) + u3*Cos((2*Pi*uf3*z)/L) + u2*Sin((2*Pi*uf2*y)/L) + u4*Sin((2*Pi*uf4*x*y*z)/L**3))**2 + &
        &         (w0 + w1*Cos((2*Pi*wf1*x)/L) + w3*Cos((2*Pi*wf3*z)/L) + w2*Sin((2*Pi*wf2*y)/L) + w4*Sin((2*Pi*wf4*x*y*z)/L**3))**2))/2.);
         
    !! H function
    r(DENS_VAR) = r(DENS_VAR) + &
        & (d0 + d2*Cos((2*df2*Pi*y)/L) + d4*Cos((2*df4*Pi*x*y*z)/L**3) + d1*Sin((2*df1*Pi*x)/L) + d3*Sin((2*df3*Pi*z)/L))*((2*Pi*w4*wf4*x*y*Cos((2*Pi*wf4*x*y*z)/L**3))/L**3 - (2*Pi*w3*wf3*Sin((2*Pi*wf3*z)/L))/L) + &
        & ((2*d3*df3*Pi*Cos((2*df3*Pi*z)/L))/L - (2*d4*df4*Pi*x*y*Sin((2*df4*Pi*x*y*z)/L**3))/L**3)*(w0 + w1*Cos((2*Pi*wf1*x)/L) + w3*Cos((2*Pi*wf3*z)/L) + w2*Sin((2*Pi*wf2*y)/L) + w4*Sin((2*Pi*wf4*x*y*z)/L**3));
         
    r(MOMX_VAR) = r(MOMX_VAR) + &
        & (d0 + d2*Cos((2*df2*Pi*y)/L) + d4*Cos((2*df4*Pi*x*y*z)/L**3) + d1*Sin((2*df1*Pi*x)/L) + d3*Sin((2*df3*Pi*z)/L))*((2*Pi*w4*wf4*x*y*Cos((2*Pi*wf4*x*y*z)/L**3))/L**3 - (2*Pi*w3*wf3*Sin((2*Pi*wf3*z)/L))/L)*&
        &  (u0 + u1*Cos((2*Pi*uf1*x)/L) + u3*Cos((2*Pi*uf3*z)/L) + u2*Sin((2*Pi*uf2*y)/L) + u4*Sin((2*Pi*uf4*x*y*z)/L**3)) + &
        & (d0 + d2*Cos((2*df2*Pi*y)/L) + d4*Cos((2*df4*Pi*x*y*z)/L**3) + d1*Sin((2*df1*Pi*x)/L) + d3*Sin((2*df3*Pi*z)/L))*((2*Pi*u4*uf4*x*y*Cos((2*Pi*uf4*x*y*z)/L**3))/L**3 - (2*Pi*u3*uf3*Sin((2*Pi*uf3*z)/L))/L)*&
        &  (w0 + w1*Cos((2*Pi*wf1*x)/L) + w3*Cos((2*Pi*wf3*z)/L) + w2*Sin((2*Pi*wf2*y)/L) + w4*Sin((2*Pi*wf4*x*y*z)/L**3)) + &
        & ((2*d3*df3*Pi*Cos((2*df3*Pi*z)/L))/L - (2*d4*df4*Pi*x*y*Sin((2*df4*Pi*x*y*z)/L**3))/L**3)*(u0 + u1*Cos((2*Pi*uf1*x)/L) + u3*Cos((2*Pi*uf3*z)/L) + u2*Sin((2*Pi*uf2*y)/L) + u4*Sin((2*Pi*uf4*x*y*z)/L**3))*&
        &  (w0 + w1*Cos((2*Pi*wf1*x)/L) + w3*Cos((2*Pi*wf3*z)/L) + w2*Sin((2*Pi*wf2*y)/L) + w4*Sin((2*Pi*wf4*x*y*z)/L**3));
         
    r(MOMY_VAR) = r(MOMY_VAR) + &
        & (d0 + d2*Cos((2*df2*Pi*y)/L) + d4*Cos((2*df4*Pi*x*y*z)/L**3) + d1*Sin((2*df1*Pi*x)/L) + d3*Sin((2*df3*Pi*z)/L))*&
        &  (v0 + v2*Cos((2*Pi*vf2*y)/L) + v4*Cos((2*Pi*vf4*x*y*z)/L**3) + v1*Sin((2*Pi*vf1*x)/L) + v3*Sin((2*Pi*vf3*z)/L))*((2*Pi*w4*wf4*x*y*Cos((2*Pi*wf4*x*y*z)/L**3))/L**3 - (2*Pi*w3*wf3*Sin((2*Pi*wf3*z)/L))/L) + &
        & (v0 + v2*Cos((2*Pi*vf2*y)/L) + v4*Cos((2*Pi*vf4*x*y*z)/L**3) + v1*Sin((2*Pi*vf1*x)/L) + v3*Sin((2*Pi*vf3*z)/L))*((2*d3*df3*Pi*Cos((2*df3*Pi*z)/L))/L - (2*d4*df4*Pi*x*y*Sin((2*df4*Pi*x*y*z)/L**3))/L**3)*&
        &  (w0 + w1*Cos((2*Pi*wf1*x)/L) + w3*Cos((2*Pi*wf3*z)/L) + w2*Sin((2*Pi*wf2*y)/L) + w4*Sin((2*Pi*wf4*x*y*z)/L**3)) + &
        & (d0 + d2*Cos((2*df2*Pi*y)/L) + d4*Cos((2*df4*Pi*x*y*z)/L**3) + d1*Sin((2*df1*Pi*x)/L) + d3*Sin((2*df3*Pi*z)/L))*((2*Pi*v3*vf3*Cos((2*Pi*vf3*z)/L))/L - (2*Pi*v4*vf4*x*y*Sin((2*Pi*vf4*x*y*z)/L**3))/L**3)*&
        &  (w0 + w1*Cos((2*Pi*wf1*x)/L) + w3*Cos((2*Pi*wf3*z)/L) + w2*Sin((2*Pi*wf2*y)/L) + w4*Sin((2*Pi*wf4*x*y*z)/L**3));
         
    ! r(MOMZ_VAR) = r(MOMZ_VAR) + &
    !     & (2*p3*pf3*Pi*Cos((2*pf3*Pi*z)/L))/L + (2*p4*pf4*Pi*x*y*Cos((2*pf4*Pi*x*y*z)/L**3))/L**3 + &
    !     &  2*(d0 + d2*Cos((2*df2*Pi*y)/L) + d4*Cos((2*df4*Pi*x*y*z)/L**3) + d1*Sin((2*df1*Pi*x)/L) + d3*Sin((2*df3*Pi*z)/L))*((2*Pi*w4*wf4*x*y*Cos((2*Pi*wf4*x*y*z)/L**3))/L**3 - (2*Pi*w3*wf3*Sin((2*Pi*wf3*z)/L))/L)*&
    !     &   (w0 + w1*Cos((2*Pi*wf1*x)/L) + w3*Cos((2*Pi*wf3*z)/L) + w2*Sin((2*Pi*wf2*y)/L) + w4*Sin((2*Pi*wf4*x*y*z)/L**3)) + &
    !     &  ((2*d3*df3*Pi*Cos((2*df3*Pi*z)/L))/L - (2*d4*df4*Pi*x*y*Sin((2*df4*Pi*x*y*z)/L**3))/L**3)*(w0 + w1*Cos((2*Pi*wf1*x)/L) + w3*Cos((2*Pi*wf3*z)/L) + w2*Sin((2*Pi*wf2*y)/L) + w4*Sin((2*Pi*wf4*x*y*z)/L**3))**2;
         
    r(ENER_VAR) = r(ENER_VAR) + &
        & ((2*Pi*w4*wf4*x*y*Cos((2*Pi*wf4*x*y*z)/L**3))/L**3 - (2*Pi*w3*wf3*Sin((2*Pi*wf3*z)/L))/L)*&
        &  (p0 + p2*Cos((2*pf2*Pi*y)/L) + p1*Sin((2*pf1*Pi*x)/L) + p3*Sin((2*pf3*Pi*z)/L) + p4*Sin((2*pf4*Pi*x*y*z)/L**3) + &
        &    (p0 + p2*Cos((2*pf2*Pi*y)/L) + p1*Sin((2*pf1*Pi*x)/L) + p3*Sin((2*pf3*Pi*z)/L) + p4*Sin((2*pf4*Pi*x*y*z)/L**3))/(-1 + kappa) + &
        &    ((d0 + d2*Cos((2*df2*Pi*y)/L) + d4*Cos((2*df4*Pi*x*y*z)/L**3) + d1*Sin((2*df1*Pi*x)/L) + d3*Sin((2*df3*Pi*z)/L))*&
        &       ((v0 + v2*Cos((2*Pi*vf2*y)/L) + v4*Cos((2*Pi*vf4*x*y*z)/L**3) + v1*Sin((2*Pi*vf1*x)/L) + v3*Sin((2*Pi*vf3*z)/L))**2 + &
        &         (u0 + u1*Cos((2*Pi*uf1*x)/L) + u3*Cos((2*Pi*uf3*z)/L) + u2*Sin((2*Pi*uf2*y)/L) + u4*Sin((2*Pi*uf4*x*y*z)/L**3))**2 + &
        &         (w0 + w1*Cos((2*Pi*wf1*x)/L) + w3*Cos((2*Pi*wf3*z)/L) + w2*Sin((2*Pi*wf2*y)/L) + w4*Sin((2*Pi*wf4*x*y*z)/L**3))**2))/2.) + &
        & (w0 + w1*Cos((2*Pi*wf1*x)/L) + w3*Cos((2*Pi*wf3*z)/L) + w2*Sin((2*Pi*wf2*y)/L) + w4*Sin((2*Pi*wf4*x*y*z)/L**3))*&
        &  ((2*p3*pf3*Pi*Cos((2*pf3*Pi*z)/L))/L + (2*p4*pf4*Pi*x*y*Cos((2*pf4*Pi*x*y*z)/L**3))/L**3 + ((2*p3*pf3*Pi*Cos((2*pf3*Pi*z)/L))/L + (2*p4*pf4*Pi*x*y*Cos((2*pf4*Pi*x*y*z)/L**3))/L**3)/(-1 + kappa) + &
        &    ((d0 + d2*Cos((2*df2*Pi*y)/L) + d4*Cos((2*df4*Pi*x*y*z)/L**3) + d1*Sin((2*df1*Pi*x)/L) + d3*Sin((2*df3*Pi*z)/L))*&
        &       (2*((2*Pi*u4*uf4*x*y*Cos((2*Pi*uf4*x*y*z)/L**3))/L**3 - (2*Pi*u3*uf3*Sin((2*Pi*uf3*z)/L))/L)*&
        &          (u0 + u1*Cos((2*Pi*uf1*x)/L) + u3*Cos((2*Pi*uf3*z)/L) + u2*Sin((2*Pi*uf2*y)/L) + u4*Sin((2*Pi*uf4*x*y*z)/L**3)) + &
        &         2*(v0 + v2*Cos((2*Pi*vf2*y)/L) + v4*Cos((2*Pi*vf4*x*y*z)/L**3) + v1*Sin((2*Pi*vf1*x)/L) + v3*Sin((2*Pi*vf3*z)/L))*&
        &          ((2*Pi*v3*vf3*Cos((2*Pi*vf3*z)/L))/L - (2*Pi*v4*vf4*x*y*Sin((2*Pi*vf4*x*y*z)/L**3))/L**3) + &
        &         2*((2*Pi*w4*wf4*x*y*Cos((2*Pi*wf4*x*y*z)/L**3))/L**3 - (2*Pi*w3*wf3*Sin((2*Pi*wf3*z)/L))/L)*&
        &          (w0 + w1*Cos((2*Pi*wf1*x)/L) + w3*Cos((2*Pi*wf3*z)/L) + w2*Sin((2*Pi*wf2*y)/L) + w4*Sin((2*Pi*wf4*x*y*z)/L**3))))/2. + &
        &    (((2*d3*df3*Pi*Cos((2*df3*Pi*z)/L))/L - (2*d4*df4*Pi*x*y*Sin((2*df4*Pi*x*y*z)/L**3))/L**3)*&
        &       ((v0 + v2*Cos((2*Pi*vf2*y)/L) + v4*Cos((2*Pi*vf4*x*y*z)/L**3) + v1*Sin((2*Pi*vf1*x)/L) + v3*Sin((2*Pi*vf3*z)/L))**2 + &
        &         (u0 + u1*Cos((2*Pi*uf1*x)/L) + u3*Cos((2*Pi*uf3*z)/L) + u2*Sin((2*Pi*uf2*y)/L) + u4*Sin((2*Pi*uf4*x*y*z)/L**3))**2 + &
        &         (w0 + w1*Cos((2*Pi*wf1*x)/L) + w3*Cos((2*Pi*wf3*z)/L) + w2*Sin((2*Pi*wf2*y)/L) + w4*Sin((2*Pi*wf4*x*y*z)/L**3))**2))/2.);

end function

end module
