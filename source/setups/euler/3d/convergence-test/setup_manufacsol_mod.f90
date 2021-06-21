module setup_manufacsol_mod

    use constants_mod
    use setup_config_mod
    use equations_const_mod

contains

pure function init_manufacsol(t,x,y,z) result(u)

    real(dp), intent(in)    :: t,x,y,z
    real(dp)                :: u(N_VARS)

    !! u(DENS_VAR) = d0 + d1*Sin((2*Pi*df1*x)/L) + d2*Cos((2*Pi*df2*y)/L) + d3*Sin((2*Pi*df3*z)/L) + d4*Cos((2*Pi*df4*x*y*z)/L**3);
    !! u(MOMX_VAR) = u0 + u1*Cos((2*Pi*uf1*x)/L) + u2*Sin((2*Pi*uf2*y)/L) + u3*Cos((2*Pi*uf3*z)/L) + u4*Sin((2*Pi*uf4*x*y*z)/L**3);
    !! u(MOMY_VAR) = v0 + v1*Sin((2*Pi*vf1*x)/L) + v2*Cos((2*Pi*vf2*y)/L) + v3*Sin((2*Pi*vf3*z)/L) + v4*Cos((2*Pi*vf4*x*y*z)/L**3);
    !! u(MOMZ_VAR) = w0 + w1*Cos((2*Pi*wf1*x)/L) + w2*Sin((2*Pi*wf2*y)/L) + w3*Cos((2*Pi*wf3*z)/L) + w4*Sin((2*Pi*wf4*x*y*z)/L**3);
    !! u(ENER_VAR) = p0 + p1*Sin((2*Pi*pf1*x)/L) + p2*Cos((2*Pi*pf2*y)/L) + p3*Sin((2*Pi*pf3*z)/L) + p4*Sin((2*Pi*pf4*x*y*z)/L**3);

    ! u(DENS_VAR) = d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L)
    ! u(MOMX_VAR) = d0/10 + d1*sin(2*pi*(-t + x)/L)/10 + d2*cos(2*pi*(-t + y)/L)/10 + d3*sin(2*pi*(-t + z)/L)/10
    ! u(MOMY_VAR) = d0/5 + d1*sin(2*pi*(-t + x)/L)/5 + d2*cos(2*pi*(-t + y)/L)/5 + d3*sin(2*pi*(-t + z)/L)/5
    ! u(MOMZ_VAR) = 3*d0/10 + 3*d1*sin(2*pi*(-t + x)/L)/10 + 3*d2*cos(2*pi*(-t + y)/L)/10 + 3*d3*sin(2*pi*(-t + z)/L)/10
    ! u(ENER_VAR) = 7*d0/100 + 7*d1*sin(2*pi*(-t + x)/L)/100 + 7*d2*cos(2*pi*(-t + y)/L)/100 + 7*d3*sin(2*pi*(-t + z)/L)/100 + (p0 + p1*sin(2*pi*(-t + x)/L) + p2*cos(2*pi*(-t + y)/L) + p3*sin(2*pi*(-t + z)/L))/(kappa - 1)

    u(DENS_VAR) = d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L)
    u(MOMX_VAR) = u0*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))
    u(MOMY_VAR) = v0*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))
    u(MOMZ_VAR) = w0*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))
    u(ENER_VAR) = (u0**2 + v0**2 + w0**2)*(d0/2 + d1*sin(2*pi*(-t + x)/L)/2 + d2*cos(2*pi*(-t + y)/L)/2 + d3*sin(2*pi*(-t + z)/L)/2) + (p0 + p1*sin(2*pi*(-t + x)/L) + p2*cos(2*pi*(-t + y)/L) + p3*sin(2*pi*(-t + z)/L))/(kappa - 1)

end function

pure function calc_manufacsol(t,x,y,z) result(r)

    real(dp), intent(in)    :: t,x,y,z
    real(dp)                :: r(N_VARS)

    r = 0.0

    !! u_t !!

# if 1
    r(DENS_VAR) = r(DENS_VAR) + (-2*pi*d1*cos(2*pi*(-t + x)/L)/L + 2*pi*d2*sin(2*pi*(-t + y)/L)/L - 2*pi*d3*cos(2*pi*(-t + z)/L)/L)

    r(MOMX_VAR) = r(MOMX_VAR) + (u0*(-2*pi*d1*cos(2*pi*(-t + x)/L)/L + 2*pi*d2*sin(2*pi*(-t + y)/L)/L - 2*pi*d3*cos(2*pi*(-t + z)/L)/L))

    r(MOMY_VAR) = r(MOMY_VAR) + (v0*(-2*pi*d1*cos(2*pi*(-t + x)/L)/L + 2*pi*d2*sin(2*pi*(-t + y)/L)/L - 2*pi*d3*cos(2*pi*(-t + z)/L)/L))

    r(MOMZ_VAR) = r(MOMZ_VAR) + (w0*(-2*pi*d1*cos(2*pi*(-t + x)/L)/L + 2*pi*d2*sin(2*pi*(-t + y)/L)/L - 2*pi*d3*cos(2*pi*(-t + z)/L)/L))

    r(ENER_VAR) = r(ENER_VAR) + ((u0**2 + v0**2 + w0**2)*(-pi*d1*cos(2*pi*(-t + x)/L)/L + pi*d2*sin(2*pi*(-t + y)/L)/L - pi*d3*cos(2*pi*(-t + z)/L)/L) + (-2*pi*p1*cos(2*pi*(-t + x)/L)/L + 2*pi*p2*sin(2*pi*(-t + y)/L)/L - 2*pi*p3*cos(2*pi*(-t + z)/L)/L)/(kappa - 1))
# endif

    !! f_x !!

    r(DENS_VAR) = r(DENS_VAR) + (2*pi*d1*u0*cos(2*pi*(-t + x)/L)/L)

    r(MOMX_VAR) = r(MOMX_VAR) + ((kappa - 1)*(-0.5*(4*pi*d1*u0**2*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))*cos(2*pi*(-t + x)/L)/L + 4*pi*d1*v0**2*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))*cos(2*pi*(-t + x)/L)/L + 4*pi*d1*w0**2*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))*cos(2*pi*(-t + x)/L)/L)/(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L)) + pi*d1*(u0**2 + v0**2 + w0**2)*cos(2*pi*(-t + x)/L)/L + 1.0*pi*d1*(u0**2*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))**2 + v0**2*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))**2 + w0**2*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))**2)*cos(2*pi*(-t + x)/L)/(L*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))**2) + 2*pi*p1*cos(2*pi*(-t + x)/L)/(L*(kappa - 1))) + 2*pi*d1*u0**2*cos(2*pi*(-t + x)/L)/L)

    r(MOMY_VAR) = r(MOMY_VAR) + (2*pi*d1*u0*v0*cos(2*pi*(-t + x)/L)/L)

    r(MOMZ_VAR) = r(MOMZ_VAR) + (2*pi*d1*u0*w0*cos(2*pi*(-t + x)/L)/L)

    r(ENER_VAR) = r(ENER_VAR) + (u0*((kappa - 1)*(-0.5*(4*pi*d1*u0**2*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))*cos(2*pi*(-t + x)/L)/L + 4*pi*d1*v0**2*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))*cos(2*pi*(-t + x)/L)/L + 4*pi*d1*w0**2*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))*cos(2*pi*(-t + x)/L)/L)/(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L)) + pi*d1*(u0**2 + v0**2 + w0**2)*cos(2*pi*(-t + x)/L)/L + 1.0*pi*d1*(u0**2*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))**2 + v0**2*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))**2 + w0**2*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))**2)*cos(2*pi*(-t + x)/L)/(L*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))**2) + 2*pi*p1*cos(2*pi*(-t + x)/L)/(L*(kappa - 1))) + pi*d1*(u0**2 + v0**2 + w0**2)*cos(2*pi*(-t + x)/L)/L + 2*pi*p1*cos(2*pi*(-t + x)/L)/(L*(kappa - 1))))

    !! g_y !!

    r(DENS_VAR) = r(DENS_VAR) + (-2*pi*d2*v0*sin(2*pi*(-t + y)/L)/L)

    r(MOMX_VAR) = r(MOMX_VAR) + (-2*pi*d2*u0*v0*sin(2*pi*(-t + y)/L)/L)

    r(MOMY_VAR) = r(MOMY_VAR) + ((kappa - 1)*(-0.5*(-4*pi*d2*u0**2*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))*sin(2*pi*(-t + y)/L)/L - 4*pi*d2*v0**2*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))*sin(2*pi*(-t + y)/L)/L - 4*pi*d2*w0**2*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))*sin(2*pi*(-t + y)/L)/L)/(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L)) - pi*d2*(u0**2 + v0**2 + w0**2)*sin(2*pi*(-t + y)/L)/L - 1.0*pi*d2*(u0**2*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))**2 + v0**2*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))**2 + w0**2*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))**2)*sin(2*pi*(-t + y)/L)/(L*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))**2) - 2*pi*p2*sin(2*pi*(-t + y)/L)/(L*(kappa - 1))) - 2*pi*d2*v0**2*sin(2*pi*(-t + y)/L)/L)

    r(MOMZ_VAR) = r(MOMZ_VAR) + (-2*pi*d2*v0*w0*sin(2*pi*(-t + y)/L)/L)

    r(ENER_VAR) = r(ENER_VAR) + (v0*((kappa - 1)*(-0.5*(-4*pi*d2*u0**2*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))*sin(2*pi*(-t + y)/L)/L - 4*pi*d2*v0**2*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))*sin(2*pi*(-t + y)/L)/L - 4*pi*d2*w0**2*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))*sin(2*pi*(-t + y)/L)/L)/(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L)) - pi*d2*(u0**2 + v0**2 + w0**2)*sin(2*pi*(-t + y)/L)/L - 1.0*pi*d2*(u0**2*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))**2 + v0**2*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))**2 + w0**2*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))**2)*sin(2*pi*(-t + y)/L)/(L*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))**2) - 2*pi*p2*sin(2*pi*(-t + y)/L)/(L*(kappa - 1))) - pi*d2*(u0**2 + v0**2 + w0**2)*sin(2*pi*(-t + y)/L)/L - 2*pi*p2*sin(2*pi*(-t + y)/L)/(L*(kappa - 1))))

    !! h_z !!

    r(DENS_VAR) = r(DENS_VAR) + (2*pi*d3*w0*cos(2*pi*(-t + z)/L)/L)

    r(MOMX_VAR) = r(MOMX_VAR) + (2*pi*d3*u0*w0*cos(2*pi*(-t + z)/L)/L)

    r(MOMY_VAR) = r(MOMY_VAR) + (2*pi*d3*v0*w0*cos(2*pi*(-t + z)/L)/L)

    r(MOMZ_VAR) = r(MOMZ_VAR) + ((kappa - 1)*(-0.5*(4*pi*d3*u0**2*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))*cos(2*pi*(-t + z)/L)/L + 4*pi*d3*v0**2*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))*cos(2*pi*(-t + z)/L)/L + 4*pi*d3*w0**2*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))*cos(2*pi*(-t + z)/L)/L)/(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L)) + pi*d3*(u0**2 + v0**2 + w0**2)*cos(2*pi*(-t + z)/L)/L + 1.0*pi*d3*(u0**2*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))**2 + v0**2*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))**2 + w0**2*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))**2)*cos(2*pi*(-t + z)/L)/(L*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))**2) + 2*pi*p3*cos(2*pi*(-t + z)/L)/(L*(kappa - 1))) + 2*pi*d3*w0**2*cos(2*pi*(-t + z)/L)/L)

    r(ENER_VAR) = r(ENER_VAR) + (w0*((kappa - 1)*(-0.5*(4*pi*d3*u0**2*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))*cos(2*pi*(-t + z)/L)/L + 4*pi*d3*v0**2*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))*cos(2*pi*(-t + z)/L)/L + 4*pi*d3*w0**2*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))*cos(2*pi*(-t + z)/L)/L)/(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L)) + pi*d3*(u0**2 + v0**2 + w0**2)*cos(2*pi*(-t + z)/L)/L + 1.0*pi*d3*(u0**2*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))**2 + v0**2*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))**2 + w0**2*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))**2)*cos(2*pi*(-t + z)/L)/(L*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))**2) + 2*pi*p3*cos(2*pi*(-t + z)/L)/(L*(kappa - 1))) + pi*d3*(u0**2 + v0**2 + w0**2)*cos(2*pi*(-t + z)/L)/L + 2*pi*p3*cos(2*pi*(-t + z)/L)/(L*(kappa - 1))))

    !! !! u_t !!

    !! r(DENS_VAR) = r(DENS_VAR) + (-2*pi*(d1*cos(2*pi*(t - x)/L) + d2*sin(2*pi*(t - y)/L) + d3*cos(2*pi*(t - z)/L))/L)

    !! r(MOMX_VAR) = r(MOMX_VAR) + (-pi*(d1*cos(2*pi*(t - x)/L) + d2*sin(2*pi*(t - y)/L) + d3*cos(2*pi*(t - z)/L))/(5*L))

    !! r(MOMY_VAR) = r(MOMY_VAR) + (-2*pi*(d1*cos(2*pi*(t - x)/L) + d2*sin(2*pi*(t - y)/L) + d3*cos(2*pi*(t - z)/L))/(5*L))

    !! r(MOMZ_VAR) = r(MOMZ_VAR) + (-3*pi*(d1*cos(2*pi*(t - x)/L) + d2*sin(2*pi*(t - y)/L) + d3*cos(2*pi*(t - z)/L))/(5*L))

    !! r(ENER_VAR) = r(ENER_VAR) + (pi*(-100*p1*cos(2*pi*(t - x)/L) - 100*p2*sin(2*pi*(t - y)/L) - 100*p3*cos(2*pi*(t - z)/L) + (7 - 7*kappa)*(d1*cos(2*pi*(t - x)/L) + d2*sin(2*pi*(t - y)/L) + d3*cos(2*pi*(t - z)/L)))/(50*L*(kappa - 1)))

    !! !! f_x !!

    !! r(DENS_VAR) = r(DENS_VAR) + (pi*d1*cos(2*pi*(-t + x)/L)/(5*L))

    !! r(MOMX_VAR) = r(MOMX_VAR) + ((kappa - 1)*(-0.5*(2*pi*d1*(d0/10 + d1*sin(2*pi*(-t + x)/L)/10 + d2*cos(2*pi*(-t + y)/L)/10 + d3*sin(2*pi*(-t + z)/L)/10)*cos(2*pi*(-t + x)/L)/(5*L) + 4*pi*d1*(d0/5 + d1*sin(2*pi*(-t + x)/L)/5 + d2*cos(2*pi*(-t + y)/L)/5 + d3*sin(2*pi*(-t + z)/L)/5)*cos(2*pi*(-t + x)/L)/(5*L) + 6*pi*d1*(3*d0/10 + 3*d1*sin(2*pi*(-t + x)/L)/10 + 3*d2*cos(2*pi*(-t + y)/L)/10 + 3*d3*sin(2*pi*(-t + z)/L)/10)*cos(2*pi*(-t + x)/L)/(5*L))/(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L)) + 1.0*pi*d1*((d0/10 + d1*sin(2*pi*(-t + x)/L)/10 + d2*cos(2*pi*(-t + y)/L)/10 + d3*sin(2*pi*(-t + z)/L)/10)**2 + (d0/5 + d1*sin(2*pi*(-t + x)/L)/5 + d2*cos(2*pi*(-t + y)/L)/5 + d3*sin(2*pi*(-t + z)/L)/5)**2 + (3*d0/10 + 3*d1*sin(2*pi*(-t + x)/L)/10 + 3*d2*cos(2*pi*(-t + y)/L)/10 + 3*d3*sin(2*pi*(-t + z)/L)/10)**2)*cos(2*pi*(-t + x)/L)/(L*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))**2) + 7*pi*d1*cos(2*pi*(-t + x)/L)/(50*L) + 2*pi*p1*cos(2*pi*(-t + x)/L)/(L*(kappa - 1))) - 2*pi*d1*(d0/10 + d1*sin(2*pi*(-t + x)/L)/10 + d2*cos(2*pi*(-t + y)/L)/10 + d3*sin(2*pi*(-t + z)/L)/10)**2*cos(2*pi*(-t + x)/L)/(L*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))**2) + 2*pi*d1*(d0/10 + d1*sin(2*pi*(-t + x)/L)/10 + d2*cos(2*pi*(-t + y)/L)/10 + d3*sin(2*pi*(-t + z)/L)/10)*cos(2*pi*(-t + x)/L)/(5*L*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))))

    !! r(MOMY_VAR) = r(MOMY_VAR) + (-2*pi*d1*(d0/10 + d1*sin(2*pi*(-t + x)/L)/10 + d2*cos(2*pi*(-t + y)/L)/10 + d3*sin(2*pi*(-t + z)/L)/10)*(d0/5 + d1*sin(2*pi*(-t + x)/L)/5 + d2*cos(2*pi*(-t + y)/L)/5 + d3*sin(2*pi*(-t + z)/L)/5)*cos(2*pi*(-t + x)/L)/(L*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))**2) + 2*pi*d1*(d0/10 + d1*sin(2*pi*(-t + x)/L)/10 + d2*cos(2*pi*(-t + y)/L)/10 + d3*sin(2*pi*(-t + z)/L)/10)*cos(2*pi*(-t + x)/L)/(5*L*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))) + pi*d1*(d0/5 + d1*sin(2*pi*(-t + x)/L)/5 + d2*cos(2*pi*(-t + y)/L)/5 + d3*sin(2*pi*(-t + z)/L)/5)*cos(2*pi*(-t + x)/L)/(5*L*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))))

    !! r(MOMZ_VAR) = r(MOMZ_VAR) + (-2*pi*d1*(d0/10 + d1*sin(2*pi*(-t + x)/L)/10 + d2*cos(2*pi*(-t + y)/L)/10 + d3*sin(2*pi*(-t + z)/L)/10)*(3*d0/10 + 3*d1*sin(2*pi*(-t + x)/L)/10 + 3*d2*cos(2*pi*(-t + y)/L)/10 + 3*d3*sin(2*pi*(-t + z)/L)/10)*cos(2*pi*(-t + x)/L)/(L*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))**2) + 3*pi*d1*(d0/10 + d1*sin(2*pi*(-t + x)/L)/10 + d2*cos(2*pi*(-t + y)/L)/10 + d3*sin(2*pi*(-t + z)/L)/10)*cos(2*pi*(-t + x)/L)/(5*L*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))) + pi*d1*(3*d0/10 + 3*d1*sin(2*pi*(-t + x)/L)/10 + 3*d2*cos(2*pi*(-t + y)/L)/10 + 3*d3*sin(2*pi*(-t + z)/L)/10)*cos(2*pi*(-t + x)/L)/(5*L*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))))

    !! r(ENER_VAR) = r(ENER_VAR) + (((kappa - 1)*(-0.5*(2*pi*d1*(d0/10 + d1*sin(2*pi*(-t + x)/L)/10 + d2*cos(2*pi*(-t + y)/L)/10 + d3*sin(2*pi*(-t + z)/L)/10)*cos(2*pi*(-t + x)/L)/(5*L) + 4*pi*d1*(d0/5 + d1*sin(2*pi*(-t + x)/L)/5 + d2*cos(2*pi*(-t + y)/L)/5 + d3*sin(2*pi*(-t + z)/L)/5)*cos(2*pi*(-t + x)/L)/(5*L) + 6*pi*d1*(3*d0/10 + 3*d1*sin(2*pi*(-t + x)/L)/10 + 3*d2*cos(2*pi*(-t + y)/L)/10 + 3*d3*sin(2*pi*(-t + z)/L)/10)*cos(2*pi*(-t + x)/L)/(5*L))/(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L)) + 1.0*pi*d1*((d0/10 + d1*sin(2*pi*(-t + x)/L)/10 + d2*cos(2*pi*(-t + y)/L)/10 + d3*sin(2*pi*(-t + z)/L)/10)**2 + (d0/5 + d1*sin(2*pi*(-t + x)/L)/5 + d2*cos(2*pi*(-t + y)/L)/5 + d3*sin(2*pi*(-t + z)/L)/5)**2 + (3*d0/10 + 3*d1*sin(2*pi*(-t + x)/L)/10 + 3*d2*cos(2*pi*(-t + y)/L)/10 + 3*d3*sin(2*pi*(-t + z)/L)/10)**2)*cos(2*pi*(-t + x)/L)/(L*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))**2) + 7*pi*d1*cos(2*pi*(-t + x)/L)/(50*L) + 2*pi*p1*cos(2*pi*(-t + x)/L)/(L*(kappa - 1))) + 7*pi*d1*cos(2*pi*(-t + x)/L)/(50*L) + 2*pi*p1*cos(2*pi*(-t + x)/L)/(L*(kappa - 1)))*(d0/10 + d1*sin(2*pi*(-t + x)/L)/10 + d2*cos(2*pi*(-t + y)/L)/10 + d3*sin(2*pi*(-t + z)/L)/10)/(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L)) - 2*pi*d1*(d0/10 + d1*sin(2*pi*(-t + x)/L)/10 + d2*cos(2*pi*(-t + y)/L)/10 + d3*sin(2*pi*(-t + z)/L)/10)*(7*d0/100 + 7*d1*sin(2*pi*(-t + x)/L)/100 + 7*d2*cos(2*pi*(-t + y)/L)/100 + 7*d3*sin(2*pi*(-t + z)/L)/100 + (kappa - 1)*(7*d0/100 + 7*d1*sin(2*pi*(-t + x)/L)/100 + 7*d2*cos(2*pi*(-t + y)/L)/100 + 7*d3*sin(2*pi*(-t + z)/L)/100 - 0.5*((d0/10 + d1*sin(2*pi*(-t + x)/L)/10 + d2*cos(2*pi*(-t + y)/L)/10 + d3*sin(2*pi*(-t + z)/L)/10)**2 + (d0/5 + d1*sin(2*pi*(-t + x)/L)/5 + d2*cos(2*pi*(-t + y)/L)/5 + d3*sin(2*pi*(-t + z)/L)/5)**2 + (3*d0/10 + 3*d1*sin(2*pi*(-t + x)/L)/10 + 3*d2*cos(2*pi*(-t + y)/L)/10 + 3*d3*sin(2*pi*(-t + z)/L)/10)**2)/(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L)) + (p0 + p1*sin(2*pi*(-t + x)/L) + p2*cos(2*pi*(-t + y)/L) + p3*sin(2*pi*(-t + z)/L))/(kappa - 1)) + (p0 + p1*sin(2*pi*(-t + x)/L) + p2*cos(2*pi*(-t + y)/L) + p3*sin(2*pi*(-t + z)/L))/(kappa - 1))*cos(2*pi*(-t + x)/L)/(L*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))**2) + pi*d1*(7*d0/100 + 7*d1*sin(2*pi*(-t + x)/L)/100 + 7*d2*cos(2*pi*(-t + y)/L)/100 + 7*d3*sin(2*pi*(-t + z)/L)/100 + (kappa - 1)*(7*d0/100 + 7*d1*sin(2*pi*(-t + x)/L)/100 + 7*d2*cos(2*pi*(-t + y)/L)/100 + 7*d3*sin(2*pi*(-t + z)/L)/100 - 0.5*((d0/10 + d1*sin(2*pi*(-t + x)/L)/10 + d2*cos(2*pi*(-t + y)/L)/10 + d3*sin(2*pi*(-t + z)/L)/10)**2 + (d0/5 + d1*sin(2*pi*(-t + x)/L)/5 + d2*cos(2*pi*(-t + y)/L)/5 + d3*sin(2*pi*(-t + z)/L)/5)**2 + (3*d0/10 + 3*d1*sin(2*pi*(-t + x)/L)/10 + 3*d2*cos(2*pi*(-t + y)/L)/10 + 3*d3*sin(2*pi*(-t + z)/L)/10)**2)/(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L)) + (p0 + p1*sin(2*pi*(-t + x)/L) + p2*cos(2*pi*(-t + y)/L) + p3*sin(2*pi*(-t + z)/L))/(kappa - 1)) + (p0 + p1*sin(2*pi*(-t + x)/L) + p2*cos(2*pi*(-t + y)/L) + p3*sin(2*pi*(-t + z)/L))/(kappa - 1))*cos(2*pi*(-t + x)/L)/(5*L*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))))

    !! !! g_y !!

    !! r(DENS_VAR) = r(DENS_VAR) + (-2*pi*d2*sin(2*pi*(-t + y)/L)/(5*L))

    !! r(MOMX_VAR) = r(MOMX_VAR) + (2*pi*d2*(d0/10 + d1*sin(2*pi*(-t + x)/L)/10 + d2*cos(2*pi*(-t + y)/L)/10 + d3*sin(2*pi*(-t + z)/L)/10)*(d0/5 + d1*sin(2*pi*(-t + x)/L)/5 + d2*cos(2*pi*(-t + y)/L)/5 + d3*sin(2*pi*(-t + z)/L)/5)*sin(2*pi*(-t + y)/L)/(L*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))**2) - 2*pi*d2*(d0/10 + d1*sin(2*pi*(-t + x)/L)/10 + d2*cos(2*pi*(-t + y)/L)/10 + d3*sin(2*pi*(-t + z)/L)/10)*sin(2*pi*(-t + y)/L)/(5*L*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))) - pi*d2*(d0/5 + d1*sin(2*pi*(-t + x)/L)/5 + d2*cos(2*pi*(-t + y)/L)/5 + d3*sin(2*pi*(-t + z)/L)/5)*sin(2*pi*(-t + y)/L)/(5*L*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))))

    !! r(MOMY_VAR) = r(MOMY_VAR) + ((kappa - 1)*(-0.5*(-2*pi*d2*(d0/10 + d1*sin(2*pi*(-t + x)/L)/10 + d2*cos(2*pi*(-t + y)/L)/10 + d3*sin(2*pi*(-t + z)/L)/10)*sin(2*pi*(-t + y)/L)/(5*L) - 4*pi*d2*(d0/5 + d1*sin(2*pi*(-t + x)/L)/5 + d2*cos(2*pi*(-t + y)/L)/5 + d3*sin(2*pi*(-t + z)/L)/5)*sin(2*pi*(-t + y)/L)/(5*L) - 6*pi*d2*(3*d0/10 + 3*d1*sin(2*pi*(-t + x)/L)/10 + 3*d2*cos(2*pi*(-t + y)/L)/10 + 3*d3*sin(2*pi*(-t + z)/L)/10)*sin(2*pi*(-t + y)/L)/(5*L))/(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L)) - 1.0*pi*d2*((d0/10 + d1*sin(2*pi*(-t + x)/L)/10 + d2*cos(2*pi*(-t + y)/L)/10 + d3*sin(2*pi*(-t + z)/L)/10)**2 + (d0/5 + d1*sin(2*pi*(-t + x)/L)/5 + d2*cos(2*pi*(-t + y)/L)/5 + d3*sin(2*pi*(-t + z)/L)/5)**2 + (3*d0/10 + 3*d1*sin(2*pi*(-t + x)/L)/10 + 3*d2*cos(2*pi*(-t + y)/L)/10 + 3*d3*sin(2*pi*(-t + z)/L)/10)**2)*sin(2*pi*(-t + y)/L)/(L*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))**2) - 7*pi*d2*sin(2*pi*(-t + y)/L)/(50*L) - 2*pi*p2*sin(2*pi*(-t + y)/L)/(L*(kappa - 1))) + 2*pi*d2*(d0/5 + d1*sin(2*pi*(-t + x)/L)/5 + d2*cos(2*pi*(-t + y)/L)/5 + d3*sin(2*pi*(-t + z)/L)/5)**2*sin(2*pi*(-t + y)/L)/(L*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))**2) - 4*pi*d2*(d0/5 + d1*sin(2*pi*(-t + x)/L)/5 + d2*cos(2*pi*(-t + y)/L)/5 + d3*sin(2*pi*(-t + z)/L)/5)*sin(2*pi*(-t + y)/L)/(5*L*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))))

    !! r(MOMZ_VAR) = r(MOMZ_VAR) + (2*pi*d2*(d0/5 + d1*sin(2*pi*(-t + x)/L)/5 + d2*cos(2*pi*(-t + y)/L)/5 + d3*sin(2*pi*(-t + z)/L)/5)*(3*d0/10 + 3*d1*sin(2*pi*(-t + x)/L)/10 + 3*d2*cos(2*pi*(-t + y)/L)/10 + 3*d3*sin(2*pi*(-t + z)/L)/10)*sin(2*pi*(-t + y)/L)/(L*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))**2) - 3*pi*d2*(d0/5 + d1*sin(2*pi*(-t + x)/L)/5 + d2*cos(2*pi*(-t + y)/L)/5 + d3*sin(2*pi*(-t + z)/L)/5)*sin(2*pi*(-t + y)/L)/(5*L*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))) - 2*pi*d2*(3*d0/10 + 3*d1*sin(2*pi*(-t + x)/L)/10 + 3*d2*cos(2*pi*(-t + y)/L)/10 + 3*d3*sin(2*pi*(-t + z)/L)/10)*sin(2*pi*(-t + y)/L)/(5*L*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))))

    !! r(ENER_VAR) = r(ENER_VAR) + (((kappa - 1)*(-0.5*(-2*pi*d2*(d0/10 + d1*sin(2*pi*(-t + x)/L)/10 + d2*cos(2*pi*(-t + y)/L)/10 + d3*sin(2*pi*(-t + z)/L)/10)*sin(2*pi*(-t + y)/L)/(5*L) - 4*pi*d2*(d0/5 + d1*sin(2*pi*(-t + x)/L)/5 + d2*cos(2*pi*(-t + y)/L)/5 + d3*sin(2*pi*(-t + z)/L)/5)*sin(2*pi*(-t + y)/L)/(5*L) - 6*pi*d2*(3*d0/10 + 3*d1*sin(2*pi*(-t + x)/L)/10 + 3*d2*cos(2*pi*(-t + y)/L)/10 + 3*d3*sin(2*pi*(-t + z)/L)/10)*sin(2*pi*(-t + y)/L)/(5*L))/(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L)) - 1.0*pi*d2*((d0/10 + d1*sin(2*pi*(-t + x)/L)/10 + d2*cos(2*pi*(-t + y)/L)/10 + d3*sin(2*pi*(-t + z)/L)/10)**2 + (d0/5 + d1*sin(2*pi*(-t + x)/L)/5 + d2*cos(2*pi*(-t + y)/L)/5 + d3*sin(2*pi*(-t + z)/L)/5)**2 + (3*d0/10 + 3*d1*sin(2*pi*(-t + x)/L)/10 + 3*d2*cos(2*pi*(-t + y)/L)/10 + 3*d3*sin(2*pi*(-t + z)/L)/10)**2)*sin(2*pi*(-t + y)/L)/(L*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))**2) - 7*pi*d2*sin(2*pi*(-t + y)/L)/(50*L) - 2*pi*p2*sin(2*pi*(-t + y)/L)/(L*(kappa - 1))) - 7*pi*d2*sin(2*pi*(-t + y)/L)/(50*L) - 2*pi*p2*sin(2*pi*(-t + y)/L)/(L*(kappa - 1)))*(d0/5 + d1*sin(2*pi*(-t + x)/L)/5 + d2*cos(2*pi*(-t + y)/L)/5 + d3*sin(2*pi*(-t + z)/L)/5)/(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L)) + 2*pi*d2*(d0/5 + d1*sin(2*pi*(-t + x)/L)/5 + d2*cos(2*pi*(-t + y)/L)/5 + d3*sin(2*pi*(-t + z)/L)/5)*(7*d0/100 + 7*d1*sin(2*pi*(-t + x)/L)/100 + 7*d2*cos(2*pi*(-t + y)/L)/100 + 7*d3*sin(2*pi*(-t + z)/L)/100 + (kappa - 1)*(7*d0/100 + 7*d1*sin(2*pi*(-t + x)/L)/100 + 7*d2*cos(2*pi*(-t + y)/L)/100 + 7*d3*sin(2*pi*(-t + z)/L)/100 - 0.5*((d0/10 + d1*sin(2*pi*(-t + x)/L)/10 + d2*cos(2*pi*(-t + y)/L)/10 + d3*sin(2*pi*(-t + z)/L)/10)**2 + (d0/5 + d1*sin(2*pi*(-t + x)/L)/5 + d2*cos(2*pi*(-t + y)/L)/5 + d3*sin(2*pi*(-t + z)/L)/5)**2 + (3*d0/10 + 3*d1*sin(2*pi*(-t + x)/L)/10 + 3*d2*cos(2*pi*(-t + y)/L)/10 + 3*d3*sin(2*pi*(-t + z)/L)/10)**2)/(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L)) + (p0 + p1*sin(2*pi*(-t + x)/L) + p2*cos(2*pi*(-t + y)/L) + p3*sin(2*pi*(-t + z)/L))/(kappa - 1)) + (p0 + p1*sin(2*pi*(-t + x)/L) + p2*cos(2*pi*(-t + y)/L) + p3*sin(2*pi*(-t + z)/L))/(kappa - 1))*sin(2*pi*(-t + y)/L)/(L*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))**2) - 2*pi*d2*(7*d0/100 + 7*d1*sin(2*pi*(-t + x)/L)/100 + 7*d2*cos(2*pi*(-t + y)/L)/100 + 7*d3*sin(2*pi*(-t + z)/L)/100 + (kappa - 1)*(7*d0/100 + 7*d1*sin(2*pi*(-t + x)/L)/100 + 7*d2*cos(2*pi*(-t + y)/L)/100 + 7*d3*sin(2*pi*(-t + z)/L)/100 - 0.5*((d0/10 + d1*sin(2*pi*(-t + x)/L)/10 + d2*cos(2*pi*(-t + y)/L)/10 + d3*sin(2*pi*(-t + z)/L)/10)**2 + (d0/5 + d1*sin(2*pi*(-t + x)/L)/5 + d2*cos(2*pi*(-t + y)/L)/5 + d3*sin(2*pi*(-t + z)/L)/5)**2 + (3*d0/10 + 3*d1*sin(2*pi*(-t + x)/L)/10 + 3*d2*cos(2*pi*(-t + y)/L)/10 + 3*d3*sin(2*pi*(-t + z)/L)/10)**2)/(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L)) + (p0 + p1*sin(2*pi*(-t + x)/L) + p2*cos(2*pi*(-t + y)/L) + p3*sin(2*pi*(-t + z)/L))/(kappa - 1)) + (p0 + p1*sin(2*pi*(-t + x)/L) + p2*cos(2*pi*(-t + y)/L) + p3*sin(2*pi*(-t + z)/L))/(kappa - 1))*sin(2*pi*(-t + y)/L)/(5*L*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))))

    !! !! h_z !!

    !! r(DENS_VAR) = r(DENS_VAR) + (3*pi*d3*cos(2*pi*(-t + z)/L)/(5*L))

    !! r(MOMX_VAR) = r(MOMX_VAR) + (-2*pi*d3*(d0/10 + d1*sin(2*pi*(-t + x)/L)/10 + d2*cos(2*pi*(-t + y)/L)/10 + d3*sin(2*pi*(-t + z)/L)/10)*(3*d0/10 + 3*d1*sin(2*pi*(-t + x)/L)/10 + 3*d2*cos(2*pi*(-t + y)/L)/10 + 3*d3*sin(2*pi*(-t + z)/L)/10)*cos(2*pi*(-t + z)/L)/(L*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))**2) + 3*pi*d3*(d0/10 + d1*sin(2*pi*(-t + x)/L)/10 + d2*cos(2*pi*(-t + y)/L)/10 + d3*sin(2*pi*(-t + z)/L)/10)*cos(2*pi*(-t + z)/L)/(5*L*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))) + pi*d3*(3*d0/10 + 3*d1*sin(2*pi*(-t + x)/L)/10 + 3*d2*cos(2*pi*(-t + y)/L)/10 + 3*d3*sin(2*pi*(-t + z)/L)/10)*cos(2*pi*(-t + z)/L)/(5*L*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))))

    !! r(MOMY_VAR) = r(MOMY_VAR) + (-2*pi*d3*(d0/5 + d1*sin(2*pi*(-t + x)/L)/5 + d2*cos(2*pi*(-t + y)/L)/5 + d3*sin(2*pi*(-t + z)/L)/5)*(3*d0/10 + 3*d1*sin(2*pi*(-t + x)/L)/10 + 3*d2*cos(2*pi*(-t + y)/L)/10 + 3*d3*sin(2*pi*(-t + z)/L)/10)*cos(2*pi*(-t + z)/L)/(L*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))**2) + 3*pi*d3*(d0/5 + d1*sin(2*pi*(-t + x)/L)/5 + d2*cos(2*pi*(-t + y)/L)/5 + d3*sin(2*pi*(-t + z)/L)/5)*cos(2*pi*(-t + z)/L)/(5*L*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))) + 2*pi*d3*(3*d0/10 + 3*d1*sin(2*pi*(-t + x)/L)/10 + 3*d2*cos(2*pi*(-t + y)/L)/10 + 3*d3*sin(2*pi*(-t + z)/L)/10)*cos(2*pi*(-t + z)/L)/(5*L*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))))

    !! r(MOMZ_VAR) = r(MOMZ_VAR) + ((kappa - 1)*(-0.5*(2*pi*d3*(d0/10 + d1*sin(2*pi*(-t + x)/L)/10 + d2*cos(2*pi*(-t + y)/L)/10 + d3*sin(2*pi*(-t + z)/L)/10)*cos(2*pi*(-t + z)/L)/(5*L) + 4*pi*d3*(d0/5 + d1*sin(2*pi*(-t + x)/L)/5 + d2*cos(2*pi*(-t + y)/L)/5 + d3*sin(2*pi*(-t + z)/L)/5)*cos(2*pi*(-t + z)/L)/(5*L) + 6*pi*d3*(3*d0/10 + 3*d1*sin(2*pi*(-t + x)/L)/10 + 3*d2*cos(2*pi*(-t + y)/L)/10 + 3*d3*sin(2*pi*(-t + z)/L)/10)*cos(2*pi*(-t + z)/L)/(5*L))/(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L)) + 1.0*pi*d3*((d0/10 + d1*sin(2*pi*(-t + x)/L)/10 + d2*cos(2*pi*(-t + y)/L)/10 + d3*sin(2*pi*(-t + z)/L)/10)**2 + (d0/5 + d1*sin(2*pi*(-t + x)/L)/5 + d2*cos(2*pi*(-t + y)/L)/5 + d3*sin(2*pi*(-t + z)/L)/5)**2 + (3*d0/10 + 3*d1*sin(2*pi*(-t + x)/L)/10 + 3*d2*cos(2*pi*(-t + y)/L)/10 + 3*d3*sin(2*pi*(-t + z)/L)/10)**2)*cos(2*pi*(-t + z)/L)/(L*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))**2) + 7*pi*d3*cos(2*pi*(-t + z)/L)/(50*L) + 2*pi*p3*cos(2*pi*(-t + z)/L)/(L*(kappa - 1))) - 2*pi*d3*(3*d0/10 + 3*d1*sin(2*pi*(-t + x)/L)/10 + 3*d2*cos(2*pi*(-t + y)/L)/10 + 3*d3*sin(2*pi*(-t + z)/L)/10)**2*cos(2*pi*(-t + z)/L)/(L*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))**2) + 6*pi*d3*(3*d0/10 + 3*d1*sin(2*pi*(-t + x)/L)/10 + 3*d2*cos(2*pi*(-t + y)/L)/10 + 3*d3*sin(2*pi*(-t + z)/L)/10)*cos(2*pi*(-t + z)/L)/(5*L*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))))

    !! r(ENER_VAR) = r(ENER_VAR) + (((kappa - 1)*(-0.5*(2*pi*d3*(d0/10 + d1*sin(2*pi*(-t + x)/L)/10 + d2*cos(2*pi*(-t + y)/L)/10 + d3*sin(2*pi*(-t + z)/L)/10)*cos(2*pi*(-t + z)/L)/(5*L) + 4*pi*d3*(d0/5 + d1*sin(2*pi*(-t + x)/L)/5 + d2*cos(2*pi*(-t + y)/L)/5 + d3*sin(2*pi*(-t + z)/L)/5)*cos(2*pi*(-t + z)/L)/(5*L) + 6*pi*d3*(3*d0/10 + 3*d1*sin(2*pi*(-t + x)/L)/10 + 3*d2*cos(2*pi*(-t + y)/L)/10 + 3*d3*sin(2*pi*(-t + z)/L)/10)*cos(2*pi*(-t + z)/L)/(5*L))/(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L)) + 1.0*pi*d3*((d0/10 + d1*sin(2*pi*(-t + x)/L)/10 + d2*cos(2*pi*(-t + y)/L)/10 + d3*sin(2*pi*(-t + z)/L)/10)**2 + (d0/5 + d1*sin(2*pi*(-t + x)/L)/5 + d2*cos(2*pi*(-t + y)/L)/5 + d3*sin(2*pi*(-t + z)/L)/5)**2 + (3*d0/10 + 3*d1*sin(2*pi*(-t + x)/L)/10 + 3*d2*cos(2*pi*(-t + y)/L)/10 + 3*d3*sin(2*pi*(-t + z)/L)/10)**2)*cos(2*pi*(-t + z)/L)/(L*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))**2) + 7*pi*d3*cos(2*pi*(-t + z)/L)/(50*L) + 2*pi*p3*cos(2*pi*(-t + z)/L)/(L*(kappa - 1))) + 7*pi*d3*cos(2*pi*(-t + z)/L)/(50*L) + 2*pi*p3*cos(2*pi*(-t + z)/L)/(L*(kappa - 1)))*(3*d0/10 + 3*d1*sin(2*pi*(-t + x)/L)/10 + 3*d2*cos(2*pi*(-t + y)/L)/10 + 3*d3*sin(2*pi*(-t + z)/L)/10)/(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L)) - 2*pi*d3*(3*d0/10 + 3*d1*sin(2*pi*(-t + x)/L)/10 + 3*d2*cos(2*pi*(-t + y)/L)/10 + 3*d3*sin(2*pi*(-t + z)/L)/10)*(7*d0/100 + 7*d1*sin(2*pi*(-t + x)/L)/100 + 7*d2*cos(2*pi*(-t + y)/L)/100 + 7*d3*sin(2*pi*(-t + z)/L)/100 + (kappa - 1)*(7*d0/100 + 7*d1*sin(2*pi*(-t + x)/L)/100 + 7*d2*cos(2*pi*(-t + y)/L)/100 + 7*d3*sin(2*pi*(-t + z)/L)/100 - 0.5*((d0/10 + d1*sin(2*pi*(-t + x)/L)/10 + d2*cos(2*pi*(-t + y)/L)/10 + d3*sin(2*pi*(-t + z)/L)/10)**2 + (d0/5 + d1*sin(2*pi*(-t + x)/L)/5 + d2*cos(2*pi*(-t + y)/L)/5 + d3*sin(2*pi*(-t + z)/L)/5)**2 + (3*d0/10 + 3*d1*sin(2*pi*(-t + x)/L)/10 + 3*d2*cos(2*pi*(-t + y)/L)/10 + 3*d3*sin(2*pi*(-t + z)/L)/10)**2)/(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L)) + (p0 + p1*sin(2*pi*(-t + x)/L) + p2*cos(2*pi*(-t + y)/L) + p3*sin(2*pi*(-t + z)/L))/(kappa - 1)) + (p0 + p1*sin(2*pi*(-t + x)/L) + p2*cos(2*pi*(-t + y)/L) + p3*sin(2*pi*(-t + z)/L))/(kappa - 1))*cos(2*pi*(-t + z)/L)/(L*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))**2) + 3*pi*d3*(7*d0/100 + 7*d1*sin(2*pi*(-t + x)/L)/100 + 7*d2*cos(2*pi*(-t + y)/L)/100 + 7*d3*sin(2*pi*(-t + z)/L)/100 + (kappa - 1)*(7*d0/100 + 7*d1*sin(2*pi*(-t + x)/L)/100 + 7*d2*cos(2*pi*(-t + y)/L)/100 + 7*d3*sin(2*pi*(-t + z)/L)/100 - 0.5*((d0/10 + d1*sin(2*pi*(-t + x)/L)/10 + d2*cos(2*pi*(-t + y)/L)/10 + d3*sin(2*pi*(-t + z)/L)/10)**2 + (d0/5 + d1*sin(2*pi*(-t + x)/L)/5 + d2*cos(2*pi*(-t + y)/L)/5 + d3*sin(2*pi*(-t + z)/L)/5)**2 + (3*d0/10 + 3*d1*sin(2*pi*(-t + x)/L)/10 + 3*d2*cos(2*pi*(-t + y)/L)/10 + 3*d3*sin(2*pi*(-t + z)/L)/10)**2)/(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L)) + (p0 + p1*sin(2*pi*(-t + x)/L) + p2*cos(2*pi*(-t + y)/L) + p3*sin(2*pi*(-t + z)/L))/(kappa - 1)) + (p0 + p1*sin(2*pi*(-t + x)/L) + p2*cos(2*pi*(-t + y)/L) + p3*sin(2*pi*(-t + z)/L))/(kappa - 1))*cos(2*pi*(-t + z)/L)/(5*L*(d0 + d1*sin(2*pi*(-t + x)/L) + d2*cos(2*pi*(-t + y)/L) + d3*sin(2*pi*(-t + z)/L))))

    !! F function

!!     r(DENS_VAR) = r(DENS_VAR) + &
!!         ((2*Pi*u4*uf4*y*z*Cos((2*Pi*uf4*x*y*z)/L**3))/L**3 - (2*Pi*u1*uf1*Sin((2*Pi*uf1*x)/L))/L)*(d0 + d2*Cos((2*df2*Pi*y)/L) + d4*Cos((2*df4*Pi*x*y*z)/L**3) + d1*Sin((2*df1*Pi*x)/L) + d3*Sin((2*df3*Pi*z)/L)) + &
!!         &  ((2*d1*df1*Pi*Cos((2*df1*Pi*x)/L))/L - (2*d4*df4*Pi*y*z*Sin((2*df4*Pi*x*y*z)/L**3))/L**3)*(u0 + u1*Cos((2*Pi*uf1*x)/L) + u3*Cos((2*Pi*uf3*z)/L) + u2*Sin((2*Pi*uf2*y)/L) + u4*Sin((2*Pi*uf4*x*y*z)/L**3));
!! 
!!     r(MOMX_VAR) = r(MOMX_VAR) + &
!!         & (2*p1*pf1*Pi*Cos((2*pf1*Pi*x)/L))/L + (2*p4*pf4*Pi*y*z*Cos((2*pf4*Pi*x*y*z)/L**3))/L**3 + &
!!         & 2*((2*Pi*u4*uf4*y*z*Cos((2*Pi*uf4*x*y*z)/L**3))/L**3 - (2*Pi*u1*uf1*Sin((2*Pi*uf1*x)/L))/L)*(d0 + d2*Cos((2*df2*Pi*y)/L) + d4*Cos((2*df4*Pi*x*y*z)/L**3) + d1*Sin((2*df1*Pi*x)/L) + d3*Sin((2*df3*Pi*z)/L))* &
!!         &  (u0 + u1*Cos((2*Pi*uf1*x)/L) + u3*Cos((2*Pi*uf3*z)/L) + u2*Sin((2*Pi*uf2*y)/L) + u4*Sin((2*Pi*uf4*x*y*z)/L**3)) + &
!!         & ((2*d1*df1*Pi*Cos((2*df1*Pi*x)/L))/L - (2*d4*df4*Pi*y*z*Sin((2*df4*Pi*x*y*z)/L**3))/L**3)*(u0 + u1*Cos((2*Pi*uf1*x)/L) + u3*Cos((2*Pi*uf3*z)/L) + u2*Sin((2*Pi*uf2*y)/L) + u4*Sin((2*Pi*uf4*x*y*z)/L**3))**2;
!! 
!!     r(MOMY_VAR) = r(MOMY_VAR) + &
!!         & ((2*Pi*u4*uf4*y*z*Cos((2*Pi*uf4*x*y*z)/L**3))/L**3 - (2*Pi*u1*uf1*Sin((2*Pi*uf1*x)/L))/L)*(d0 + d2*Cos((2*df2*Pi*y)/L) + d4*Cos((2*df4*Pi*x*y*z)/L**3) + d1*Sin((2*df1*Pi*x)/L) + d3*Sin((2*df3*Pi*z)/L))* &
!!         &  (v0 + v2*Cos((2*Pi*vf2*y)/L) + v4*Cos((2*Pi*vf4*x*y*z)/L**3) + v1*Sin((2*Pi*vf1*x)/L) + v3*Sin((2*Pi*vf3*z)/L)) + &
!!         & (v0 + v2*Cos((2*Pi*vf2*y)/L) + v4*Cos((2*Pi*vf4*x*y*z)/L**3) + v1*Sin((2*Pi*vf1*x)/L) + v3*Sin((2*Pi*vf3*z)/L))*((2*d1*df1*Pi*Cos((2*df1*Pi*x)/L))/L - (2*d4*df4*Pi*y*z*Sin((2*df4*Pi*x*y*z)/L**3))/L**3)* &
!!         &  (u0 + u1*Cos((2*Pi*uf1*x)/L) + u3*Cos((2*Pi*uf3*z)/L) + u2*Sin((2*Pi*uf2*y)/L) + u4*Sin((2*Pi*uf4*x*y*z)/L**3)) + &
!!         & (d0 + d2*Cos((2*df2*Pi*y)/L) + d4*Cos((2*df4*Pi*x*y*z)/L**3) + d1*Sin((2*df1*Pi*x)/L) + d3*Sin((2*df3*Pi*z)/L))* &
!!         &  (u0 + u1*Cos((2*Pi*uf1*x)/L) + u3*Cos((2*Pi*uf3*z)/L) + u2*Sin((2*Pi*uf2*y)/L) + u4*Sin((2*Pi*uf4*x*y*z)/L**3))*((2*Pi*v1*vf1*Cos((2*Pi*vf1*x)/L))/L - (2*Pi*v4*vf4*y*z*Sin((2*Pi*vf4*x*y*z)/L**3))/L**3);
!! 
!!     r(MOMZ_VAR) = r(MOMZ_VAR) + &
!!         & ((2*Pi*w4*wf4*y*z*Cos((2*Pi*wf4*x*y*z)/L**3))/L**3 - (2*Pi*w1*wf1*Sin((2*Pi*wf1*x)/L))/L)*(d0 + d2*Cos((2*df2*Pi*y)/L) + d4*Cos((2*df4*Pi*x*y*z)/L**3) + d1*Sin((2*df1*Pi*x)/L) + d3*Sin((2*df3*Pi*z)/L))* &
!!         &  (u0 + u1*Cos((2*Pi*uf1*x)/L) + u3*Cos((2*Pi*uf3*z)/L) + u2*Sin((2*Pi*uf2*y)/L) + u4*Sin((2*Pi*uf4*x*y*z)/L**3)) + &
!!         & ((2*Pi*u4*uf4*y*z*Cos((2*Pi*uf4*x*y*z)/L**3))/L**3 - (2*Pi*u1*uf1*Sin((2*Pi*uf1*x)/L))/L)*(d0 + d2*Cos((2*df2*Pi*y)/L) + d4*Cos((2*df4*Pi*x*y*z)/L**3) + d1*Sin((2*df1*Pi*x)/L) + d3*Sin((2*df3*Pi*z)/L))* &
!!         &  (w0 + w1*Cos((2*Pi*wf1*x)/L) + w3*Cos((2*Pi*wf3*z)/L) + w2*Sin((2*Pi*wf2*y)/L) + w4*Sin((2*Pi*wf4*x*y*z)/L**3)) + &
!!         & ((2*d1*df1*Pi*Cos((2*df1*Pi*x)/L))/L - (2*d4*df4*Pi*y*z*Sin((2*df4*Pi*x*y*z)/L**3))/L**3)*(u0 + u1*Cos((2*Pi*uf1*x)/L) + u3*Cos((2*Pi*uf3*z)/L) + u2*Sin((2*Pi*uf2*y)/L) + u4*Sin((2*Pi*uf4*x*y*z)/L**3))* &
!!         &  (w0 + w1*Cos((2*Pi*wf1*x)/L) + w3*Cos((2*Pi*wf3*z)/L) + w2*Sin((2*Pi*wf2*y)/L) + w4*Sin((2*Pi*wf4*x*y*z)/L**3));
!!          
!!     r(ENER_VAR) = r(ENER_VAR) + &
!!         & ((2*Pi*u4*uf4*y*z*Cos((2*Pi*uf4*x*y*z)/L**3))/L**3 - (2*Pi*u1*uf1*Sin((2*Pi*uf1*x)/L))/L)* &
!!         &    (p0 + p2*Cos((2*pf2*Pi*y)/L) + p1*Sin((2*pf1*Pi*x)/L) + p3*Sin((2*pf3*Pi*z)/L) + p4*Sin((2*pf4*Pi*x*y*z)/L**3) +  &
!!         &      (p0 + p2*Cos((2*pf2*Pi*y)/L) + p1*Sin((2*pf1*Pi*x)/L) + p3*Sin((2*pf3*Pi*z)/L) + p4*Sin((2*pf4*Pi*x*y*z)/L**3))/(-1 + kappa) +  &
!!         &      ((d0 + d2*Cos((2*df2*Pi*y)/L) + d4*Cos((2*df4*Pi*x*y*z)/L**3) + d1*Sin((2*df1*Pi*x)/L) + d3*Sin((2*df3*Pi*z)/L))*&
!!         &         ((v0 + v2*Cos((2*Pi*vf2*y)/L) + v4*Cos((2*Pi*vf4*x*y*z)/L**3) + v1*Sin((2*Pi*vf1*x)/L) + v3*Sin((2*Pi*vf3*z)/L))**2 +  &
!!         &           (u0 + u1*Cos((2*Pi*uf1*x)/L) + u3*Cos((2*Pi*uf3*z)/L) + u2*Sin((2*Pi*uf2*y)/L) + u4*Sin((2*Pi*uf4*x*y*z)/L**3))**2 +  &
!!         &           (w0 + w1*Cos((2*Pi*wf1*x)/L) + w3*Cos((2*Pi*wf3*z)/L) + w2*Sin((2*Pi*wf2*y)/L) + w4*Sin((2*Pi*wf4*x*y*z)/L**3))**2))/2.) +  &
!!         &   (u0 + u1*Cos((2*Pi*uf1*x)/L) + u3*Cos((2*Pi*uf3*z)/L) + u2*Sin((2*Pi*uf2*y)/L) + u4*Sin((2*Pi*uf4*x*y*z)/L**3))*&
!!         &    ((2*p1*pf1*Pi*Cos((2*pf1*Pi*x)/L))/L + (2*p4*pf4*Pi*y*z*Cos((2*pf4*Pi*x*y*z)/L**3))/L**3 + ((2*p1*pf1*Pi*Cos((2*pf1*Pi*x)/L))/L + (2*p4*pf4*Pi*y*z*Cos((2*pf4*Pi*x*y*z)/L**3))/L**3)/(-1 + kappa) +  &
!!         &      ((d0 + d2*Cos((2*df2*Pi*y)/L) + d4*Cos((2*df4*Pi*x*y*z)/L**3) + d1*Sin((2*df1*Pi*x)/L) + d3*Sin((2*df3*Pi*z)/L))*&
!!         &         (2*((2*Pi*u4*uf4*y*z*Cos((2*Pi*uf4*x*y*z)/L**3))/L**3 - (2*Pi*u1*uf1*Sin((2*Pi*uf1*x)/L))/L)*&
!!         &            (u0 + u1*Cos((2*Pi*uf1*x)/L) + u3*Cos((2*Pi*uf3*z)/L) + u2*Sin((2*Pi*uf2*y)/L) + u4*Sin((2*Pi*uf4*x*y*z)/L**3)) +  &
!!         &           2*(v0 + v2*Cos((2*Pi*vf2*y)/L) + v4*Cos((2*Pi*vf4*x*y*z)/L**3) + v1*Sin((2*Pi*vf1*x)/L) + v3*Sin((2*Pi*vf3*z)/L))*&
!!         &            ((2*Pi*v1*vf1*Cos((2*Pi*vf1*x)/L))/L - (2*Pi*v4*vf4*y*z*Sin((2*Pi*vf4*x*y*z)/L**3))/L**3) +  &
!!         &           2*((2*Pi*w4*wf4*y*z*Cos((2*Pi*wf4*x*y*z)/L**3))/L**3 - (2*Pi*w1*wf1*Sin((2*Pi*wf1*x)/L))/L)*&
!!         &            (w0 + w1*Cos((2*Pi*wf1*x)/L) + w3*Cos((2*Pi*wf3*z)/L) + w2*Sin((2*Pi*wf2*y)/L) + w4*Sin((2*Pi*wf4*x*y*z)/L**3))))/2. +  &
!!         &      (((2*d1*df1*Pi*Cos((2*df1*Pi*x)/L))/L - (2*d4*df4*Pi*y*z*Sin((2*df4*Pi*x*y*z)/L**3))/L**3)*&
!!         &         ((v0 + v2*Cos((2*Pi*vf2*y)/L) + v4*Cos((2*Pi*vf4*x*y*z)/L**3) + v1*Sin((2*Pi*vf1*x)/L) + v3*Sin((2*Pi*vf3*z)/L))**2 +  &
!!         &           (u0 + u1*Cos((2*Pi*uf1*x)/L) + u3*Cos((2*Pi*uf3*z)/L) + u2*Sin((2*Pi*uf2*y)/L) + u4*Sin((2*Pi*uf4*x*y*z)/L**3))**2 +  &
!!         &           (w0 + w1*Cos((2*Pi*wf1*x)/L) + w3*Cos((2*Pi*wf3*z)/L) + w2*Sin((2*Pi*wf2*y)/L) + w4*Sin((2*Pi*wf4*x*y*z)/L**3))**2))/2.);
!!          
!!     !! G function
!! 
!!     r(DENS_VAR) = r(DENS_VAR) + &
!!          & (v0 + v2*Cos((2*Pi*vf2*y)/L) + v4*Cos((2*Pi*vf4*x*y*z)/L**3) + v1*Sin((2*Pi*vf1*x)/L) + v3*Sin((2*Pi*vf3*z)/L))*((-2*d2*df2*Pi*Sin((2*df2*Pi*y)/L))/L - (2*d4*df4*Pi*x*z*Sin((2*df4*Pi*x*y*z)/L**3))/L**3) + &
!!          & (d0 + d2*Cos((2*df2*Pi*y)/L) + d4*Cos((2*df4*Pi*x*y*z)/L**3) + d1*Sin((2*df1*Pi*x)/L) + d3*Sin((2*df3*Pi*z)/L))*((-2*Pi*v2*vf2*Sin((2*Pi*vf2*y)/L))/L - (2*Pi*v4*vf4*x*z*Sin((2*Pi*vf4*x*y*z)/L**3))/L**3);
!! 
!!     r(MOMX_VAR) = r(MOMX_VAR) + &
!!         & ((2*Pi*u2*uf2*Cos((2*Pi*uf2*y)/L))/L + (2*Pi*u4*uf4*x*z*Cos((2*Pi*uf4*x*y*z)/L**3))/L**3)*(d0 + d2*Cos((2*df2*Pi*y)/L) + d4*Cos((2*df4*Pi*x*y*z)/L**3) + d1*Sin((2*df1*Pi*x)/L) + d3*Sin((2*df3*Pi*z)/L))*&
!!         &   (v0 + v2*Cos((2*Pi*vf2*y)/L) + v4*Cos((2*Pi*vf4*x*y*z)/L**3) + v1*Sin((2*Pi*vf1*x)/L) + v3*Sin((2*Pi*vf3*z)/L)) + &
!!         &  (v0 + v2*Cos((2*Pi*vf2*y)/L) + v4*Cos((2*Pi*vf4*x*y*z)/L**3) + v1*Sin((2*Pi*vf1*x)/L) + v3*Sin((2*Pi*vf3*z)/L))*((-2*d2*df2*Pi*Sin((2*df2*Pi*y)/L))/L - (2*d4*df4*Pi*x*z*Sin((2*df4*Pi*x*y*z)/L**3))/L**3)*&
!!         &   (u0 + u1*Cos((2*Pi*uf1*x)/L) + u3*Cos((2*Pi*uf3*z)/L) + u2*Sin((2*Pi*uf2*y)/L) + u4*Sin((2*Pi*uf4*x*y*z)/L**3)) + &
!!         &  (d0 + d2*Cos((2*df2*Pi*y)/L) + d4*Cos((2*df4*Pi*x*y*z)/L**3) + d1*Sin((2*df1*Pi*x)/L) + d3*Sin((2*df3*Pi*z)/L))*&
!!         &   (u0 + u1*Cos((2*Pi*uf1*x)/L) + u3*Cos((2*Pi*uf3*z)/L) + u2*Sin((2*Pi*uf2*y)/L) + u4*Sin((2*Pi*uf4*x*y*z)/L**3))*((-2*Pi*v2*vf2*Sin((2*Pi*vf2*y)/L))/L - (2*Pi*v4*vf4*x*z*Sin((2*Pi*vf4*x*y*z)/L**3))/L**3);
!!          
!!     r(MOMY_VAR) = r(MOMY_VAR) + &
!!         & (2*p4*pf4*Pi*x*z*Cos((2*pf4*Pi*x*y*z)/L**3))/L**3 - (2*p2*pf2*Pi*Sin((2*pf2*Pi*y)/L))/L + &
!!         &  (v0 + v2*Cos((2*Pi*vf2*y)/L) + v4*Cos((2*Pi*vf4*x*y*z)/L**3) + v1*Sin((2*Pi*vf1*x)/L) + v3*Sin((2*Pi*vf3*z)/L))**2*&
!!         &   ((-2*d2*df2*Pi*Sin((2*df2*Pi*y)/L))/L - (2*d4*df4*Pi*x*z*Sin((2*df4*Pi*x*y*z)/L**3))/L**3) + &
!!         &  2*(d0 + d2*Cos((2*df2*Pi*y)/L) + d4*Cos((2*df4*Pi*x*y*z)/L**3) + d1*Sin((2*df1*Pi*x)/L) + d3*Sin((2*df3*Pi*z)/L))*&
!!         &   (v0 + v2*Cos((2*Pi*vf2*y)/L) + v4*Cos((2*Pi*vf4*x*y*z)/L**3) + v1*Sin((2*Pi*vf1*x)/L) + v3*Sin((2*Pi*vf3*z)/L))*((-2*Pi*v2*vf2*Sin((2*Pi*vf2*y)/L))/L - (2*Pi*v4*vf4*x*z*Sin((2*Pi*vf4*x*y*z)/L**3))/L**3);
!!          
!!     r(MOMZ_VAR) = r(MOMZ_VAR) + &
!!         & ((2*Pi*w2*wf2*Cos((2*Pi*wf2*y)/L))/L + (2*Pi*w4*wf4*x*z*Cos((2*Pi*wf4*x*y*z)/L**3))/L**3)*(d0 + d2*Cos((2*df2*Pi*y)/L) + d4*Cos((2*df4*Pi*x*y*z)/L**3) + d1*Sin((2*df1*Pi*x)/L) + d3*Sin((2*df3*Pi*z)/L))*&
!!         &  (v0 + v2*Cos((2*Pi*vf2*y)/L) + v4*Cos((2*Pi*vf4*x*y*z)/L**3) + v1*Sin((2*Pi*vf1*x)/L) + v3*Sin((2*Pi*vf3*z)/L)) + &
!!         & (v0 + v2*Cos((2*Pi*vf2*y)/L) + v4*Cos((2*Pi*vf4*x*y*z)/L**3) + v1*Sin((2*Pi*vf1*x)/L) + v3*Sin((2*Pi*vf3*z)/L))*((-2*d2*df2*Pi*Sin((2*df2*Pi*y)/L))/L - (2*d4*df4*Pi*x*z*Sin((2*df4*Pi*x*y*z)/L**3))/L**3)*&
!!         &  (w0 + w1*Cos((2*Pi*wf1*x)/L) + w3*Cos((2*Pi*wf3*z)/L) + w2*Sin((2*Pi*wf2*y)/L) + w4*Sin((2*Pi*wf4*x*y*z)/L**3)) + &
!!         & (d0 + d2*Cos((2*df2*Pi*y)/L) + d4*Cos((2*df4*Pi*x*y*z)/L**3) + d1*Sin((2*df1*Pi*x)/L) + d3*Sin((2*df3*Pi*z)/L))*((-2*Pi*v2*vf2*Sin((2*Pi*vf2*y)/L))/L - (2*Pi*v4*vf4*x*z*Sin((2*Pi*vf4*x*y*z)/L**3))/L**3)*&
!!         &  (w0 + w1*Cos((2*Pi*wf1*x)/L) + w3*Cos((2*Pi*wf3*z)/L) + w2*Sin((2*Pi*wf2*y)/L) + w4*Sin((2*Pi*wf4*x*y*z)/L**3));
!! 
!!     r(ENER_VAR) = r(ENER_VAR) + &
!!         & ((-2*Pi*v2*vf2*Sin((2*Pi*vf2*y)/L))/L - (2*Pi*v4*vf4*x*z*Sin((2*Pi*vf4*x*y*z)/L**3))/L**3)*&
!!         &  (p0 + p2*Cos((2*pf2*Pi*y)/L) + p1*Sin((2*pf1*Pi*x)/L) + p3*Sin((2*pf3*Pi*z)/L) + p4*Sin((2*pf4*Pi*x*y*z)/L**3) + &
!!         &    (p0 + p2*Cos((2*pf2*Pi*y)/L) + p1*Sin((2*pf1*Pi*x)/L) + p3*Sin((2*pf3*Pi*z)/L) + p4*Sin((2*pf4*Pi*x*y*z)/L**3))/(-1 + kappa) + &
!!         &    ((d0 + d2*Cos((2*df2*Pi*y)/L) + d4*Cos((2*df4*Pi*x*y*z)/L**3) + d1*Sin((2*df1*Pi*x)/L) + d3*Sin((2*df3*Pi*z)/L))*&
!!         &       ((v0 + v2*Cos((2*Pi*vf2*y)/L) + v4*Cos((2*Pi*vf4*x*y*z)/L**3) + v1*Sin((2*Pi*vf1*x)/L) + v3*Sin((2*Pi*vf3*z)/L))**2 + &
!!         &         (u0 + u1*Cos((2*Pi*uf1*x)/L) + u3*Cos((2*Pi*uf3*z)/L) + u2*Sin((2*Pi*uf2*y)/L) + u4*Sin((2*Pi*uf4*x*y*z)/L**3))**2 + &
!!         &         (w0 + w1*Cos((2*Pi*wf1*x)/L) + w3*Cos((2*Pi*wf3*z)/L) + w2*Sin((2*Pi*wf2*y)/L) + w4*Sin((2*Pi*wf4*x*y*z)/L**3))**2))/2.) + &
!!         & (v0 + v2*Cos((2*Pi*vf2*y)/L) + v4*Cos((2*Pi*vf4*x*y*z)/L**3) + v1*Sin((2*Pi*vf1*x)/L) + v3*Sin((2*Pi*vf3*z)/L))*&
!!         &  ((2*p4*pf4*Pi*x*z*Cos((2*pf4*Pi*x*y*z)/L**3))/L**3 - (2*p2*pf2*Pi*Sin((2*pf2*Pi*y)/L))/L + ((2*p4*pf4*Pi*x*z*Cos((2*pf4*Pi*x*y*z)/L**3))/L**3 - (2*p2*pf2*Pi*Sin((2*pf2*Pi*y)/L))/L)/(-1 + kappa) + &
!!         &    ((d0 + d2*Cos((2*df2*Pi*y)/L) + d4*Cos((2*df4*Pi*x*y*z)/L**3) + d1*Sin((2*df1*Pi*x)/L) + d3*Sin((2*df3*Pi*z)/L))*&
!!         &       (2*((2*Pi*u2*uf2*Cos((2*Pi*uf2*y)/L))/L + (2*Pi*u4*uf4*x*z*Cos((2*Pi*uf4*x*y*z)/L**3))/L**3)*&
!!         &          (u0 + u1*Cos((2*Pi*uf1*x)/L) + u3*Cos((2*Pi*uf3*z)/L) + u2*Sin((2*Pi*uf2*y)/L) + u4*Sin((2*Pi*uf4*x*y*z)/L**3)) + &
!!         &         2*(v0 + v2*Cos((2*Pi*vf2*y)/L) + v4*Cos((2*Pi*vf4*x*y*z)/L**3) + v1*Sin((2*Pi*vf1*x)/L) + v3*Sin((2*Pi*vf3*z)/L))*&
!!         &          ((-2*Pi*v2*vf2*Sin((2*Pi*vf2*y)/L))/L - (2*Pi*v4*vf4*x*z*Sin((2*Pi*vf4*x*y*z)/L**3))/L**3) + &
!!         &         2*((2*Pi*w2*wf2*Cos((2*Pi*wf2*y)/L))/L + (2*Pi*w4*wf4*x*z*Cos((2*Pi*wf4*x*y*z)/L**3))/L**3)*&
!!         &          (w0 + w1*Cos((2*Pi*wf1*x)/L) + w3*Cos((2*Pi*wf3*z)/L) + w2*Sin((2*Pi*wf2*y)/L) + w4*Sin((2*Pi*wf4*x*y*z)/L**3))))/2. + &
!!         &    (((-2*d2*df2*Pi*Sin((2*df2*Pi*y)/L))/L - (2*d4*df4*Pi*x*z*Sin((2*df4*Pi*x*y*z)/L**3))/L**3)*&
!!         &       ((v0 + v2*Cos((2*Pi*vf2*y)/L) + v4*Cos((2*Pi*vf4*x*y*z)/L**3) + v1*Sin((2*Pi*vf1*x)/L) + v3*Sin((2*Pi*vf3*z)/L))**2 + &
!!         &         (u0 + u1*Cos((2*Pi*uf1*x)/L) + u3*Cos((2*Pi*uf3*z)/L) + u2*Sin((2*Pi*uf2*y)/L) + u4*Sin((2*Pi*uf4*x*y*z)/L**3))**2 + &
!!         &         (w0 + w1*Cos((2*Pi*wf1*x)/L) + w3*Cos((2*Pi*wf3*z)/L) + w2*Sin((2*Pi*wf2*y)/L) + w4*Sin((2*Pi*wf4*x*y*z)/L**3))**2))/2.);
!!          
!!     !! H function
!!     r(DENS_VAR) = r(DENS_VAR) + &
!!         & (d0 + d2*Cos((2*df2*Pi*y)/L) + d4*Cos((2*df4*Pi*x*y*z)/L**3) + d1*Sin((2*df1*Pi*x)/L) + d3*Sin((2*df3*Pi*z)/L))*((2*Pi*w4*wf4*x*y*Cos((2*Pi*wf4*x*y*z)/L**3))/L**3 - (2*Pi*w3*wf3*Sin((2*Pi*wf3*z)/L))/L) + &
!!         & ((2*d3*df3*Pi*Cos((2*df3*Pi*z)/L))/L - (2*d4*df4*Pi*x*y*Sin((2*df4*Pi*x*y*z)/L**3))/L**3)*(w0 + w1*Cos((2*Pi*wf1*x)/L) + w3*Cos((2*Pi*wf3*z)/L) + w2*Sin((2*Pi*wf2*y)/L) + w4*Sin((2*Pi*wf4*x*y*z)/L**3));
!!          
!!     r(MOMX_VAR) = r(MOMX_VAR) + &
!!         & (d0 + d2*Cos((2*df2*Pi*y)/L) + d4*Cos((2*df4*Pi*x*y*z)/L**3) + d1*Sin((2*df1*Pi*x)/L) + d3*Sin((2*df3*Pi*z)/L))*((2*Pi*w4*wf4*x*y*Cos((2*Pi*wf4*x*y*z)/L**3))/L**3 - (2*Pi*w3*wf3*Sin((2*Pi*wf3*z)/L))/L)*&
!!         &  (u0 + u1*Cos((2*Pi*uf1*x)/L) + u3*Cos((2*Pi*uf3*z)/L) + u2*Sin((2*Pi*uf2*y)/L) + u4*Sin((2*Pi*uf4*x*y*z)/L**3)) + &
!!         & (d0 + d2*Cos((2*df2*Pi*y)/L) + d4*Cos((2*df4*Pi*x*y*z)/L**3) + d1*Sin((2*df1*Pi*x)/L) + d3*Sin((2*df3*Pi*z)/L))*((2*Pi*u4*uf4*x*y*Cos((2*Pi*uf4*x*y*z)/L**3))/L**3 - (2*Pi*u3*uf3*Sin((2*Pi*uf3*z)/L))/L)*&
!!         &  (w0 + w1*Cos((2*Pi*wf1*x)/L) + w3*Cos((2*Pi*wf3*z)/L) + w2*Sin((2*Pi*wf2*y)/L) + w4*Sin((2*Pi*wf4*x*y*z)/L**3)) + &
!!         & ((2*d3*df3*Pi*Cos((2*df3*Pi*z)/L))/L - (2*d4*df4*Pi*x*y*Sin((2*df4*Pi*x*y*z)/L**3))/L**3)*(u0 + u1*Cos((2*Pi*uf1*x)/L) + u3*Cos((2*Pi*uf3*z)/L) + u2*Sin((2*Pi*uf2*y)/L) + u4*Sin((2*Pi*uf4*x*y*z)/L**3))*&
!!         &  (w0 + w1*Cos((2*Pi*wf1*x)/L) + w3*Cos((2*Pi*wf3*z)/L) + w2*Sin((2*Pi*wf2*y)/L) + w4*Sin((2*Pi*wf4*x*y*z)/L**3));
!!          
!!     r(MOMY_VAR) = r(MOMY_VAR) + &
!!         & (d0 + d2*Cos((2*df2*Pi*y)/L) + d4*Cos((2*df4*Pi*x*y*z)/L**3) + d1*Sin((2*df1*Pi*x)/L) + d3*Sin((2*df3*Pi*z)/L))*&
!!         &  (v0 + v2*Cos((2*Pi*vf2*y)/L) + v4*Cos((2*Pi*vf4*x*y*z)/L**3) + v1*Sin((2*Pi*vf1*x)/L) + v3*Sin((2*Pi*vf3*z)/L))*((2*Pi*w4*wf4*x*y*Cos((2*Pi*wf4*x*y*z)/L**3))/L**3 - (2*Pi*w3*wf3*Sin((2*Pi*wf3*z)/L))/L) + &
!!         & (v0 + v2*Cos((2*Pi*vf2*y)/L) + v4*Cos((2*Pi*vf4*x*y*z)/L**3) + v1*Sin((2*Pi*vf1*x)/L) + v3*Sin((2*Pi*vf3*z)/L))*((2*d3*df3*Pi*Cos((2*df3*Pi*z)/L))/L - (2*d4*df4*Pi*x*y*Sin((2*df4*Pi*x*y*z)/L**3))/L**3)*&
!!         &  (w0 + w1*Cos((2*Pi*wf1*x)/L) + w3*Cos((2*Pi*wf3*z)/L) + w2*Sin((2*Pi*wf2*y)/L) + w4*Sin((2*Pi*wf4*x*y*z)/L**3)) + &
!!         & (d0 + d2*Cos((2*df2*Pi*y)/L) + d4*Cos((2*df4*Pi*x*y*z)/L**3) + d1*Sin((2*df1*Pi*x)/L) + d3*Sin((2*df3*Pi*z)/L))*((2*Pi*v3*vf3*Cos((2*Pi*vf3*z)/L))/L - (2*Pi*v4*vf4*x*y*Sin((2*Pi*vf4*x*y*z)/L**3))/L**3)*&
!!         &  (w0 + w1*Cos((2*Pi*wf1*x)/L) + w3*Cos((2*Pi*wf3*z)/L) + w2*Sin((2*Pi*wf2*y)/L) + w4*Sin((2*Pi*wf4*x*y*z)/L**3));
!!          
!!     r(MOMZ_VAR) = r(MOMZ_VAR) + &
!!         & (2*p3*pf3*Pi*Cos((2*pf3*Pi*z)/L))/L + (2*p4*pf4*Pi*x*y*Cos((2*pf4*Pi*x*y*z)/L**3))/L**3 + &
!!         &  2*(d0 + d2*Cos((2*df2*Pi*y)/L) + d4*Cos((2*df4*Pi*x*y*z)/L**3) + d1*Sin((2*df1*Pi*x)/L) + d3*Sin((2*df3*Pi*z)/L))*((2*Pi*w4*wf4*x*y*Cos((2*Pi*wf4*x*y*z)/L**3))/L**3 - (2*Pi*w3*wf3*Sin((2*Pi*wf3*z)/L))/L)*&
!!         &   (w0 + w1*Cos((2*Pi*wf1*x)/L) + w3*Cos((2*Pi*wf3*z)/L) + w2*Sin((2*Pi*wf2*y)/L) + w4*Sin((2*Pi*wf4*x*y*z)/L**3)) + &
!!         &  ((2*d3*df3*Pi*Cos((2*df3*Pi*z)/L))/L - (2*d4*df4*Pi*x*y*Sin((2*df4*Pi*x*y*z)/L**3))/L**3)*(w0 + w1*Cos((2*Pi*wf1*x)/L) + w3*Cos((2*Pi*wf3*z)/L) + w2*Sin((2*Pi*wf2*y)/L) + w4*Sin((2*Pi*wf4*x*y*z)/L**3))**2;
!!          
!!     r(ENER_VAR) = r(ENER_VAR) + &
!!         & ((2*Pi*w4*wf4*x*y*Cos((2*Pi*wf4*x*y*z)/L**3))/L**3 - (2*Pi*w3*wf3*Sin((2*Pi*wf3*z)/L))/L)*&
!!         &  (p0 + p2*Cos((2*pf2*Pi*y)/L) + p1*Sin((2*pf1*Pi*x)/L) + p3*Sin((2*pf3*Pi*z)/L) + p4*Sin((2*pf4*Pi*x*y*z)/L**3) + &
!!         &    (p0 + p2*Cos((2*pf2*Pi*y)/L) + p1*Sin((2*pf1*Pi*x)/L) + p3*Sin((2*pf3*Pi*z)/L) + p4*Sin((2*pf4*Pi*x*y*z)/L**3))/(-1 + kappa) + &
!!         &    ((d0 + d2*Cos((2*df2*Pi*y)/L) + d4*Cos((2*df4*Pi*x*y*z)/L**3) + d1*Sin((2*df1*Pi*x)/L) + d3*Sin((2*df3*Pi*z)/L))*&
!!         &       ((v0 + v2*Cos((2*Pi*vf2*y)/L) + v4*Cos((2*Pi*vf4*x*y*z)/L**3) + v1*Sin((2*Pi*vf1*x)/L) + v3*Sin((2*Pi*vf3*z)/L))**2 + &
!!         &         (u0 + u1*Cos((2*Pi*uf1*x)/L) + u3*Cos((2*Pi*uf3*z)/L) + u2*Sin((2*Pi*uf2*y)/L) + u4*Sin((2*Pi*uf4*x*y*z)/L**3))**2 + &
!!         &         (w0 + w1*Cos((2*Pi*wf1*x)/L) + w3*Cos((2*Pi*wf3*z)/L) + w2*Sin((2*Pi*wf2*y)/L) + w4*Sin((2*Pi*wf4*x*y*z)/L**3))**2))/2.) + &
!!         & (w0 + w1*Cos((2*Pi*wf1*x)/L) + w3*Cos((2*Pi*wf3*z)/L) + w2*Sin((2*Pi*wf2*y)/L) + w4*Sin((2*Pi*wf4*x*y*z)/L**3))*&
!!         &  ((2*p3*pf3*Pi*Cos((2*pf3*Pi*z)/L))/L + (2*p4*pf4*Pi*x*y*Cos((2*pf4*Pi*x*y*z)/L**3))/L**3 + ((2*p3*pf3*Pi*Cos((2*pf3*Pi*z)/L))/L + (2*p4*pf4*Pi*x*y*Cos((2*pf4*Pi*x*y*z)/L**3))/L**3)/(-1 + kappa) + &
!!         &    ((d0 + d2*Cos((2*df2*Pi*y)/L) + d4*Cos((2*df4*Pi*x*y*z)/L**3) + d1*Sin((2*df1*Pi*x)/L) + d3*Sin((2*df3*Pi*z)/L))*&
!!         &       (2*((2*Pi*u4*uf4*x*y*Cos((2*Pi*uf4*x*y*z)/L**3))/L**3 - (2*Pi*u3*uf3*Sin((2*Pi*uf3*z)/L))/L)*&
!!         &          (u0 + u1*Cos((2*Pi*uf1*x)/L) + u3*Cos((2*Pi*uf3*z)/L) + u2*Sin((2*Pi*uf2*y)/L) + u4*Sin((2*Pi*uf4*x*y*z)/L**3)) + &
!!         &         2*(v0 + v2*Cos((2*Pi*vf2*y)/L) + v4*Cos((2*Pi*vf4*x*y*z)/L**3) + v1*Sin((2*Pi*vf1*x)/L) + v3*Sin((2*Pi*vf3*z)/L))*&
!!         &          ((2*Pi*v3*vf3*Cos((2*Pi*vf3*z)/L))/L - (2*Pi*v4*vf4*x*y*Sin((2*Pi*vf4*x*y*z)/L**3))/L**3) + &
!!         &         2*((2*Pi*w4*wf4*x*y*Cos((2*Pi*wf4*x*y*z)/L**3))/L**3 - (2*Pi*w3*wf3*Sin((2*Pi*wf3*z)/L))/L)*&
!!         &          (w0 + w1*Cos((2*Pi*wf1*x)/L) + w3*Cos((2*Pi*wf3*z)/L) + w2*Sin((2*Pi*wf2*y)/L) + w4*Sin((2*Pi*wf4*x*y*z)/L**3))))/2. + &
!!         &    (((2*d3*df3*Pi*Cos((2*df3*Pi*z)/L))/L - (2*d4*df4*Pi*x*y*Sin((2*df4*Pi*x*y*z)/L**3))/L**3)*&
!!         &       ((v0 + v2*Cos((2*Pi*vf2*y)/L) + v4*Cos((2*Pi*vf4*x*y*z)/L**3) + v1*Sin((2*Pi*vf1*x)/L) + v3*Sin((2*Pi*vf3*z)/L))**2 + &
!!         &         (u0 + u1*Cos((2*Pi*uf1*x)/L) + u3*Cos((2*Pi*uf3*z)/L) + u2*Sin((2*Pi*uf2*y)/L) + u4*Sin((2*Pi*uf4*x*y*z)/L**3))**2 + &
!!         &         (w0 + w1*Cos((2*Pi*wf1*x)/L) + w3*Cos((2*Pi*wf3*z)/L) + w2*Sin((2*Pi*wf2*y)/L) + w4*Sin((2*Pi*wf4*x*y*z)/L**3))**2))/2.);

end function

end module
