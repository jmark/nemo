module kernel_evolve_mod

use constants_mod

contains

subroutine kernel_evolve(quad,flux,state)

    use mesh_types_mod, only: quad_t
    use kernel_flux_mod, only: kernel_calc_flux
    use kernel_validate_mod, only: kernel_repair

    use timedisc_probe_mod, only: timedisc_probe

# if MODULE_SOURCETERM
    use setup_sourceterm_mod, only: calc_sourceterm
# endif

    type(quad_t), intent(inout) :: quad
    real(dp), intent(out)       :: flux(N_NODES,N_NODES,N_NODES,N_VARS)
    real(dp), intent(out)       :: state(N_NODES,N_NODES,N_NODES,N_VARS)

    logical :: ok

# if MODULE_SOURCETERM
    real(dp) :: source(N_NODES,N_NODES,N_NODES,N_VARS)
# endif

    call kernel_calc_flux(quad,flux)

# if MODULE_SOURCETERM
    call calc_sourceterm(quad,source)
    flux = flux + source
# endif

    call timedisc_probe(quad,flux,state)
    ok = kernel_repair(quad,state)

end subroutine

end module
