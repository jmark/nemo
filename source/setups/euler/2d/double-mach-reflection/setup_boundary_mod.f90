module setup_boundary_mod

use constants_mod

contains

subroutine setup_side_boundary(side,inner,outer)

    use setup_config_mod
    use equations_const_mod
    use mesh_types_mod, only: side_t
    use quad_utils_mod, only: side_get_coords
    use equations_mod, only: pressure

    type(side_t), intent(in)    :: side
    real(dp), intent(in)        :: inner(N_NODES,N_VARS)
    real(dp), intent(out)       :: outer(N_NODES,N_VARS)

    real(dp) :: coords(N_NODES,2), center(2)

    real(dp) :: dens(N_NODES)
    real(dp) :: velx(N_NODES)
    real(dp) :: vely(N_NODES)
    real(dp) :: pres(N_NODES)

    if (side%boundary == NOR) then
        outer(:,DENS_VAR) = dens_L
        outer(:,MOMX_VAR) = momx_L
        outer(:,MOMY_VAR) = momy_L
        outer(:,ENER_VAR) = ener_L

    else if (side%boundary == SOU) then
        outer(:,DENS_VAR) = dens_R
        outer(:,MOMX_VAR) = momx_R
        outer(:,MOMY_VAR) = momy_R
        outer(:,ENER_VAR) = ener_R

    else
        outer = inner

    end if

end subroutine

end module
