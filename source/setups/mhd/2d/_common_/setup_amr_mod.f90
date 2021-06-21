module setup_amr_mod

use constants_mod

real(dp), parameter :: loehner_refine = 0.8_dp
real(dp), parameter :: loehner_coarse = 0.2_dp

contains

subroutine hook_pre_timestep_300

    use setup_config_mod, only: amr_each_ntimesteps
    use share_mod, only: rt => runtime

    if (mod(rt%timesteps,amr_each_ntimesteps) == 0) then
        call setup_amr_adapt()
    end if

end subroutine

subroutine setup_amr_adapt

    use share_mod, only: mesh
    use mesh_utils_mod, only: mesh_adapt
    use mesh_utils_mod, only: mesh_partition

    return

    call mesh_adapt(mesh,setup_amr_probe_refine,setup_amr_replace_refine,setup_amr_probe_coarsen,setup_amr_replace_coarsen)
    call mesh_partition(mesh)

end subroutine

function setup_amr_probe_refine(quad) result(ok)

    use setup_config_mod
    use equations_const_mod
    use mesh_types_mod, only: quad_t
    use loehner_mod, only: calc_loehner
    use equations_mod, only: pressure

    type(quad_t), intent(in)    :: quad
    logical                     :: ok

    ok = calc_loehner(quad%pld%hydro%state(:,:,DENS_VAR)) > loehner_refine

    if (.not.ok) then
        ok = calc_loehner(pressure(quad%pld%hydro%state)) > loehner_refine
    end if

end function

function setup_amr_probe_coarsen(family) result(ok)

    use setup_config_mod
    use equations_const_mod
    use mesh_types_mod, only: quad_t
    use mesh_const_mod, only: MESH_CHILDREN
    use kernel_amr_mod, only: kernel_amr_coarsen
    use loehner_mod, only: calc_loehner
    use equations_mod, only: pressure

    type(quad_t), intent(in)    :: family(MESH_CHILDREN)
    logical                     :: ok

    integer :: f

    real(dp) :: temp1(N_NODES,N_NODES,MESH_CHILDREN)
    real(dp) :: temp2(N_NODES,N_NODES)

    do f = 1,MESH_CHILDREN
        temp1(:,:,f) = family(f)%pld%hydro%state(:,:,DENS_VAR)
    end do

    call kernel_amr_coarsen(temp1,temp2)

    ok = calc_loehner(temp2) < loehner_coarse

    if (.not.ok) then
        ok = calc_loehner(pressure(temp2)) < loehner_coarse
    end if

end function

subroutine setup_amr_replace_refine(incoming,outgoing)

    use mesh_types_mod, only: quad_t
    use mesh_const_mod, only: MESH_CHILDREN

    use kernel_amr_mod, only: kernel_amr_refine

    type(quad_t), intent(in)    :: incoming
    type(quad_t), intent(inout) :: outgoing(MESH_CHILDREN)

    integer :: f

    real(dp) :: temp(N_NODES,N_NODES,N_VARS,MESH_CHILDREN)

    call kernel_amr_refine(incoming%pld%hydro%state,temp)

    do f = 1,MESH_CHILDREN
        outgoing(f)%pld%hydro%state = temp(:,:,:,f)
    end do

end subroutine

subroutine setup_amr_replace_coarsen(incoming,outgoing)

    use mesh_types_mod, only: quad_t
    use mesh_const_mod, only: MESH_CHILDREN

    use kernel_amr_mod, only: kernel_amr_coarsen

    type(quad_t), intent(in)    :: incoming(MESH_CHILDREN)
    type(quad_t), intent(inout) :: outgoing

    integer :: f

    real(dp) :: temp(N_NODES,N_NODES,N_VARS,MESH_CHILDREN)

    do f = 1,MESH_CHILDREN
        temp(:,:,:,f) = incoming(f)%pld%hydro%state
    end do

    call kernel_amr_coarsen(temp,outgoing%pld%hydro%state)

end subroutine

end module
