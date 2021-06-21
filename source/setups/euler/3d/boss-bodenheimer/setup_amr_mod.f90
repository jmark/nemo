module setup_amr_mod

use constants_mod

contains

subroutine hook_pre_timestep_300

    use setup_config_mod
    use share_mod, only: mesh
    use share_mod, only: rt => runtime

    use mesh_utils_mod, only: mesh_adapt
    use mesh_utils_mod, only: mesh_partition

    if (mod(rt%timesteps,amr_each_ntimesteps) == 0) then
        call mesh_adapt(mesh,setup_amr_probe_refine,setup_amr_replace_refine,setup_amr_probe_coarsen,setup_amr_replace_coarsen)
        call mesh_partition(mesh)
    end if

end subroutine

function setup_amr_probe_refine(quad) result(ok)

    use setup_config_mod
    use equations_const_mod
    use mesh_types_mod, only: quad_t

    use equations_const_mod
    use equations_mod, only: pressure
    use quad_utils_mod, only: quad_get_delta
    use quad_utils_mod, only: quad_get_center

    type(quad_t), intent(in)    :: quad
    logical                     :: ok

    real(dp) :: lambda2,delta2

    ! ok = .false.
    ! if (any(quad%pld%hydro%state(:,:,:,DENS_VAR) > 2.5*dens_cloud)) then
    !     ok = .true.
    ! end if

    ! if (SUM(quad_get_center(quad)**2) > (1.25_dp*cloudradius)**2) then
    !     ok = .false.
    !     return
    ! end if

    delta2 = maxval(quad_get_delta(quad))**2 * (1.0_dp/REAL(N_NODES**2,dp) * ncells_refine**2)
    lambda2 = PI/grav*minval(pressure(quad%pld%hydro%state)/quad%pld%hydro%state(:,:,:,DENS_VAR)**2)

    ok = lambda2 < delta2

end function

function setup_amr_probe_coarsen(family) result(ok)

    use setup_config_mod
    use mesh_types_mod, only: quad_t
    use mesh_const_mod, only: MESH_CHILDREN

    use equations_const_mod
    use equations_mod, only: pressure
    use quad_utils_mod, only: quad_get_delta

    type(quad_t), intent(in)    :: family(MESH_CHILDREN)
    logical                     :: ok

    real(dp) :: lambda2,delta2
    integer  :: f

    ! ok = .true.
    ! do f = 1,MESH_CHILDREN
    !     if (.not.all(family(f)%pld%hydro%state(:,:,:,DENS_VAR) < dens_cloud)) then
    !         ok = .false.
    !         return
    !     end if
    ! end do

    ! ok = .false.
    ! return

    delta2 = maxval(2.0_dp*quad_get_delta(family(1)))**2 * (1.0_dp/REAL(N_NODES**2,dp) * ncells_coarse**2)

    ok = .true.
    do f = 1,MESH_CHILDREN
        lambda2 = PI/grav*minval(pressure(family(f)%pld%hydro%state)/family(f)%pld%hydro%state(:,:,:,DENS_VAR)**2)
        ok = lambda2 > delta2
        if (.not.ok) exit
    end do

end function

subroutine setup_amr_replace_refine(incoming,outgoing)

    use mesh_types_mod, only: quad_t
    use mesh_const_mod, only: MESH_CHILDREN
    use kernel_share_mod, only: refineMat
    use kernel_const_mod, only: DGFV_QMAX
    use tensor_mod, only: blockmulXYZ
    use equations_mod, only: isvalid

    type(quad_t), intent(in)    :: incoming
    type(quad_t), intent(inout) :: outgoing(MESH_CHILDREN)

    integer :: q,h,f
    logical :: ok

    do q = DGFV_QMAX-2,0,-1
    !do q = DGFV_QMAX,0,-1
        do f = 1,MESH_CHILDREN
            do h = 1,N_VARS
                outgoing(f)%pld%hydro%state(:,:,:,h) = &
                    blockmulXYZ(refineMat(:,:,X_DIR,f,q),refineMat(:,:,Y_DIR,f,q),refineMat(:,:,Z_DIR,f,q),&
                        incoming%pld%hydro%state(:,:,:,h))
            end do
            ok = isvalid(outgoing(f)%pld%hydro%state)
            if (.not.ok) exit
        end do
        if (ok) exit
    end do

end subroutine

subroutine setup_amr_replace_coarsen(incoming,outgoing)

    use mesh_types_mod, only: quad_t
    use mesh_const_mod, only: MESH_CHILDREN
    use kernel_share_mod, only: coarseMat
    use kernel_const_mod, only: DGFV_QMAX
    use tensor_mod, only: blockmulXYZ
    use equations_mod, only: isvalid

    type(quad_t), intent(in)    :: incoming(MESH_CHILDREN)
    type(quad_t), intent(inout) :: outgoing

    integer :: q,h,f

    outgoing%pld%hydro%state = 0.0_dp
    !do q = DGFV_QMAX,0,-1
    do q = DGFV_QMAX-2,0,-1
        do f = 1,MESH_CHILDREN
            do h = 1,N_VARS
                outgoing%pld%hydro%state(:,:,:,h) = outgoing%pld%hydro%state(:,:,:,h) + &
                    blockmulXYZ(coarseMat(:,:,X_DIR,f,q),coarseMat(:,:,Y_DIR,f,q),coarseMat(:,:,Z_DIR,f,q),&
                        incoming(f)%pld%hydro%state(:,:,:,h))
            end do
        end do
        if (isvalid(outgoing%pld%hydro%state)) exit
    end do

end subroutine

end module
