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

    type(quad_t), intent(in)    :: quad
    logical                     :: ok

    if (any(quad%pld%hydro%state(:,:,:,DENS_VAR) > dens_refine_threshold)) then
        ok = .true.
    else
        ok = .false.
    end if

end function

function setup_amr_probe_coarsen(family) result(ok)

    use setup_config_mod
    use equations_const_mod
    use mesh_types_mod, only: quadptr_t
    use mesh_types_mod, only: quad_t
    use mesh_const_mod, only: MESH_CHILDREN

    !type(quadptr_t), intent(in) :: family(MESH_CHILDREN)
    type(quad_t), intent(in) :: family(MESH_CHILDREN)
    logical                     :: ok

    integer :: i

    ok = .true.

    do i = 1,MESH_CHILDREN
        ! write (*,'((a12),(99(i7)))') 'probe_coarsen', i, family(i)%p%locid,family(i)%p%mesh%nquads
        !if (any(family(i)%p%pld%hydro%state(:,:,:,DENS_VAR) > dens_coarse_threshold)) then
        if (any(family(i)%pld%hydro%state(:,:,:,DENS_VAR) > dens_coarse_threshold)) then
            ok = .false.
            return
        end if
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
                    blockmulXYZ(refineMat(:,:,1,f,q),refineMat(:,:,2,f,q),refineMat(:,:,3,f,q),&
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
                    blockmulXYZ(coarseMat(:,:,1,f,q),coarseMat(:,:,2,f,q),coarseMat(:,:,3,f,q),&
                        incoming(f)%pld%hydro%state(:,:,:,h))
            end do
        end do
        if (isvalid(outgoing%pld%hydro%state)) exit
    end do

end subroutine

end module
