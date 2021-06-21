module setup_amr_mod

use constants_mod

real(dp), save :: r_blast
real(dp), save :: r_inner
real(dp), save :: r_outer

contains

subroutine hook_pre_timestep_300

    use setup_config_mod
    use share_mod, only: mesh
    use share_mod, only: rt => runtime

    use mesh_utils_mod, only: mesh_adapt
    use mesh_utils_mod, only: mesh_partition

    integer(kind=8) :: oldsize,newsize
    integer(kind=8) :: goldsize,gnewsize
 
    ! R_contact = lambda t: R_blast(t)/1.181
    ! R_reverse = lambda t: R_contact(t)*0.935

    r_blast = 1.06_dp * (blast_energy**2/ejecta_mass/dens0)**(1.0_dp/7.0_dp) * rt%simtime**(4.0_dp/7.0_dp)
    r_outer = max(0.2_dp*diameter,1.15_dp*r_blast)

    if (rt%simtime < switch_time) then
        r_inner = -99.0_dp
    else
        r_inner = 0.7_dp*r_blast
    end if

    if (mod(rt%timesteps,amr_each_ntimesteps) == 0) then
        oldsize = mesh%nquads
        goldsize = mesh%nghosts

        call mesh_adapt(mesh,setup_amr_probe_refine,setup_amr_replace_refine,setup_amr_probe_coarsen,setup_amr_replace_coarsen)
        call mesh_partition(mesh)

        newsize = mesh%nquads
        gnewsize = mesh%nghosts

        ! if (rt%mpi%root .and. mod(rt%timesteps,128) == 0) then
        !     write (*,'(a10,1(I8),(ES12.4),9(I8))') 'AMR', rt%timesteps,r_blast, rt%mpi%rank,oldsize,newsize,goldsize,gnewsize
        ! end if
            
        if (mod(rt%timesteps,128) == 0) then
            write (*,'(a10,99(I8))') 'AMR', rt%timesteps,rt%mpi%rank,oldsize,newsize,goldsize,gnewsize
        end if
    end if

end subroutine

function setup_amr_probe_refine(quad) result(ok)

    use setup_config_mod
    use equations_const_mod
    use mesh_types_mod, only: quad_t
    use quad_utils_mod, only: quad_get_center

    type(quad_t), intent(in)    :: quad
    logical                     :: ok

    real(dp) :: center(N_DIMS), r_cell

    center = quad_get_center(quad)
    r_cell = SQRT(SUM(center**2))

    if (r_inner < r_cell .and. r_cell < r_outer) then
        ok = .true.
    else
        ok = .false.
    end if

end function

function setup_amr_probe_coarsen(family) result(ok)

    use setup_config_mod
    use mesh_types_mod, only: quad_t
    use mesh_const_mod, only: MESH_CHILDREN
    use quad_utils_mod, only: quad_get_center

    type(quad_t), intent(in)    :: family(MESH_CHILDREN)
    logical                     :: ok

    real(dp) :: center(N_DIMS), r_cell

    center = quad_get_center(family(MESH_CHILDREN))
    r_cell = SQRT(SUM(center**2))

    if (r_cell < r_inner) then
        ok = .true.
    else
        ok = .false.
    end if

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

    !do q = DGFV_QMAX,0,-1
    do q = 0,0,-1
        do f = 1,MESH_CHILDREN
            do h = 1,N_VARS
                outgoing(f)%pld%hydro%state(:,:,:,h) = &
                    blockmulXYZ(refineMat(:,:,1,f,q),refineMat(:,:,2,f,q),refineMat(:,:,3,f,q),&
                        incoming%pld%hydro%state(:,:,:,h))
            end do
            ! ok = isvalid(outgoing(f)%pld%hydro%state)
            ! if (.not.ok) exit
        end do
        ! if (ok) exit
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
    do q = 0,0,-1
        do f = 1,MESH_CHILDREN
            do h = 1,N_VARS
                outgoing%pld%hydro%state(:,:,:,h) = outgoing%pld%hydro%state(:,:,:,h) + &
                    blockmulXYZ(coarseMat(:,:,1,f,q),coarseMat(:,:,2,f,q),coarseMat(:,:,3,f,q),&
                        incoming(f)%pld%hydro%state(:,:,:,h))
            end do
        end do
        ! if (isvalid(outgoing%pld%hydro%state)) exit
    end do

end subroutine

end module
