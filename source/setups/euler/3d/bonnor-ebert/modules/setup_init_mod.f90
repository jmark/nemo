module setup_init_mod

use constants_mod

real(dp), allocatable, save :: profile(:,:)

contains

subroutine fill(quad)

    use setup_config_mod
    use equations_const_mod

    use mesh_types_mod, only: quad_t
    use quad_utils_mod, only: quad_get_coords

    type(quad_t), intent(inout) :: quad

    real(dp) :: dens(N_NODES,N_NODES,N_NODES)
    real(dp) :: velx(N_NODES,N_NODES,N_NODES)
    real(dp) :: vely(N_NODES,N_NODES,N_NODES)
    real(dp) :: velz(N_NODES,N_NODES,N_NODES)
    real(dp) :: momx(N_NODES,N_NODES,N_NODES)
    real(dp) :: momy(N_NODES,N_NODES,N_NODES)
    real(dp) :: momz(N_NODES,N_NODES,N_NODES)
    real(dp) :: pres(N_NODES,N_NODES,N_NODES)
    real(dp) :: ener(N_NODES,N_NODES,N_NODES)

    real(dp) :: rpot(N_NODES,N_NODES,N_NODES)
    real(dp) :: racc(N_NODES,N_NODES,N_NODES)

    real(dp) :: coords(N_NODES,N_NODES,N_NODES,N_DIMS)
    real(dp) :: radius(N_NODES,N_NODES,N_NODES)

    integer :: i,j,k,l,r,m

    !write (*,'(i4,99(ES12.4))') quad%locid,quad%delta,quad%center

    coords = quad_get_coords(quad)
    radius = SQRT(SUM(coords**2,dim=N_DIMS+1))

    ! write (*,*) coords(1,1,1,:)

    dens = dens0
    velx = velx0
    vely = vely0
    velz = velz0
    pres = pres0

# if 1
    do i = 1,N_NODES
        do j = 1,N_NODES
            do k = 1,N_NODES
                L = 1
                R = size(profile,dim=2)

                do while (L <= R)
                    m = floor(0.5*(L + R))
                    if (profile(1,m) < radius(i,j,k)) then
                        L = m + 1
                    else if (profile(1,m) > radius(i,j,k)) then
                        R = m - 1
                    else
                        exit
                    end if
                end do

                dens(i,j,k) = profile(2,m)
                pres(i,j,k) = profile(3,m)
                rpot(i,j,k) = profile(4,m)
                racc(i,j,k) = profile(5,m)
            end do
        end do
    end do
# endif

    momx = dens * velx
    momy = dens * vely
    momz = dens * velz

    ener = pres/(kappa-1) + 0.5*dens*(velx**2+vely**2+velz**2)

    ! quad%pld%hydro%state(:,:,:,DENS_VAR) = 0.5 + coords(:,:,:,3) ! + coords(:,:,:,2) + coords(:,:,:,3)

    quad%pld%hydro%state(:,:,:,DENS_VAR) = dens
    quad%pld%hydro%state(:,:,:,MOMX_VAR) = momx
    quad%pld%hydro%state(:,:,:,MOMY_VAR) = momy
    quad%pld%hydro%state(:,:,:,MOMZ_VAR) = momz
    quad%pld%hydro%state(:,:,:,ENER_VAR) = ener

    quad%aux%gravity%rpot = rpot
    quad%aux%gravity%racc(:,:,:,1) = racc * coords(:,:,:,1) / radius
    quad%aux%gravity%racc(:,:,:,2) = racc * coords(:,:,:,2) / radius
    quad%aux%gravity%racc(:,:,:,3) = racc * coords(:,:,:,3) / radius

end subroutine

function setup_init_probe_refine(quad) result(ok)

    use setup_config_mod
    use mesh_types_mod, only: quad_t
    use mesh_const_mod, only: MESH_CORNERS
    use quad_utils_mod, only: quad_get_center
    use quad_utils_mod, only: quad_get_vertices

    type(quad_t), intent(in)        :: quad
    logical                         :: ok

    real(dp) :: verts(N_DIMS,MESH_CORNERS)
    real(dp) :: center(N_DIMS)

    verts = quad_get_vertices(quad)

    center = quad_get_center(quad)

    ok = .false.
    !if (rand() < 0.3) then
    if (any(abs(SQRT(SUM((verts)**2,dim=1))) < 0.5 * 0.5*diameter)) then
    !if (any(abs(SQRT(SUM((verts)**2,dim=1))) < 0.25 * 0.5*diameter)) then
    !if (center(1) + center(2) + center(3) < 0) then
    !if (all(center < 0)) then
        ok = .true.
    end if

end function

function setup_init_probe_coarsen(quads) result(ok)

    use mesh_types_mod, only: quad_t
    use mesh_const_mod, only: MESH_CHILDREN
    use setup_config_mod, only: setup_amr_threshold_coarsen

    type(quad_t), intent(in)    :: quads(MESH_CHILDREN)
    logical                     :: ok

    real(dp)                    :: sensor(MESH_CHILDREN)

    integer :: i

    ok = .false.

end function

subroutine hook_init_810

    use share_mod, only: mesh

    use mesh_utils_mod, only: mesh_iterate_quads
    use mesh_utils_mod, only: mesh_refine

    integer :: level, i, nlines,io
    
    nlines = 0
    OPEN (55, file = '../bes.dat', status = 'old')
    DO
        READ (55,*, iostat=io)
        IF (io /=0 ) EXIT
        nlines = nlines + 1
    END DO
    CLOSE (55)

    write (*,*) '# lines', nlines

    allocate(profile(5,nlines))

    open(unit = 55, file = '../bes.dat', status = 'old')

    do i = 1,nlines
        read(55,*) profile(:,i)
    end do

    close (55)

    ! do i = 1,nlines
    !     write (*,*) i,profile(:,i)
    ! end do

    do level = mesh%minlevel,mesh%maxlevel-1
        call mesh_refine(mesh,setup_init_probe_refine)
    end do
    call mesh_iterate_quads(mesh,fill)

    deallocate(profile)

end subroutine

end module
