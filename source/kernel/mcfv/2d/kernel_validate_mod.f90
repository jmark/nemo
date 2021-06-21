module kernel_validate_mod

use constants_mod

private
public :: kernel_repair
public :: kernel_validate

interface kernel_repair
    procedure kernel_repair_quad_only
    procedure kernel_repair_quad_and_state
end interface

contains

function kernel_repair_quad_only(quad) result(ok)

    use mesh_types_mod, only: quad_t

    type(quad_t), intent(inout)    :: quad
    logical                         :: ok

    ok = kernel_repair_quad_and_state(quad,quad%pld%hydro%state)

end function

function kernel_repair_quad_and_state(quad,state) result(ok)

    use share_mod, only: rt => runtime
    use share_mod, only: tr => traceback

    use mesh_types_mod, only: quad_t
    use equations_mod, only: pressure
    use mesh_const_mod, only: MESH_FACES
    use utils_mod, only: utils_get_newunit
    use quad_utils_mod, only: quad_get_center

    use status_mod, only: raise_error

# ifdef _OPENMP
    USE OMP_LIB, only: omp_get_thread_num
    USE OMP_LIB, only: omp_get_num_threads
# endif

    type(quad_t), intent(inout) :: quad
    real(dp), intent(inout)     :: state(N_NODES,N_NODES,N_VARS)
    logical                     :: ok

    real(dp) :: sfac
    real(dp) :: weight
    real(dp) :: avgs(N_VARS)
    real(dp) :: center(N_DIMS)
    real(dp) :: tiles(2,2,N_VARS)
    real(dp) :: pres(N_NODES,N_NODES)
    real(dp) :: fpres(N_NODES)

    integer :: ioUnit
    integer :: openStat
    integer :: i,j,k,h,ntrial, f
    integer :: q, p, qq, ii ,jj, mm, nn

    character(len=32) :: pvarStr
    character(len=32) :: outputfp
    character(len=32) :: countStr
    character(len=32) :: ompStr
    character(len=32) :: mpiStr

    ok = dgfv_repair_2d(quad,state)

    if (.not.ok) then
        tr%ok = .false.
        tr%errcode = 1
        tr%message = 'PP limiter failed!'

        center = quad_get_center(quad)

# ifdef _OPENMP
        write(ompStr,'(I0.6)') omp_get_thread_num()
# else
        write(ompStr,'(I0.6)') 0
# endif

        write(mpiStr,'(I0.6)') rt%mpi%rank

        outputfp = 'PP_LIMITER_' // trim(mpiStr) // '_' // trim(ompStr) // '.log'

        OPEN(UNIT     = utils_get_newunit(ioUnit) , &
             FILE     = TRIM(outputfp)     , &
             FORM     = 'FORMATTED'        , &
             STATUS   = 'UNKNOWN'          , &
             RECL     = 50000              , &
             IOSTAT   = openStat             )

        IF (openStat.NE.0) THEN
            write (*,*) 'ERROR: cannot open '// TRIM(outputfp)

        else
            weight = 1./real(N_NODES**2,dp)
            do j = 1,N_NODES; do i = 1,N_NODES
                pres(i,j) = pressure(state(i,j,:))
            end do; end do

            write (ioUnit,'(a12,1x,2(ES12.4))') 'quad center', center
            write (ioUnit,*)

            write (ioUnit,*) 'density', weight*sum(state(:,:,1))
            do i = 1,N_NODES
                write (ioUnit,101) state(i,:,1)
            end do
            write (ioUnit,*)

            write (ioUnit,*) 'momx', weight*sum(state(:,:,2))
            do i = 1,N_NODES
                write (ioUnit,101) state(i,:,2)
            end do
            write (ioUnit,*)

            write (ioUnit,*) 'momy', weight*sum(state(:,:,3))
            do i = 1,N_NODES
                write (ioUnit,101) state(i,:,3)
            end do
            write (ioUnit,*)

            write (ioUnit,*) 'energy', weight*sum(state(:,:,4))
            do i = 1,N_NODES
                write (ioUnit,101) state(i,:,4)
            end do
            write (ioUnit,*)

            write (ioUnit,*) 'pressure', weight*sum(pres)
            do i = 1,N_NODES
                write (ioUnit,101) pres(i,:)
            end do
            write (ioUnit,*)
        
            !! do f = 1,MESH_FACES
            !!     forall (i = 1:N_NODES)
            !!         fpres(i) = pressure(quad%aux%hydro%facevars(i,:,f))
            !!     end forall

            !!     write (ioUnit,*) 'facevalues', f
            !!     do h = 1,N_VARS
            !!         write (ioUnit,101) quad%aux%hydro%facevars(:,h,f)
            !!     end do
            !!     write (ioUnit,101) fpres
            !!     write (ioUnit,*)
            !! end do

            do f = 1,MESH_FACES
                write (ioUnit,*) 'faceflux', f
                do h = 1,N_VARS
                    write (ioUnit,101) quad%aux%hydro%faceflux(:,h,f)
                end do
                write (ioUnit,*)
            end do

            CLOSE(ioUnit)
        end if

        call raise_error('Kernel repair failed!')
    end if

101 format(32(ES12.4))
end function

function dgfv_repair_2d(quad,state) result(ok)

    use equations_mod, only: isvalid

    use mesh_types_mod, only: quad_t

    type(quad_t)            :: quad
    real(dp), intent(inout) :: state(N_NODES,N_NODES,N_VARS)

    logical :: ok

    real(dp) :: ws
    real(dp) :: tiles(2,2,N_VARS)
    real(dp) :: hiles(2,2,N_VARS)

    real(dp) :: means(N_VARS)
    real(dp) :: highs(N_NODES,N_NODES,N_VARS)

    integer :: level
    integer :: i,j,k,h,ntrial
    integer :: q, p, ii ,jj, mm, nn

    real(dp), parameter :: decrement = 0.98_dp
    real(dp), parameter :: repair_decr = 0.05_dp
    integer, parameter :: nlevels = ${int(math.log(PP_N_NODES,2))}$

    real(dp) :: squeeze

    ok = isvalid(state)
    if (ok) return

# if 0
    do h = 1,N_VARS
        means(h) = 1.0_dp/REAL(N_NODES**N_DIMS,dp)*sum(state(:,:,h))
    end do

    do h = 1,N_VARS
        highs(:,:,h) = state(:,:,h) - means(h)
    end do

    squeeze = 1.0_dp
    do while (squeeze > 0.0_dp)
        squeeze = max(0.0_dp,squeeze - 0.1)

        do h = 1,N_VARS
            state(:,:,h) = means(h) + squeeze*highs(:,:,h)
        end do

        ok = isvalid(state)
        if (ok) exit
    end do

# else

    do level = 1,nlevels
        mm = 2**(level-1)
        nn = 2**(nlevels-level)
        ws = 1.0_dp/REAL(nn**2,dp)

        do p = 1,mm; do q = 1,mm

            tiles = 0.0_dp
            do h = 1,N_VARS
                do jj = 1,2; do ii = 1,2
                    do j = 1 + (2*(p-1) + jj-1)*nn, (2*(p-1) + jj)*nn
                    do i = 1 + (2*(q-1) + ii-1)*nn, (2*(q-1) + ii)*nn
                        tiles(ii,jj,h) = tiles(ii,jj,h) + ws * state(i,j,h)
                    end do
                    end do
                end do; end do
            end do

            ok = isvalid(tiles)
            if (ok) cycle

            do h = 1,N_VARS
                means(h) = 0.25_dp*sum(tiles(:,:,h))
            end do

            do h = 1,N_VARS
                hiles(:,:,h) = tiles(:,:,h) - means(h)
            end do

            squeeze = 1.0_dp
            do while (squeeze > 0.0_dp)
                squeeze = max(0.0_dp,squeeze - repair_decr)

                do h = 1,N_VARS
                    tiles(:,:,h) = means(h) + squeeze*hiles(:,:,h)
                end do

                ok = isvalid(tiles)
                if (ok) exit
            end do

            do h = 1,N_VARS
                do jj = 1,2; do ii = 1,2
                    do j = 1 + (2*(p-1) + jj-1)*nn, (2*(p-1) + jj)*nn
                    do i = 1 + (2*(q-1) + ii-1)*nn, (2*(q-1) + ii)*nn
                        state(i,j,h) = means(h) + squeeze * (state(i,j,h) - means(h))
                    end do
                    end do
                end do; end do
            end do

        end do; end do !! tiles

    end do !! levels

    ok = isvalid(state)

# endif

end function

pure function kernel_validate(state) result(ok)

    use equations_mod, only: isvalid

    real(dp), intent(in) :: state(N_NODES,N_NODES,N_VARS)

    logical :: ok

    ok = isvalid(state)
        
end function

end module
