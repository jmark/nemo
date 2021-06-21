module trafo_mod

use constants_mod

contains

pure function minmod(uL, uC, uR) result(slope)

    real(dp), intent(in) :: uL, uC, uR
    real(dp) :: slope, dL, dR

    dL = uC-uL
    dR = uR-uC

    slope = sign(1.0_dp,dL)*merge(min(abs(dL),abs(dR)), 0.0_dp, dL*dR > 0)

end function

pure subroutine linear_scale_up_1d(N, src, snk)
    
    integer, intent(in) :: N
    real(dp), intent(in)    :: src(0:N+1)
    real(dp), intent(out)   :: snk(2*N)

    integer :: i
    real(dp) :: slope
    real(dp), parameter :: beta = 0.5_dp

    do i = 2,2*N,2
        associate(uL => src(i/2-1), uC => src(i/2), uR => src(i/2+1))
        slope = minmod(uL, uC, uR)
        snk(i-1) = uC - beta*slope
        snk(i+0) = uC + beta*slope
        end associate
    end do

end subroutine

pure subroutine linear_scale_up_2d(N, src, snk)
    
    integer, intent(in) :: N
    real(dp), intent(in)    :: src(0:N+1,0:N+1)
    real(dp), intent(out)   :: snk(2*N,2*N)

    integer :: i,j
    real(dp) :: xslope, yslope
    real(dp), parameter :: beta = 0.75_dp

    do i = 2,2*N,2
    do j = 2,2*N,2
        associate(uN => src(i/2-1,j/2), uS => src(i/2+1,j/2))
        associate(uW => src(i/2,j/2-1), uE => src(i/2,j/2+1), uC => src(i/2,j/2))
        xslope = minmod(uN, uC, uS)
        yslope = minmod(uW, uC, uE)
        snk(i-1,j-1) = uC - 0.5_dp*beta*(xslope + yslope)
        snk(i-0,j-1) = uC - 0.5_dp*beta*(yslope - xslope)
        snk(i-1,j+0) = uC - 0.5_dp*beta*(xslope - yslope)
        snk(i+0,j+0) = uC + 0.5_dp*beta*(xslope + yslope)
        end associate
        end associate
    end do
    end do

end subroutine

pure subroutine scale_up_2d(bar, foo)
    
    !! R^(N x N) -> R^(2*N x 2*N)

    use constants_mod

    real(dp), intent(in)   :: bar(:,:)
    real(dp), intent(out)  :: foo(:,:)

    integer :: i,j,k,l, ndims(2)

    ndims = shape(bar)

    do j = 1,ndims(2); do i = 1,ndims(1)
        do l = 0,1; do k = 0,1
            foo(2*i-k,2*j-l) = bar(i,j)
        end do; end do
    end do; end do

end subroutine

pure subroutine scale_down_2d(foo, bar)

    !! R^(2*N x 2*N) -> R^(N x N)

    use constants_mod

    real(dp), intent(in)    :: foo(:,:)
    real(dp), intent(out)   :: bar(:,:)

    integer :: i,j,k,l, ndims(2)

    ndims = shape(bar)
    bar   = 0.0

    do j = 1,ndims(2); do i = 1,ndims(1)
        do l = 0,1; do k = 0,1
            bar(i,j) = bar(i,j) + 0.25*foo(2*i-k,2*j-l)
        end do; end do
    end do; end do

end subroutine

pure subroutine scale_up_3d(bar, foo)
    
    !! R^(N x N x N) -> R^(2*N x 2*N x 2*N)

    use constants_mod

    real(dp), intent(in)   :: bar(:,:,:)
    real(dp), intent(out)  :: foo(:,:,:)

    integer :: i,j,k
    integer :: q,r,s
    integer :: ndims(3)

    ndims = shape(bar)

    do k = 1,ndims(3)
    do j = 1,ndims(2)
    do i = 1,ndims(1)
        do s = 0,1
        do r = 0,1
        do q = 0,1
            foo(2*i-q,2*j-r,2*k-s) = bar(i,j,k)
        end do
        end do
        end do
    end do
    end do
    end do

end subroutine

pure subroutine scale_down_3d(foo, bar)

    !! R^(2*N x 2*N x 2*N) -> R^(N x N x N)

    use constants_mod

    real(dp), intent(in)    :: foo(:,:,:)
    real(dp), intent(out)   :: bar(:,:,:)

    integer :: i,j,k
    integer :: q,r,s
    integer :: ndims(3)

    ndims = shape(bar)
    bar   = 0.0

    do k = 1,ndims(3)
    do j = 1,ndims(2)
    do i = 1,ndims(1)
        do s = 0,1
        do r = 0,1
        do q = 0,1
             bar(i,j,k) = bar(i,j,k) + 0.125_dp*foo(2*i-q,2*j-r,2*k-s)
        end do
        end do
        end do
    end do
    end do
    end do

end subroutine

end module
