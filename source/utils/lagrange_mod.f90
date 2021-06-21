module lagrange_mod

use constants_mod

!integer, parameter :: RP = selected_real_kind(15)
!real, parameter :: PI = 4.0_rp * atan(1.0_rp)
!real, parameter :: TOL = 4.0e-16

CONTAINS

subroutine lagrange_diffMatrix(x,n,d)

    integer, intent(in)                             :: n

    real(dp), dimension(n), intent(in)       :: x
    real(dp), dimension(n,n), intent(out)  :: d

    integer                                         :: i,j
    real(dp), dimension(n)                   :: w

    call lagrange_barycentricweights(x,n,w)

    do i=1,n
        do j=1,n
            if (i.eq.j) cycle
            d(i,j) = w(j) / ( w(i) * (x(i)-x(j)) )
            d(i,i) = d(i,i) - d(i,j)
        end do
    end do

end subroutine

subroutine lagrange_barycentricWeights(x,n,w)

    integer, intent(in)         :: N
    real(dp), intent(in)   :: x(N)
    real(dp), intent(out)  :: w(N)

    integer :: i,j

    w = 1

    do j=1,n
        do i=1,n
            if (j.eq.i) cycle
            w(j) = w(j) / (x(j)-x(i))
        end do
    end do

end subroutine

pure function lagrange_polynomial(N, nodes, j, x) result(lp)
    
    integer, intent(in) :: N
    real(dp), intent(in)    :: nodes(N)
    integer, intent(in) :: j
    real(dp), intent(in)    :: x

    real(dp)    :: lp
    integer :: i

    lp = 1.0_dp

    if (abs(x-nodes(j)) < 10.0_dp*TOL) return !! Kronecker property

    do i = 1,N
        if (i == j) cycle
        lp = lp * (x - nodes(i))/(nodes(j) - nodes(i))
    end do
     
end function

pure function lagrange_basis(N, nodes, x) result(lp)
    
    integer, intent(in) :: N
    real(dp), intent(in)    :: nodes(N)
    real(dp), intent(in)    :: x

    real(dp)    :: lp(N)
    integer :: i

    !! Kronecker property
    if (any(abs(x-nodes) < 10*TOL)) then
        lp = merge(1.0_dp,0.0_dp,abs(x-nodes) < 10.0_dp*TOL)
        return
    end if

    do i = 1,N
        lp(i) = Lagrange_Polynomial(N, nodes, i, x)
    end do
     
end function

pure subroutine lagrange_VanderMondeMatrix(N, M, xs, ys, Vdm)

    integer, intent(in)     :: N, m
    real(dp), intent(in)        :: xs(N), ys(M)
    real(dp), intent(out)       :: Vdm(M,N)

    integer :: i

    do i = 1,M
        Vdm(i,:) = Lagrange_Basis(N, xs, ys(i))
    end do

end subroutine

end module
