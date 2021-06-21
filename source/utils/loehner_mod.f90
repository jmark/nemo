module loehner_mod

use constants_mod

interface calc_loehner
    procedure calc_loehner_1d
    procedure calc_loehner_2d
    procedure calc_loehner_3d
end interface

real(dp), parameter :: eps = 0.01_dp

contains

pure function calc_loehner_1d(u) result(r)

    real(dp), intent(in)    :: u(N_NODES)
    real(dp)                :: r

    integer :: i

    real(dp) :: xtemp(2:N_NODES-1)

    do i = 2,N_NODES-1
        xtemp(i) = abs(u(i-1) - 2.0_dp*u(i) + u(i+1)) / &
            (abs(u(i+1)-u(i)) + abs(u(i)-u(i-1)) + &
            eps*(abs(u(i+1)) + 2.0_dp*abs(u(i)) + abs(u(i-1))))
    end do

    r = maxval(xtemp)

end function

pure function calc_loehner_2d(u) result(r)

    real(dp), intent(in)    :: u(N_NODES,N_NODES)
    real(dp)                :: r

    integer :: i,j

    real(dp) :: xtemp(2:N_NODES-1,1:N_NODES)
    real(dp) :: ytemp(1:N_NODES,2:N_NODES-1)

    do j = 1,N_NODES; do i = 2,N_NODES-1;
        xtemp(i,j) = abs(u(i-1,j) - 2.0_dp*u(i,j) + u(i+1,j)) / &
            (abs(u(i+1,j)-u(i,j)) + abs(u(i,j)-u(i-1,j)) + &
            eps*(abs(u(i+1,j)) + 2.0_dp*abs(u(i,j)) + abs(u(i-1,j))))
    end do; end do;

    do j = 2,N_NODES-1; do i = 1,N_NODES;
        ytemp(i,j) = abs(u(i,j-1) - 2.0_dp*u(i,j) + u(i,j+1)) / &
            (abs(u(i,j+1)-u(i,j)) + abs(u(i,j)-u(i,j-1)) + &
            eps*(abs(u(i,j+1)) + 2.0_dp*abs(u(i,j)) + abs(u(i,j-1))))
    end do; end do;

    r = max(maxval(xtemp),maxval(ytemp))

end function

pure function calc_loehner_3d(u) result(r)

    real(dp), intent(in)    :: u(N_NODES,N_NODES,N_NODES)
    real(dp)                :: r

    integer :: i,j,k

    real(dp) :: xtemp(2:N_NODES-1,1:N_NODES,1:N_NODES)
    real(dp) :: ytemp(1:N_NODES,2:N_NODES-1,1:N_NODES)
    real(dp) :: ztemp(1:N_NODES,1:N_NODES,2:N_NODES-1)

    do k = 1,N_NODES; do j = 1,N_NODES; do i = 2,N_NODES-1;
        xtemp(i,j,k) = abs(u(i-1,j,k) - 2.0_dp*u(i,j,k) + u(i+1,j,k)) / &
            (abs(u(i+1,j,k)-u(i,j,k)) + abs(u(i,j,k)-u(i-1,j,k)) + &
            eps*(abs(u(i+1,j,k)) + 2.0_dp*abs(u(i,j,k)) + abs(u(i-1,j,k))))
    end do; end do; end do;

    do k = 1,N_NODES; do j = 2,N_NODES-1; do i = 1,N_NODES;
        ytemp(i,j,k) = abs(u(i,j-1,k) - 2.0_dp*u(i,j,k) + u(i,j+1,k)) / &
            (abs(u(i,j+1,k)-u(i,j,k)) + abs(u(i,j,k)-u(i,j-1,k)) + &
            eps*(abs(u(i,j+1,k)) + 2.0_dp*abs(u(i,j,k)) + abs(u(i,j-1,k))))
    end do; end do; end do;

    do k = 2,N_NODES-1; do j = 1,N_NODES; do i = 1,N_NODES;
        ztemp(i,j,k) = abs(u(i,j,k-1) - 2.0_dp*u(i,j,k) + u(i,j,k+1)) / &
            (abs(u(i,j,k+1)-u(i,j,k)) + abs(u(i,j,k)-u(i,j,k-1)) + &
            eps*(abs(u(i,j,k+1)) + 2.0_dp*abs(u(i,j,k)) + abs(u(i,j,k-1))))
    end do; end do; end do;

    r = max(maxval(xtemp),maxval(ytemp),maxval(ztemp))

end function

!! Some experimental stuff.
!! pure function calc_loehner(u)
!! 
!!     real(dp), intent(in)    :: u(N_NODES,N_NODES)
!!     real(dp)                :: calc_loehner
!! 
!!     real(dp), parameter :: eps = 0.2_dp
!! 
!!     real(dp) :: stencil(-1:1,-1:1)
!!     real(dp) :: xtemp(-1:1)
!!     real(dp) :: ytemp(-1:1)
!! 
!!     integer :: i,j
!! 
!!     stencil = matmul(fv8fv3Mat,matmul(u,fv8fv3Tam))
!! 
!!     forall (j = -1:1)
!!         xtemp(j) = abs(stencil(-1,j) - 2.0_dp*stencil(0,j) + stencil(1,j)) / &
!!                  (abs(stencil(1,j)-stencil(0,j)) + abs(stencil(0,j)-stencil(-1,j)) + &
!!                  eps*(abs(stencil(1,j)) + 2.0_dp*abs(stencil(0,j)) + abs(stencil(-1,j))))
!!     end forall
!! 
!!     forall (i = -1:1)
!!          ytemp(i) = abs(stencil(i,-1) - 2.0_dp*stencil(i,0) + stencil(i,1)) / &
!!              (abs(stencil(i,1)-stencil(i,0)) + abs(stencil(i,0)-stencil(i,-1)) + &
!!              eps*(abs(stencil(i,1)) + 2.0_dp*abs(stencil(i,0)) + abs(stencil(i,-1))))
!!      end forall
!! 
!!     ! calc_loehner = max(maxval(xtemp),maxval(ytemp))
!!     calc_loehner = max(sum(xtemp)/3,sum(ytemp)/3)
!! 
!! end function


end module
