module subcellfv_mod

use constants_mod

!integer, parameter  :: RP = selected_real_kind(15)
!real(dp), parameter     :: PI = 4.0_rp * atan(1.0_rp)
!real(dp), parameter     :: TOL = 4.0e-16

contains

subroutine subcellfv_midpoints(N, sfvw, sfvx, sfvb)

    integer,intent(in) :: N             !! number of sub-cells
    real(dp),intent(out)   :: sfvw          !! width of sub-volumes in reference space
    real(dp),intent(out)   :: sfvx(N)       !! cell-centers of sub-cells in reference space
    real(dp),intent(out)   :: sfvb(N+1)     !! positions of boundaries of sub-cells in reference space 

    integer :: i

    sfvw  = 2.0_dp/N

    sfvb(1) = -1.0_dp

    do i = 1,N
      sfvb(i+1) =         sfvb(i) + sfvw
      sfvx(i)   = 0.5_dp*(sfvb(i) + sfvb(i+1))
    end do

end subroutine

subroutine subcellfv_vanderMondeMatrix(N, M, nodes, weights, Vdm)

    use lagrange_mod, only: lagrange_VanderMondeMatrix

    integer, intent(in)     :: N, M
    real(dp), intent(in)        :: nodes(N)
    real(dp), intent(in)        :: weights(N)
    real(dp), intent(out)       :: Vdm(M,N)

    real(dp) :: sfvw
    real(dp) :: sfvx(M)
    real(dp) :: sfvb(M+1)

    real(dp) :: sfvVdm(N,N)
    real(dp) :: sfvnodes(N)

    integer :: i,j,k

    Vdm = 0.0_dp

    CALL subcellfv_midpoints(M, sfvw, sfvx, sfvb)

    do i = 1,M

        sfvnodes = sfvb(i) + 0.5_dp*(nodes+1.0_dp) * (sfvb(i+1) - sfvb(i))

        call lagrange_VanderMondeMatrix(N, N, nodes, sfvnodes, sfvVdm)

        do j = 1,N
            do k = 1,N
                Vdm(i,j) = Vdm(i,j) + 0.5_dp * weights(k) * sfvVdm(k,j)
            end do
        end do

    end do

end subroutine

# if 0
subroutine subcellfv_dgmidpoints(N, weights, dgb)

    integer,intent(in) :: N        
    real(dp),intent(in)    :: weights  
    real(dp),intent(out)   :: dgx(N)
    real(dp),intent(out)   :: dgb(N+1)

    integer :: i

    dgb(1) = -1.0_dp

    do i = 1,N
      dgb(i+1) = dgb(i) + weights(i)
      dgx(i)   = 0.5_dp*(dgb(i) + dgb(i+1))
    end do

end subroutine

subroutine subcellfv_reconstruction(N, M, nodes, weights)

    use lagrange_mod, only: lagrange_VanderMondeMatrix

    integer, intent(in)     :: N, M
    real(dp), intent(in)        :: nodes(N)
    real(dp), intent(in)        :: weights(N)
    real(dp), intent(out)       :: Vdm(N,M)

    real(dp) :: sfvw
    real(dp) :: sfvx(M)
    real(dp) :: sfvb(M+1)

    real(dp) :: sfvVdm(N,N)
    real(dp) :: sfvnodes(N)

    real(dp) :: dgb(N)
    real(dp) :: dgb(N+1)

    integer :: i,j,k

    Vdm = 0.0_dp

    call subcellfv_midpoints(M, sfvw, sfvx, sfvb)
    call subcellfv_dgmidpoints(N, dgx, dgb)

    write (*,*) dgx
    write (*,*) dgb

    !do j = 1,M
    !    do i = 1,N
    !        Vdm(i,j) = merge(sfvw, 0.0, 
    !    end do
    !end do

end subroutine
# endif

end module
