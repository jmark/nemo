module kernel_share_mod

use constants_mod

!! ------------------------------------------------------------------------- !!
!! AMR operators

real(dp), save :: rfaceMat(N_NODES,N_NODES,4)
real(dp), save :: rfaceTam(N_NODES,N_NODES,4)

real(dp), save :: cfaceMat(N_NODES,N_NODES,4)
real(dp), save :: cfaceTam(N_NODES,N_NODES,4)

real(dp), save :: refineMat(N_NODES,N_NODES,3,8)
real(dp), save :: coarseMat(N_NODES,N_NODES,3,8)
                                             
real(dp), save :: refineTam(N_NODES,N_NODES,3,8)
real(dp), save :: coarseTam(N_NODES,N_NODES,3,8)

!! ------------------------------------------------------------------------- !!

real(dp), save :: nodes(N_NODES)
real(dp), save :: weights(N_NODES)

contains

subroutine hook_init_400

    use linalg_mod, only: shiftmatrix

    logical :: ok

    integer :: i,j,k,n,q,f,p,o, d ,qq

    real(dp) :: tempMat(N_NODES,N_NODES)

    !! --------------------------------------------------------------------- !!

    do i = 1,N_NODES
        nodes(i) = -0.5_dp + (i-1)/REAL(N_NODES,dp) + 0.5_dp/REAL(N_NODES,dp)
    end do

    weights = 1.0_dp / REAL(N_NODES,dp)

    !! --------------------------------------------------------------------- !!

    rfaceMat = 0.0_dp

    tempMat = 0.0_dp
    tempMat(1,1) = 1.0_dp
    tempMat(2,1) = 1.0_dp

    do f = 1,4
        do i = 1,N_NODES/2
            rfaceMat(:,:,f) = rfaceMat(:,:,f) + shiftMatrix(2*(i-1),(f-1)*N_NODES/2+(i-1),tempMat)
        end do
    end do

    do f = 1,4
        cfaceMat(:,:,f) = 0.5_dp*transpose(rfaceMat(:,:,f))
    end do

    rfaceMat(:,:,3) = rfaceMat(:,:,1)
    rfaceMat(:,:,4) = rfaceMat(:,:,2)

    rfaceTam(:,:,1) = transpose(rfaceMat(:,:,1))
    rfaceTam(:,:,2) = transpose(rfaceMat(:,:,1))
    rfaceTam(:,:,3) = transpose(rfaceMat(:,:,2))
    rfaceTam(:,:,4) = transpose(rfaceMat(:,:,2))

    cfaceMat(:,:,3) = cfaceMat(:,:,1)
    cfaceMat(:,:,4) = cfaceMat(:,:,2)

    cfaceTam(:,:,1) = transpose(cfaceMat(:,:,1))
    cfaceTam(:,:,2) = transpose(cfaceMat(:,:,1))
    cfaceTam(:,:,3) = transpose(cfaceMat(:,:,2))
    cfaceTam(:,:,4) = transpose(cfaceMat(:,:,2))

    refineMat(:,:,1,1) = rfaceMat(:,:,1)
    refineMat(:,:,2,1) = rfaceMat(:,:,1)
    refineMat(:,:,3,1) = rfaceMat(:,:,1)
                     
    refineMat(:,:,1,2) = rfaceMat(:,:,2)
    refineMat(:,:,2,2) = rfaceMat(:,:,1)
    refineMat(:,:,3,2) = rfaceMat(:,:,1)
                     
    refineMat(:,:,1,3) = rfaceMat(:,:,1)
    refineMat(:,:,2,3) = rfaceMat(:,:,2)
    refineMat(:,:,3,3) = rfaceMat(:,:,1)
                     
    refineMat(:,:,1,4) = rfaceMat(:,:,2)
    refineMat(:,:,2,4) = rfaceMat(:,:,2)
    refineMat(:,:,3,4) = rfaceMat(:,:,1)
                     
    refineMat(:,:,1,5) = rfaceMat(:,:,1)
    refineMat(:,:,2,5) = rfaceMat(:,:,1)
    refineMat(:,:,3,5) = rfaceMat(:,:,2)
                     
    refineMat(:,:,1,6) = rfaceMat(:,:,2)
    refineMat(:,:,2,6) = rfaceMat(:,:,1)
    refineMat(:,:,3,6) = rfaceMat(:,:,2)
                     
    refineMat(:,:,1,7) = rfaceMat(:,:,1)
    refineMat(:,:,2,7) = rfaceMat(:,:,2)
    refineMat(:,:,3,7) = rfaceMat(:,:,2)
                     
    refineMat(:,:,1,8) = rfaceMat(:,:,2)
    refineMat(:,:,2,8) = rfaceMat(:,:,2)
    refineMat(:,:,3,8) = rfaceMat(:,:,2)

    coarseMat(:,:,1,1) = cfaceMat(:,:,1)
    coarseMat(:,:,2,1) = cfaceMat(:,:,1)
    coarseMat(:,:,3,1) = cfaceMat(:,:,1)
                     
    coarseMat(:,:,1,2) = cfaceMat(:,:,2)
    coarseMat(:,:,2,2) = cfaceMat(:,:,1)
    coarseMat(:,:,3,2) = cfaceMat(:,:,1)
                     
    coarseMat(:,:,1,3) = cfaceMat(:,:,1)
    coarseMat(:,:,2,3) = cfaceMat(:,:,2)
    coarseMat(:,:,3,3) = cfaceMat(:,:,1)
                     
    coarseMat(:,:,1,4) = cfaceMat(:,:,2)
    coarseMat(:,:,2,4) = cfaceMat(:,:,2)
    coarseMat(:,:,3,4) = cfaceMat(:,:,1)
                     
    coarseMat(:,:,1,5) = cfaceMat(:,:,1)
    coarseMat(:,:,2,5) = cfaceMat(:,:,1)
    coarseMat(:,:,3,5) = cfaceMat(:,:,2)
                     
    coarseMat(:,:,1,6) = cfaceMat(:,:,2)
    coarseMat(:,:,2,6) = cfaceMat(:,:,1)
    coarseMat(:,:,3,6) = cfaceMat(:,:,2)
                     
    coarseMat(:,:,1,7) = cfaceMat(:,:,1)
    coarseMat(:,:,2,7) = cfaceMat(:,:,2)
    coarseMat(:,:,3,7) = cfaceMat(:,:,2)
                     
    coarseMat(:,:,1,8) = cfaceMat(:,:,2)
    coarseMat(:,:,2,8) = cfaceMat(:,:,2)
    coarseMat(:,:,3,8) = cfaceMat(:,:,2)

    do f = 1,8
    do d = 1,3
        refineTam(:,:,d,f) = transpose(refineMat(:,:,d,f))
        coarseTam(:,:,d,f) = transpose(coarseMat(:,:,d,f))
    end do
    end do

end subroutine

end module
