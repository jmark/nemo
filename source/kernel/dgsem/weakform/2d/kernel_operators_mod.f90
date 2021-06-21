module kernel_operators_mod

use constants_mod

contains

subroutine kernel_operators_lagrangeProlongationMatrixDGFV(N, dgnodes, fvnodes, matM, matP)

    use lagrange_mod, only: lagrange_polynomial

    integer, intent(in)     :: N
    real(dp), intent(in)    :: dgnodes(N)
    real(dp), intent(in)    :: fvnodes(N)
    real(dp), intent(out)   :: matM(N,N,N)
    real(dp), intent(out)   :: matP(N,N,N)

    integer :: i,j,k

    do k = 1,N
    do j = 1,N
    do i = 1,N
        matM(i,j,k) = lagrange_polynomial(N, dgnodes, i, -1.0_dp)*lagrange_polynomial(N, dgnodes, j, fvnodes(k))
        matP(i,j,k) = lagrange_polynomial(N, dgnodes, i,  1.0_dp)*lagrange_polynomial(N, dgnodes, j, fvnodes(k))
    end do
    end do
    end do

end subroutine

subroutine kernel_operators_DiffMatrix(N, dgnodes, mat)

    use utils_mod, only: krondelta
    use interpol_mod, only: DiffMatrix

    integer, intent(in)     :: N
    real(dp), intent(in)    :: dgnodes(N_NODES)
    real(dp), intent(out)   :: mat(N_NODES,N_NODES)

    real(dp) :: DiffMat(N_NODES,N_NODES)

    real(dp) :: p,q
    integer :: i,j

    DiffMat = 0.0
    call DiffMatrix(N-1, dgnodes(1:N), DiffMat(1:N,1:N))

    p = real(N,dp)
    q = real(N,dp)

    mat = 0.0
    do j = 1,N_NODES; do i = 1,N_NODES
        mat(i,j) = krondelta(ceiling(i/p),ceiling(j/q)) * DiffMat(i-floor((i-1)/p)*N, j-floor((j-1)/q)*N)
    end do; end do

    !! correct element length
    mat = real(N_NODES,dp)/real(N,dp)*mat

end subroutine

# if 0
subroutine kernel_operators_StrongDGDiffMatrix(N, dgnodes, dgweights, mat)

    use utils_mod, only: krondelta
    use interpol_mod, only: DiffMatrix

    integer, intent(in) :: N
    real(dp), intent(in)    :: dgnodes(N_NODES)
    real(dp), intent(in)    :: dgweights(N_NODES)
    real(dp), intent(out)   :: mat(N_NODES,N_NODES)

    real(dp) :: MassMat(N_NODES,N_NODES)
    real(dp) :: DiffMat(N_NODES,N_NODES)

    real(dp) :: p,q
    integer :: i,j

    DiffMat = 0.0
    call DiffMatrix(N-1, dgnodes(1:N), DiffMat(1:N,1:N))

    p = real(N,dp)
    q = real(N,dp)

    mat = 0.0
    do j = 1,N_NODES; do i = 1,N_NODES
        mat(i,j) = krondelta(ceiling(i/p),ceiling(j/q)) * DiffMat(i-floor((i-1)/p)*N, j-floor((j-1)/q)*N)
    end do; end do

    !! correct element length
    mat = real(N_NODES,dp)/real(N,dp)*mat

end subroutine
# endif

subroutine kernel_operators_WeakDGDiffMatrix(N, dgnodes, dgweights, mat)

    use utils_mod, only: krondelta
    use interpol_mod, only: DiffMatrix

    integer, intent(in) :: N
    real(dp), intent(in)    :: dgnodes(N)
    real(dp), intent(in)    :: dgweights(N)
    real(dp), intent(out)   :: mat(N_NODES,N_NODES)

    real(dp) :: temp(N,N)
    real(dp) :: DiffMat(N,N)

    real(dp) :: p,q
    integer :: i,j

    DiffMat = 0.0
    call DiffMatrix(N-1, dgnodes, DiffMat)

    temp = 0.0
    forall (j = 1:N, i = 1:N) temp(i,j) = -2 * dgweights(j)/dgweights(i) * DiffMat(j,i)

    p = real(N,dp)
    q = real(N,dp)

    mat = 0.0
    do j = 1,N_NODES; do i = 1,N_NODES
        mat(i,j) = krondelta(ceiling(i/p),ceiling(j/q)) * temp(i-floor((i-1)/p)*N, j-floor((j-1)/q)*N)
    end do; end do

    !! correct element length
    !mat = real(N_NODES,dp)/real(N,dp)*mat
    mat = -mat/real(N,dp)

# if 0
    write (*,*) 'nodes', qq, dgnodes
    write (*,*) 'weights', qq, dgweights

    write (*,*) qq, 'Riffdmat'
    do i = 1,N_NODES
        write (*,1001) diffrmat(i,:)
    end do
    write (*,*)


    write (*,*) qq, 'xDiffrMat'
    do i = 1,N_NODES
        write (*,1001) xDiffrMat(i,:,qq)
    end do
    write (*,*)

1001 format(PP_N_NODES(f12.3))
# endif

end subroutine

subroutine kernel_operators_WeakDGDiffMatrixV(N, dgnodes, dgweights, mat)

    use utils_mod, only: krondelta
    use interpol_mod, only: DiffMatrix

    integer, intent(in) :: N
    real(dp), intent(in)    :: dgnodes(N)
    real(dp), intent(in)    :: dgweights(N)
    real(dp), intent(out)   :: mat(N,N)

    real(dp) :: temp(N,N)
    real(dp) :: DiffMat(N,N)

    real(dp) :: p,q
    integer :: i,j

    DiffMat = 0.0
    call DiffMatrix(N-1, dgnodes, DiffMat)

    temp = 0.0
    forall (j = 1:N, i = 1:N) temp(i,j) = -2 * dgweights(j)/dgweights(i) * DiffMat(j,i)

    mat = temp

end subroutine

subroutine kernel_operators_AverageMatrix(N, mat)

    use utils_mod, only: krondelta

    integer, intent(in) :: N
    real(dp), intent(out)   :: mat(N_NODES,N_NODES)

    real(dp) :: p,q
    integer :: i,j

    p = real(N,dp)
    q = real(N,dp)

    mat = 0.0
    do j = 1,N_NODES; do i = 1,N_NODES
        mat(i,j) = 1/real(N,dp) * krondelta(ceiling(i/p),ceiling(j/q))
    end do; end do

# if 0
    write (*,*) 'AvgMat', qq
    do i = 1,N_NODES
        write (*,1001) AvgMat(i,:,qq)
    end do
    write (*,*)

1001 format(12(f7.2))
1002 format(24(f12.3))
# endif

end subroutine

!! --------------------------------------------------------------------- !!
!! Modal space operators

subroutine kernel_operators_ModalTrafoMatrix(N, dgnodes, toModal, toNodal)

    use utils_mod, only: krondelta
    use interpol_mod, only: mkVanderMondeMatrix

    integer, intent(in) :: N
    real(dp), intent(in)    :: dgnodes(N)

    real(dp), intent(out)   :: toModal(N_NODES,N_NODES)
    real(dp), intent(out)   :: toNodal(N_NODES,N_NODES)

    real(dp) :: toNodal_(N,N)
    real(dp) :: toModal_(N,N)

    real(dp) :: p,q
    integer :: i,j
    
    p = real(N,dp)
    q = real(N,dp)

    toNodal_ = 0.0
    toModal_ = 0.0
    call mkVanderMondeMatrix(n, dgnodes, toNodal_, toModal_)

    toNodal = 0
    toModal = 0
    do j = 1,N_NODES; do i = 1,N_NODES
        toNodal(i,j) = krondelta(ceiling(i/p),ceiling(j/q)) * toNodal_(i-floor((i-1)/p)*N, j-floor((j-1)/q)*N)
        toModal(i,j) = krondelta(ceiling(i/p),ceiling(j/q)) * toModal_(i-floor((i-1)/p)*N, j-floor((j-1)/q)*N)
    end do; end do
end subroutine

subroutine kernel_operators_ModalTrafoMatrix2(N, dgnodes, toModal, toNodal)

    use interpol_mod, only: mkVanderMondeMatrix

    integer, intent(in)     :: N
    real(dp), intent(in)    :: dgnodes(N)
    real(dp), intent(out)   :: toModal(N,N)
    real(dp), intent(out)   :: toNodal(N,N)

    toNodal = 0.0
    toModal = 0.0
    call mkVanderMondeMatrix(n, dgnodes, toNodal, toModal)

end subroutine

!! --------------------------------------------------------------------- !!
!! trafo operators: DG <-> FV

subroutine kernel_operators_dg2fvMatrix(N, M, dgnodes, dgweights, mat)

    use utils_mod, only: krondelta
    use subcellfv_mod, only: subcellfv_vanderMondeMatrix

    integer, intent(in)     :: N,M
    real(dp), intent(in)    :: dgnodes(N)
    real(dp), intent(in)    :: dgweights(N)
    real(dp), intent(out)   :: mat(M,M)

    real(dp) :: temp(N,N)

    real(dp) :: p,q
    integer :: i,j
    
    temp = 0.0
    call subcellfv_vanderMondeMatrix(N, N, dgnodes, dgweights, temp)
    
    p = real(N,dp)
    q = real(N,dp)

    mat = 0.0
    do j = 1,M; do i = 1,M
        mat(i,j) = krondelta(ceiling(i/p),ceiling(j/q)) * temp(i-floor((i-1)/p)*N, j-floor((j-1)/q)*N)
    end do; end do

end subroutine

subroutine kernel_operators_dg2fvMatrixV(N, M, dgnodes, dgweights, mat)

    use utils_mod, only: krondelta
    use subcellfv_mod, only: subcellfv_vanderMondeMatrix

    integer, intent(in)     :: N,M
    real(dp), intent(in)    :: dgnodes(N)
    real(dp), intent(in)    :: dgweights(N)
    real(dp), intent(out)   :: mat(M,N)

    real(dp) :: p,q
    integer :: i,j
    
    mat = 0.0
    call subcellfv_vanderMondeMatrix(N, M, dgnodes, dgweights, mat)
    
end subroutine

subroutine kernel_operators_refineMatrix(N, dgnodes, dgweights, mat)

    use utils_mod, only: krondelta
    use subcellfv_mod, only: subcellfv_vanderMondeMatrix

    integer, intent(in) :: N
    real(dp), intent(in)    :: dgnodes(N)
    real(dp), intent(in)    :: dgweights(N)
    real(dp), intent(out)   :: mat(2*N_NODES,N_NODES)

    real(dp) :: temp(2*N,N)

    real(dp) :: p,q
    integer :: i,j
    
    temp = 0.0
    call subcellfv_vanderMondeMatrix(N, 2*N, dgnodes, dgweights, temp)
    
    p = real(2*N,dp)
    q = real(N,dp)

    mat = 0.0
    do j = 1,N_NODES; do i = 1,2*N_NODES
        mat(i,j) = krondelta(ceiling(i/p),ceiling(j/q)) * temp(i-floor((i-1)/p)*2*N, j-floor((j-1)/q)*N)
    end do; end do

end subroutine

subroutine kernel_operators_coarsenMatrix(N, dgnodes, dgweights, mat)

    use utils_mod, only: krondelta
    use subcellfv_mod, only: subcellfv_vanderMondeMatrix

    integer, intent(in) :: N
    real(dp), intent(in)    :: dgnodes(N)
    real(dp), intent(in)    :: dgweights(N)
    real(dp), intent(out)   :: mat(N_NODES/2,N_NODES)

    real(dp) :: temp(N/2,N)

    real(dp) :: p,q
    integer :: i,j
    
    temp = 0.0
    call subcellfv_vanderMondeMatrix(N, N/2, dgnodes, dgweights, temp)
    
    p = real(N/2,dp)
    q = real(N,dp)

    mat = 0.0
    do j = 1,N_NODES; do i = 1,N_NODES/2
        mat(i,j) = krondelta(ceiling(i/p),ceiling(j/q)) * temp(i-floor((i-1)/p)*N/2, j-floor((j-1)/q)*N)
    end do; end do

end subroutine

subroutine kernel_operators_overMatrix(N, M, nodesA, nodesB, mat)

    use lagrange_mod, only: lagrange_VanderMondeMatrix

    integer, intent(in)     :: N, M
    real(dp), intent(in)    :: nodesA(N)
    real(dp), intent(in)    :: nodesB(M)
    real(dp), intent(out)   :: mat(M,N)

    call lagrange_VanderMondeMatrix(N, M, nodesA, nodesB, mat)
   
end subroutine

subroutine kernel_operators_revoMatrix(N, M, dgnodes, dgweights, mat)

    use utils_mod, only: krondelta
    use subcellfv_mod, only: subcellfv_vanderMondeMatrix

    integer, intent(in)     :: N, M
    real(dp), intent(in)    :: dgnodes(N)
    real(dp), intent(in)    :: dgweights(N)
    real(dp), intent(out)   :: mat(M,N)

    real(dp) :: temp(M,N)

    real(dp) :: p,q
    integer :: i,j
    
    call subcellfv_vanderMondeMatrix(N, M, dgnodes, dgweights, mat)
    
    !p = real(dp)(M)
    !q = real(dp)(N)

    !mat = 0.0
    !do j = 1,N; do i = 1,M
    !    mat(i,j) = krondelta(ceiling(i/p),ceiling(j/q)) * temp(i-floor((i-1)/p)*M, j-floor((j-1)/q)*N)
    !end do; end do

end subroutine

subroutine kernel_operators_interpolateBoundaryFVMatrix(N, dgnodes, fvnodes, mat)

    use utils_mod, only: krondelta
    use lagrange_mod, only: lagrange_polynomial

    integer, intent(in) :: N
    real(dp), intent(in)    :: dgnodes(N)
    real(dp), intent(in)    :: fvnodes(N)
    real(dp), intent(out)   :: mat(N_NODES,N_NODES)

    real(dp) :: temp(N,N)

    real(dp) :: p,q
    integer :: i,j

    temp = 0.0
    do j = 1,N; do i = 1,N
        temp(i,j) = lagrange_polynomial(n, dgnodes, i, fvnodes(j))
    end do; end do

    p = real(N,dp)
    q = real(N,dp)

    mat = 0.0
    do j = 1,N_NODES; do i = 1,N_NODES
        mat(i,j) = krondelta(ceiling(i/p),ceiling(j/q)) * temp(i-floor((i-1)/p)*n, j-floor((j-1)/q)*n)
    end do; end do

# if 0
        write (*,*) qq, 'dg2fvMat'
        do i = 1,N_NODES
            write (*,1001) dg2fvMat(i,:,qq)
        end do
        write (*,*)

        write (*,*) qq, 'fv2dgMat'
        do i = 1,N_NODES
            write (*,1001) fv2dgMat(i,:,qq)
        end do
        write (*,*)
1001 format(PP_N_NODES(f12.3))
1002 format(24(f12.3))
# endif

end subroutine

subroutine kernel_operators_interpolateBoundaryFVMatrix3(N, dgnodes, mat)

    use utils_mod, only: krondelta
    use utils_mod, only: eyeproduct
    use lagrange_mod, only: lagrange_polynomial
    use interpol_mod, only: EquidistAndInterfaceNodes_3times

    integer, intent(in)     :: N
    real(dp), intent(in)    :: dgnodes(N)
    real(dp), intent(out)   :: mat(3*N_NODES,N_NODES)

    real(dp) :: temp(3*N,N)
    real(dp) :: itfnodes(3*N)
    real(dp) :: itfweights(3*N)

    integer :: i,j

    call EquidistAndInterfaceNodes_3times(3*N, itfnodes, itfweights)

    do i = 1,3*N; do j = 1,N
        temp(i,j) = lagrange_polynomial(n, dgnodes, j, itfnodes(i))
    end do; end do

    mat = eyeproduct(N_NODES/N,temp)

# if 0
    write (*,*) 'lagriptlMat3',N
    do i = 1,3*N_NODES
        write (*,'(8(ES12.4))') mat(i,:)
    end do
    write (*,*)
# endif

end subroutine

subroutine kernel_operators_interpolateBoundaryDGMatrix(N, nodes, MatM, MatP)

    use utils_mod, only: krondelta
    use lagrange_mod, only: lagrange_polynomial
    use interpol_mod, only: InterfaceNodes

    integer, intent(in)     :: N
    real(dp), intent(in)    :: nodes(N)
    real(dp), intent(out)   :: matM(N_NODES-1,N_NODES)
    real(dp), intent(out)   :: matP(N_NODES-1,N_NODES)

    real(dp) :: tempM(N,N)
    real(dp) :: tempP(N,N)

    real(dp) :: p,q
    integer :: i,j

    real(dp), dimension(N) :: xL,xC,xR

    call InterfaceNodes(N,xL,xC,xR)

    tempM = 0.0
    tempP = 0.0
    !do i = 1,N; do j = 1,N
    do j = 1,N
        tempM(1,j) = lagrange_polynomial(N, nodes, j, -1.0_dp)
        tempP(N,j) = lagrange_polynomial(N, nodes, j,  1.0_dp)

        !tempM(i,j) = lagrange_polynomial(N, nodes, j, xL(i))
        !tempP(i,j) = lagrange_polynomial(N, nodes, j, xR(i))
    end do

    p = real(N,dp)
    q = real(N,dp)

    matM = 0.0
    matP = 0.0

    do i = 2,N_NODES; do j = 1,N_NODES
        matM(i-1,j) = krondelta(ceiling(i/p),ceiling(j/q)) * tempM(i-floor((i-1)/p)*n, j-floor((j-1)/q)*n)
    end do; end do

    do i = 1,N_NODES-1; do j = 1,N_NODES
        matP(i,j) = krondelta(ceiling(i/p),ceiling(j/q)) * tempP(i-floor((i-1)/p)*n, j-floor((j-1)/q)*n)
    end do; end do

# if 1
    write (*,*) 'lagritplMatM', N
    do i = 1,N_NODES-1
        write (*,'(8(ES12.4))') matM(i,:)
    end do
    write (*,*)

    write (*,*) 'lagritplMatP', N
    do i = 1,N_NODES-1
        write (*,'(8(ES12.4))') matP(i,:)
    end do
    write (*,*)
# endif

end subroutine

subroutine kernel_operators_interpolateBoundaryDGMatrixV(N, nodes, MatM, MatP)

    use utils_mod, only: krondelta
    use lagrange_mod, only: lagrange_polynomial

    integer, intent(in)     :: N
    real(dp), intent(in)    :: nodes(N)
    real(dp), intent(out)   :: matM(N_NODES-1,N_NODES)
    real(dp), intent(out)   :: matP(N_NODES-1,N_NODES)

    real(dp) :: tempM(N,N)
    real(dp) :: tempP(N,N)

    real(dp) :: p,q
    integer :: i,j

    tempM = 0.0
    tempP = 0.0
    do j = 1,N
        tempM(1,j) = lagrange_polynomial(N, nodes, j, -1.0_dp)
        tempP(N,j) = lagrange_polynomial(N, nodes, j,  1.0_dp)
    end do

    p = real(N,dp)
    q = real(N,dp)

    matM = 0.0
    matP = 0.0

    do i = 2,N_NODES-1; do j = 1,N_NODES
        matM(i-1,j) = krondelta(ceiling(i/p),ceiling(j/q)) * tempM(i-floor((i-1)/p)*n, j-floor((j-1)/q)*n)
    end do; end do

    do i = 1,N_NODES-1; do j = 1,N_NODES
        matP(i,j) = krondelta(ceiling(i/p),ceiling(j/q)) * tempP(i-floor((i-1)/p)*n, j-floor((j-1)/q)*n)
    end do; end do

# if 1
    write (*,*) 'lagritplMatM', N
    do i = 1,N_NODES-1
        write (*,'(8(ES12.4))') matM(i,:)
    end do
    write (*,*)

    write (*,*) 'lagritplMatP', N
    do i = 1,N_NODES-1
        write (*,'(8(ES12.4))') matP(i,:)
    end do
    write (*,*)
# endif

end subroutine

subroutine kernel_operators_interpolateBackMatrix(N, nodes, weights, MatM, MatP)

    use lagrange_mod, only: lagrange_polynomial

    integer, intent(in) :: N
    real(dp), intent(in)    :: nodes(N)
    real(dp), intent(in)    :: weights(N)
    real(dp), intent(out)   :: matM(N_NODES,N_NODES)
    real(dp), intent(out)   :: matP(N_NODES,N_NODES)

    real(dp) :: lagrVecM(N)
    real(dp) :: lagrVecP(N)

    real(dp) :: tempM(N,N)
    real(dp) :: tempP(N,N)

    real(dp) :: p,q
    integer :: i,j

    do i = 1,N
        LagrVecM(i) = lagrange_polynomial(n, nodes, i, -1.0_dp)
        lagrVecP(i) = lagrange_polynomial(n, nodes, i,  1.0_dp)
    end do

    tempM = spread(2./weights*lagrVecM,2,N)
    tempP = spread(2./weights*lagrVecP,2,N)

    p = real(N,dp)
    q = real(N,dp)
    
    matM = 0.0_dp
    matP = 0.0_dp
    do j = 1,N_NODES; do i = 1,N_NODES
        matM(i,j) = tempM(i-floor((i-1)/p)*N, j-floor((j-1)/q)*N)
        matP(i,j) = tempP(i-floor((i-1)/p)*N, j-floor((j-1)/q)*N)
    end do; end do

    matM = real(N_NODES,dp)/real(N,dp)*matM
    matP = real(N_NODES,dp)/real(N,dp)*matP

# if 0
    write (*,*) 'xmLagrM'
    do i = 1,N_NODES
        write (*,1001) xmLagrM(i,:,qq)
    end do
    write (*,*)

    write (*,*) 'xmLagrP'
    do i = 1,N_NODES
        write (*,1001) xmLagrP(i,:,qq)
    end do
    write (*,*)

1001 format(PP_N_NODES(f12.3))
1002 format(24(f12.3))
# endif
end subroutine

subroutine kernel_operators_surfaceMatrix(N, nodes, weights, surfMat)

    use lagrange_mod, only: lagrange_polynomial

    integer, intent(in)     :: N
    real(dp), intent(in)    :: nodes(N)
    real(dp), intent(in)    :: weights(N)

    real(dp), intent(out)   :: surfMat(N_NODES,N_NODES)

    real(dp) :: matM(N_NODES,N_NODES)
    real(dp) :: matP(N_NODES,N_NODES)

    real(dp) :: lagrM(N_NODES)
    real(dp) :: lagrP(N_NODES)

    integer :: i,j

    do i = 1,N_NODES
        lagrM(i) = 2.0_dp/REAL(N,dp)/weights(i)*lagrange_polynomial(N, nodes, i, -1.0_dp)
        lagrP(i) = 2.0_dp/REAL(N,dp)/weights(i)*lagrange_polynomial(N, nodes, i,  1.0_dp)
    end do

    matM = 0.0_dp
    matP = 0.0_dp

    matM(:,1      ) = lagrM
    matP(:,N_NODES) = lagrP

    surfMat = matM - matP
    
end subroutine

subroutine kernel_operators_surfaceMatrix_ext(N,M, nodes, weights, surfMat)

    use utils_mod, only: krondelta
    use share_mod, only: rt => runtime
    use lagrange_mod, only: lagrange_polynomial

    integer, intent(in)     :: N,M
    real(dp), intent(in)    :: nodes(N)
    real(dp), intent(in)    :: weights(N)

    real(dp), intent(out)   :: surfMat(M,M+1)

    real(dp) :: matM(M,M+1)
    real(dp) :: matP(M,M+1)

    real(dp) :: lagrM(M)
    real(dp) :: lagrP(M)

    integer :: i,j

    do i = 1,M
        !lagrM(i) = REAL(N_NODES,dp)/REAL(N,dp)*2./weights(mod(i-1,N)+1)*lagrange_polynomial(N, nodes, mod(i-1,N)+1, -1.0_dp)
        !lagrP(i) = REAL(N_NODES,dp)/REAL(N,dp)*2./weights(mod(i-1,N)+1)*lagrange_polynomial(N, nodes, mod(i-1,N)+1,  1.0_dp)
        !lagrM(i) = 2.0_dp/REAL(N,dp)/weights(mod(i-1,N)+1)*lagrange_polynomial(N, nodes, mod(i-1,N)+1, -1.0_dp)
        !lagrP(i) = 2.0_dp/REAL(N,dp)/weights(mod(i-1,N)+1)*lagrange_polynomial(N, nodes, mod(i-1,N)+1,  1.0_dp)

        lagrM(i) = 2.0_dp/weights(mod(i-1,N)+1)*lagrange_polynomial(N, nodes, mod(i-1,N)+1, -1.0_dp)
        lagrP(i) = 2.0_dp/weights(mod(i-1,N)+1)*lagrange_polynomial(N, nodes, mod(i-1,N)+1,  1.0_dp)
    end do

    matM = 0.0_dp
    matP = 0.0_dp

    do j = 1,M+1; do i = 1,M
        matM(i,j) = lagrM(i) * krondelta(0 + 1 + N * ((i-1)/N), j)
        matP(i,j) = lagrP(i) * krondelta(N + 1 + N * ((i-1)/N), j)
    end do; end do

    surfMat = matM - matP
    
# if 0
    if (rt%mpi%isroot) then
        !write (*,*) 'matM', N
        !do i = 1,N_NODES
        !    write (*,'(9(ES12.4))') matM(i,:)
        !end do
        !write (*,*)
        !write (*,*) 'matP', N
        !do i = 1,N_NODES
        !    write (*,'(9(ES12.4))') matP(i,:)
        !end do
        !write (*,*)

        write (*,*) 'surfMat', N
        do i = 1,N_NODES
            write (*,'(9(ES12.4))') surfMat(i,:)
        end do
        write (*,*)
    end if
# endif

end subroutine

subroutine kernel_operators_surfaceMatrix_split(N, nodes, weights, surfMatM, surfMatP)

    use utils_mod, only: krondelta
    use share_mod, only: rt => runtime
    use lagrange_mod, only: lagrange_polynomial

    integer, intent(in)     :: N
    real(dp), intent(in)    :: nodes(N)
    real(dp), intent(in)    :: weights(N)

    real(dp), intent(out)   :: surfMatM(N_NODES,N_NODES+1)
    real(dp), intent(out)   :: surfMatP(N_NODES,N_NODES+1)

    real(dp) :: matM(N_NODES,N_NODES+1)
    real(dp) :: matP(N_NODES,N_NODES+1)

    real(dp) :: lagrM(N_NODES)
    real(dp) :: lagrP(N_NODES)

    integer :: i,j

    do i = 1,N_NODES
        !lagrM(i) = REAL(N_NODES,dp)/REAL(N,dp)*2./weights(mod(i-1,N)+1)*lagrange_polynomial(N, nodes, mod(i-1,N)+1, -1.0_dp)
        !lagrP(i) = REAL(N_NODES,dp)/REAL(N,dp)*2./weights(mod(i-1,N)+1)*lagrange_polynomial(N, nodes, mod(i-1,N)+1,  1.0_dp)
        lagrM(i) = 2.0_dp/REAL(N,dp)/weights(mod(i-1,N)+1)*lagrange_polynomial(N, nodes, mod(i-1,N)+1, -1.0_dp)
        lagrP(i) = 2.0_dp/REAL(N,dp)/weights(mod(i-1,N)+1)*lagrange_polynomial(N, nodes, mod(i-1,N)+1,  1.0_dp)
    end do

    matM = 0.0_dp
    matP = 0.0_dp

    do j = 1,N_NODES+1; do i = 1,N_NODES
        matM(i,j) = lagrM(i) * krondelta(0 + 1 + N * ((i-1)/N), j)
        matP(i,j) = lagrP(i) * krondelta(N + 1 + N * ((i-1)/N), j)
    end do; end do

    surfMatM =  matM
    surfMatP = -matP
    
end subroutine

!! FV difference operator
subroutine kernel_operators_FVDifferenceMatrix(mat)

    use utils_mod, only: krondelta

    real(dp), intent(out) :: mat(1:N_NODES,0:N_NODES)

    integer :: i,j

    mat = 0.0
    do j = 0,N_NODES; do i = 1,N_NODES
        mat(i,j) = krondelta(i-1,j) - krondelta(i,j)
    end do; end do

    mat = N_NODES * mat

end subroutine

!! FV surface operator
subroutine kernel_operators_FVProlongBackMatrix(mat)

    use utils_mod, only: krondelta

    real(dp), intent(out)   :: mat(N_NODES,2*N_NODES,2)

    integer :: i,j

    mat = 0.0

    do j = 1,2*N_NODES; do i = 1,N_NODES
        mat(i,j,1) = krondelta(i,j)
        mat(i,j,2) = 0.5*(krondelta(i*2-1,j) + krondelta(i*2,j))
    end do; end do

# if 0
    write (*,*) 'prolongBack'
    do i = 1,N_NODES
        write (*,1002) prolongBack(i,:,0,2)
    end do
    write (*,*)

1001 format(PP_N_NODES(f8.3))
1002 format(16(f6.2))
# endif

end subroutine

!! FV surface operator
subroutine kernel_operators_FVProlongBackMatrixN(N,mat)

    use utils_mod, only: krondelta

    integer, intent(in)     :: N
    real(dp), intent(out)   :: mat(N,2*N)

    integer :: i,j

    mat = 0.0
    do j = 1,2*N; do i = 1,N
        mat(i,j) = 0.5*(krondelta(i*2-1,j) + krondelta(i*2,j))
    end do; end do

end subroutine

!! FV simpson rule operator
subroutine kernel_operators_SimpsonRuleMatrixN(N,mat)

    integer, intent(in)     :: N
    real(dp), intent(out)   :: mat(N,3*N)

    integer :: i,j

    real(dp), parameter :: coeffs(3) = (/1.0_dp,4.0_dp,1.0_dp/) / 6.0_dp
    !real(dp), parameter :: coeffs(3) = (/0.0_dp,1.0_dp,0.0_dp/)

    mat = 0.0
    do i = 1,N; do j = (i-1)*3 + 2,i*3,3
        mat(i,j-1) = coeffs(1)
        mat(i,j+0) = coeffs(2)
        mat(i,j+1) = coeffs(3)
    end do; end do

end subroutine

pure subroutine kernel_operators_InterfaceNodes(N,xs,xL,xC,xR,x3)

    implicit none

    integer,intent(in)      :: N
    real(dp), intent(in)    :: xs(3)
    real(dp),intent(out)    :: xL(N),xC(N),xR(N),x3(3*N)

    integer :: i

    do i = 1,N
        xL(i) = -1.0_dp + 2.0_dp*(real(i,dp)-0.5_dp)/real(N,dp) + xs(1)/real(N,dp)
        xC(i) = -1.0_dp + 2.0_dp*(real(i,dp)-0.5_dp)/real(N,dp) + xs(2)/real(N,dp)
        xR(i) = -1.0_dp + 2.0_dp*(real(i,dp)-0.5_dp)/real(N,dp) + xs(3)/real(N,dp)
    end do

    do i = 1,N
        x3(3*(i-1)+1) = xL(i)
        x3(3*(i-1)+2) = xC(i)
        x3(3*(i-1)+3) = xR(i)
    end do

end subroutine

subroutine kernel_operators_faceMatrix(N,L, M,P, xnodes, ynodes, mat)

    !! asuming nodes are within [-0.5,0.5]

    use utils_mod, only: krondelta
    use lagrange_mod, only: lagrange_polynomial

    integer, intent(in)     :: N,L, M,P
    real(dp), intent(in)    :: xnodes(L)
    real(dp), intent(in)    :: ynodes(P)
    real(dp), intent(out)   :: mat(M,N)

    real(dp) :: xbd(-1:1,N/L) !! left,center,right
    real(dp) :: ybd(-1:1,M/P) !! left,center,right

    real(dp) :: xxnodes(N)
    real(dp) :: yynodes(M)

    real(dp) :: dx, dxx, dy, prefac, ofsL, ofsR
    integer :: i,j, ii,jj,iii

    integer :: Nx,Ny

    real(dp), parameter :: tol = 1e-9

    Nx = N/L
    Ny = M/P

    dx = 1.0_dp/real(Nx,dp)
    dy = 1.0_dp/real(Ny,dp)

    !! compute boundaries of input sub-elements
    do i = 1,Nx
        xbd(-1,i) = -0.5_dp + (i-1.0_dp)*dx
        xbd( 0,i) = -0.5_dp + (i-0.5_dp)*dx
        xbd( 1,i) = -0.5_dp + (i       )*dx
    end do

    !! compute boundaries of output sub-elements
    do i = 1,Ny
        ybd(-1,i) = -0.5_dp + (i-1.0_dp)*dy
        ybd( 0,i) = -0.5_dp + (i-0.5_dp)*dy
        ybd( 1,i) = -0.5_dp + (i       )*dy
    end do

    !! compute destination nodes
    do i = 1,Ny
        do ii = 1,P
            yynodes((i-1)*P+ii) = ybd(0,i) + dy*ynodes(ii)
        end do
    end do

    !! compute interpolation matrix
    mat = 0.0_dp
    do iii = 1,Ny
        do ii= 1,P
            i = (iii-1)*P+ii

            if (ii == 1) then
                ofsL = tol
                ofsR = 0.0_dp

            else if (ii == P) then
                ofsL = 0.0_dp
                ofsR = -tol
            end if

            do j = 1,Nx

                if (abs(ybd( 1,iii) - yynodes(i)) > tol .and. abs(xbd(1,j) - yynodes(i)) <= tol) then
                    do jj = 1,L
                        mat(i,(j-1)*L + jj) = 0.5*lagrange_polynomial(L, xnodes, jj, 0.5_dp)
                    end do

                    cycle
                endif

                if (abs(ybd(-1,iii) - yynodes(i)) > tol .and. abs(xbd(-1,j) - yynodes(i)) <= tol) then
                    do jj = 1,L
                        mat(i,(j-1)*L + jj) = 0.5*lagrange_polynomial(L, xnodes, jj, -0.5_dp)
                    end do

                    cycle
                endif

                if (xbd(-1,j) < yynodes(i) + ofsL .and. yynodes(i) + ofsR < xbd(1,j)) then
                        
                    dxx = abs(xbd(-1,j) - yynodes(i))/dx
                    do jj = 1,L
                        mat(i,(j-1)*L + jj) = lagrange_polynomial(L, xnodes, jj, -0.5_dp + dxx)
                    end do

                    exit
                end if

            end do
        end do
    end do

# if 0
    write (*,*) N,L,M,P,'faceMat'
    do i = 1,M
        write (*,'(99(ES12.4))') Mat(i,:)
    end do
    write (*,*)
# endif

end subroutine

subroutine kernel_operators_interpolationMatrix(N,M,Q, xs,ws, ys,vs, mat)

    use lagrange_mod, only: lagrange_polynomial

    integer, intent(in)     :: N,M,Q        !! # of xs-nodes, ys-nodes, sub-elements
    real(dp), intent(in)    :: xs(N),ws(N)  !! nodes, weights of source element
    real(dp), intent(in)    :: ys(M),vs(M)  !! nodes, weights of target sub-elements
    real(dp), intent(out)   :: mat(M,N,Q)

    real(dp) :: xL,xR   !! left,right boundaries
    real(dp) :: mu(Q)   !! center nodes of sub-elements
    real(dp) :: xi(N,Q) !! scaled quadrature nodes of 'xs' around center nodes 'mu'

    integer :: i,j,k

    !! works for reference elements of arbitrary size
    xL = -0.5*sum(vs)
    xR = -xL
    
    !! make body-centered linear space
    forall (k = 1:Q)
        mu(k) = xL + (k-0.5)*abs(xR-xL)/Q
    end forall

    !! interpolation nodes within each sub-element of master element
    forall (k = 1:Q, i = 1:M)
        xi(i,k) = mu(k) + ys(i)/Q
    end forall

    !! interpolation matrices
    forall (k = 1:Q, j = 1:N, i = 1:M)
        mat(i,j,k) = lagrange_polynomial(N,xs,j,xi(i,k))
    end forall

# if 0
    do k = 1,Q
        write (*,*) 'interpolationMat',N,M,Q,k
        do i = 1,M
            write (*,'(99(ES12.4))') mat(i,:,k)
        end do
        write (*,*)
    end do
# endif

end subroutine

subroutine kernel_operators_projectionMatrix(N,M,Q, xs,ws, ys,vs, mat)

    use lagrange_mod, only: lagrange_polynomial

    integer, intent(in)     :: N,M,Q        !! # of xs-nodes, ys-nodes, sub-elements
    real(dp), intent(in)    :: xs(N),ws(N)  !! nodes, weights of sub-elements
    real(dp), intent(in)    :: ys(M),vs(M)  !! nodes, weights of target element
    real(dp), intent(out)   :: mat(M,N,Q)

    real(dp) :: xL,xR   !! left,right boundaries
    real(dp) :: mu(Q)   !! center nodes of sub-elements
    real(dp) :: xi(N,Q) !! scaled quadrature nodes of 'xs' around center nodes 'mu'

    integer :: i,j,k

    !! works for reference elements of arbitrary size
    xL = -0.5*sum(vs)
    xR = -xL
    
    !! make body-centered linear space
    forall (k = 1:Q)
        mu(k) = xL + (k-0.5)*abs(xR-xL)/Q
    end forall

    !! interpolation nodes within each sub-element of master element
    forall (k = 1:Q, j = 1:N)
        xi(j,k) = mu(k) + xs(j)/Q
    end forall

    !! projection matrices
    forall (k = 1:Q, j = 1:N, i = 1:M)
        mat(i,j,k) = lagrange_polynomial(M,ys,i,xi(j,k))*ws(j)/vs(i)/Q
    end forall

# if 0
    do k = 1,Q
        write (*,*) 'projectionMat',N,M,Q,k
        do i = 1,M
            write (*,'(99(ES12.4))') mat(i,:,k)
        end do
        write (*,*)
    end do
# endif

end subroutine

end module
