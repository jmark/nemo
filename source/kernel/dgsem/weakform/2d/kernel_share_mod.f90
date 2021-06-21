module kernel_share_mod

use constants_mod

use mesh_const_mod, only: MESH_HALF
use mesh_const_mod, only: MESH_CHILDREN

real(dp), parameter :: eps = 10*TOL

real(dp), save :: dtMin
real(dp), save :: lmax

!! -------------------------------------------------------------------------- !!

real(dp), save :: nodes(N_NODES)
real(dp), save :: weights(N_NODES)

!! -------------------------------------------------------------------------- !!
!! prolongate to boundaries

real(dp), save :: lagrItplVecM(N_NODES)
real(dp), save :: lagrItplVecP(N_NODES)

!! ------------------------------------------------------------------------- !!
!! weak-form DG operators

real(dp), save :: dgDiffMat(N_NODES,N_NODES)
real(dp), save :: dgDiffTam(N_NODES,N_NODES)

real(dp), save :: dgSurfMat(N_NODES,N_NODES)
real(dp), save :: dgSurfTam(N_NODES,N_NODES)

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
!! AMR operators

real(dp), save :: toMidpointsMat(N_NODES,N_NODES)
real(dp), save :: toMidpointsTam(N_NODES,N_NODES)

contains

subroutine init

    use kernel_operators_mod
    use share_mod, only: mesh
    use linalg_mod, only: invert
    use utils_mod, only: krondelta
    use interpol_mod, only: NodesAndWeights
    use lagrange_mod, only: lagrange_polynomial
    use interpol_mod, only: EquidistantNodesAndWeights

    use utils_mod, only: kronproduct
    use utils_mod, only: eyeproduct

    use interpol_mod, only: LegendreGaussNodesAndWeights
    use interpol_mod, only: LegendreGaussLobattoNodesAndWeights

    use kernel_operators_mod
    use linalg_mod, only: invert
    use interpol_mod, only: InterfaceNodes
    use interpol_mod, only: mkVanderMondeMatrix
    use interpol_mod, only: EquidistantNodesAndWeights

    use utils_mod, only: eyeproduct
    use interpol_mod, only: DiffMatrix
    use interpol_mod, only: NodesAndWeights
    use interpol_mod, only: LegendreGaussNodesAndWeights
    use interpol_mod, only: LegendreGaussLobattoNodesAndWeights

    use lagrange_mod, only: lagrange_polynomial

    use interpol_mod, only: LagrangePolynomialDerivative
    use interpol_mod, only: LagrangePolynomialDerivativeMatrix

    use linalg_mod, only: shiftMatrix

    use utils_mod, only: pp

    use status_mod, only: raise_error

    logical :: ok

    integer :: i,j,k,n,q,f,p,o, d ,qq

    real(dp) :: dgnodes(N_NODES), dgweights(N_NODES)
    real(dp) :: fvnodes(N_NODES), fvweights(N_NODES)

    !! --------------------------------------------------------------------- !!

    call LegendreGaussNodesAndWeights(N_NODES-1,dgnodes,dgweights)
    call EquidistantNodesAndWeights(N_NODES-1,fvnodes,fvweights)

    nodes = 0.5_dp*dgnodes
    weights = 0.5_dp*dgweights

    !! --------------------------------------------------------------------- !!
    !! DG operators
   
    call kernel_operators_WeakDGDiffMatrix(N_NODES,dgnodes,dgweights,dgDiffMat)
    dgDiffTam = transpose(dgDiffMat)

    ! call pp('dgdiffMat', dgdiffMat)

    call kernel_operators_SurfaceMatrix(N_NODES,dgnodes,dgweights,dgsurfMat)
    dgsurfTam = transpose(dgsurfMat)

    ! call pp('dgdiffMat', dgsurfMat)

    !! --------------------------------------------------------------------- !!
    !! Boundary Interpolation Operators

    lagrItplVecM = 0.0
    lagrItplVecP = 0.0

    do j = 1,N_NODES
        lagrItplVecM(j) = lagrange_polynomial(N_NODES,dgnodes,j,-1.0_dp)
        lagrItplVecP(j) = lagrange_polynomial(N_NODES,dgnodes,j, 1.0_dp)
    end do

    ! call pp('lagrItplVecP', lagrItplVecP)
    ! call pp('lagrItplVecM', lagrItplVecM)

    !! --------------------------------------------------------------------- !!
    !! Visulization Matrix

    toMidpointsMat = 0.0

    do i = 1,N_NODES
        do j = 1,N_NODES
            toMidpointsMat(i,j) = lagrange_polynomial(N_NODES,dgnodes,j,fvnodes(i))
        end do
    end do

    toMidpointsTam = transpose(toMidpointsMat)

    ! call pp('toMidpointsMat', toMidpointsMat)
    ! call raise_error('debugging')

end subroutine

end module
