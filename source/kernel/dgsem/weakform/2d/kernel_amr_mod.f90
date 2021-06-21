module kernel_amr_mod

use constants_mod

interface kernel_amr_coarsen
    procedure :: kernel_amr_coarsen_one_var
    procedure :: kernel_amr_coarsen_all_var
end interface

contains

subroutine kernel_amr_refine(incoming,outgoing)

    use mesh_const_mod, only: MESH_CHILDREN

    use kernel_share_mod, only: refineMat
    use kernel_share_mod, only: refineTam

    real(dp), intent(in)  :: incoming(N_NODES,N_NODES,N_VARS)
    real(dp), intent(out) :: outgoing(N_NODES,N_NODES,N_VARS,MESH_CHILDREN)

    integer :: h,f

    do f = 1,MESH_CHILDREN
        do h = 1,N_VARS
            outgoing(:,:,h,f) = &
                matmul(refineMat(:,:,X_DIR,f),matmul(incoming(:,:,h),refineTam(:,:,Y_DIR,f)))
        end do
    end do

end subroutine

subroutine kernel_amr_coarsen_one_var(incoming,outgoing)

    use mesh_const_mod, only: MESH_CHILDREN

    use kernel_share_mod, only: coarseMat
    use kernel_share_mod, only: coarseTam

    real(dp), intent(in)  :: incoming(N_NODES,N_NODES,MESH_CHILDREN)
    real(dp), intent(out) :: outgoing(N_NODES,N_NODES)

    integer :: f

    outgoing = 0.0_dp
    do f = 1,MESH_CHILDREN
        outgoing = outgoing + &
            matmul(coarseMat(:,:,X_DIR,f),matmul(incoming(:,:,f),coarseTam(:,:,Y_DIR,f)))
    end do

end subroutine

subroutine kernel_amr_coarsen_all_var(incoming,outgoing)

    use mesh_const_mod, only: MESH_CHILDREN

    use kernel_share_mod, only: coarseMat
    use kernel_share_mod, only: coarseTam

    real(dp), intent(in)  :: incoming(N_NODES,N_NODES,N_VARS,MESH_CHILDREN)
    real(dp), intent(out) :: outgoing(N_NODES,N_NODES,N_VARS)

    integer :: h,f

    outgoing = 0.0_dp
    do f = 1,MESH_CHILDREN
        do h = 1,N_VARS
            outgoing(:,:,h) = outgoing(:,:,h) + &
                matmul(coarseMat(:,:,X_DIR,f),matmul(incoming(:,:,h,f),coarseTam(:,:,Y_DIR,f)))
        end do
    end do

end subroutine

end module
