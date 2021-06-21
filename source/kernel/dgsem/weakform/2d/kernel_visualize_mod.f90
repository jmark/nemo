module kernel_visualize_mod

use constants_mod

contains

function kernel_visualize(state) result(visual)

    use kernel_share_mod, only: toMidpointsMat
    use kernel_share_mod, only: toMidpointsTam

    real(dp), intent(in) :: state(N_NODES,N_NODES)
    real(dp)             :: visual(N_NODES,N_NODES)

    visual = matmul(toMidpointsMat,matmul(state,toMidpointsTam))

end function

end module
