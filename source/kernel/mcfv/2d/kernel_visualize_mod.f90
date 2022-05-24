module kernel_visualize_mod

use constants_mod

contains

function kernel_visualize(state) result(visual)

    real(dp), intent(in) :: state(N_NODES,N_NODES)
    real(dp)             :: visual(N_NODES,N_NODES)

    visual = state

end function

end module
