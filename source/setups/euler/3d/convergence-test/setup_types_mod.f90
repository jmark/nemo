module setup_types_mod

use constants_mod

type :: auxiliary_manufacsol_t

    real(dp) :: u0(N_NODES,N_NODES,N_NODES,N_VARS)
    real(dp) :: st(N_NODES,N_NODES,N_NODES,N_VARS)

end type

end module
