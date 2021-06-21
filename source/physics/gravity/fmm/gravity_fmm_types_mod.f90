module gravity_fmm_types_mod

use constants_mod

type :: auxiliary_gravity_t
    
    !! gravitational potential
    real(dp) :: pot(N_NODES,N_NODES,N_NODES)

    !! gravitational acceleration
    real(dp) :: acc(N_NODES,N_NODES,N_NODES,N_DIMS)

end type

end module
