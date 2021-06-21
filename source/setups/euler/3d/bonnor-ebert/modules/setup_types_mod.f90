module setup_types_mod

use constants_mod
use mesh_const_mod, only: MESH_VERTS

type :: auxiliary_gravity_t
    
    !! semi-analytic
    real(dp) :: rpot(N_NODES,N_NODES,N_NODES)
    real(dp) :: racc(N_NODES,N_NODES,N_NODES,N_DIMS)

    !! gravitational potential
    real(dp) :: pot(N_NODES,N_NODES,N_NODES)

    !! gravitational acceleration
    real(dp) :: acc(N_NODES,N_NODES,N_NODES,N_DIMS)

end type

end module
