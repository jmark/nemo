module setup_analyze_types_mod

use constants_mod
use mesh_const_mod

type :: auxiliary_analyze_t

    real(dp) :: halo(N_NODES,N_NODES,N_VARS,MESH_FACES)
    ! real(dp) :: divB(N_NODES,N_NODES,N_NODES)

end type

end module
