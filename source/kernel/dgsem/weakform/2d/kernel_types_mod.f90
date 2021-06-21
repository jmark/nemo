module kernel_types_mod

use constants_mod
use mesh_const_mod

type :: payload_hydro_t

    real(dp) :: state(N_NODES,N_NODES,N_VARS)

end type

type :: auxiliary_hydro_t

    real(dp) :: facevars(N_NODES,N_VARS,MESH_FACES)
    real(dp) :: faceflux(N_NODES,N_VARS,MESH_FACES)

    !real(dp) :: faceavgs(N_NODES,N_VARS,    MESH_FACES)
 
end type

end module
