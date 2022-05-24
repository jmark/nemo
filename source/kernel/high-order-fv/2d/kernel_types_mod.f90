module kernel_types_mod

use constants_mod
use mesh_const_mod

integer, parameter :: N = N_NODES

type :: payload_hydro_t

    ! real(dp) :: states(-N+1:2*N,N_VARS)
    real(dp) :: state(N,N,N_VARS)
    real(dp) :: update(N,N,N_VARS)

end type

type :: auxiliary_hydro_t

    real(dp) :: halo(N_VARS,-N+1:2*N,-N+1:2*N)

    !real(dp) :: facevars(N_NODES,N_VARS,2,MESH_FACES)
    !real(dp) :: faceflux(N_NODES,N_VARS,  MESH_FACES)
    !real(dp) :: halovars(N_NODES,N_VARS,  MESH_FACES)
    !real(dp) :: entropy(N_NODES,N_NODES)

end type

end module
