module setup_types_mod

use constants_mod

type :: auxiliary_setup_t

    integer :: signed_distance

    logical     :: inside(N_NODES,N_NODES)

    integer(8)  :: quadid(N_NODES,N_NODES)
    real(dp)    :: lagrvecs(N_NODES,N_DIMS)
 
end type

end module
