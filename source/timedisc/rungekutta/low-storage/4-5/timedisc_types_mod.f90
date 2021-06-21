module timedisc_types_mod

use constants_mod

type :: auxiliary_timedisc_t

# if PP_N_DIMS == 2
    real(dp) ::  register0(N_NODES,N_NODES,N_VARS)

# elif PP_N_DIMS == 3
    real(dp) ::  register0(N_NODES,N_NODES,N_NODES,N_VARS)
# else
# error: unknown PP_N_DIMS
# endif

end type

end module
