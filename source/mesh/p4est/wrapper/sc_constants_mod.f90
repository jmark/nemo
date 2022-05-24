module sc_constants_mod

integer, parameter :: SC_LP_DEFAULT    = -1     !! this selects the SC default threshold
integer, parameter :: SC_LP_ALWAYS     =  0     !! this will log everything
integer, parameter :: SC_LP_TRACE      =  1     !! this will prefix file and line number
integer, parameter :: SC_LP_DEBUG      =  2     !! any information on the internal state
integer, parameter :: SC_LP_VERBOSE    =  3     !! information on conditions, decisions
integer, parameter :: SC_LP_INFO       =  4     !! the main things a function is doing
integer, parameter :: SC_LP_STATISTICS =  5     !! important for consistency/performance
integer, parameter :: SC_LP_PRODUCTION =  6     !! a few lines for a major api function
integer, parameter :: SC_LP_ESSENTIAL  =  7     !! this logs a few lines max per program
integer, parameter :: SC_LP_ERROR      =  8     !! this logs errors only
integer, parameter :: SC_LP_SILENT     =  9     !! this never logs anything

end module
