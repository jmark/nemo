module constants_mod

! use ISO_FORTRAN_ENV

integer, parameter :: N_NODES   = PP_N_NODES
integer, parameter :: N_DIMS    = PP_N_DIMS
integer, parameter :: N_VARS    = PP_N_VARS

! integer, parameter :: i1 = INT8
! integer, parameter :: i2 = INT16
! integer, parameter :: i4 = INT32
! integer, parameter :: i8 = INT64

! integer, parameter :: sp = REAL32
! integer, parameter :: dp = REAL64
! integer, parameter :: qp = REAL128

! integer, parameter :: i4 = selected_int_kind(15)
integer, parameter :: SP = selected_real_kind(7)
integer, parameter :: dp = selected_real_kind(15)
integer, parameter :: qp = selected_real_kind(31)

real(kind=dp), parameter :: PI = 4.0_dp * atan(1.0_dp)
real(kind=dp), parameter :: TOL = 4.0e-16_dp
real(kind=dp), parameter :: EPSIL = 4.0e-16_dp

!! direction ID codes
integer, parameter :: X_DIR = 1
integer, parameter :: Y_DIR = 2
integer, parameter :: Z_DIR = 3
integer, parameter :: W_DIR = 4

!! face ID codes
integer, parameter :: ZEN = 0
integer, parameter :: NOR = 1
integer, parameter :: SOU = 2
integer, parameter :: WES = 3
integer, parameter :: EAS = 4
integer, parameter :: FRO = 5
integer, parameter :: BAC = 6

!! corner ID codes
integer, parameter :: NWF = 1
integer, parameter :: SWF = 2
integer, parameter :: NEF = 3
integer, parameter :: SEF = 4

integer, parameter :: NWB = 5
integer, parameter :: SWB = 6
integer, parameter :: NEB = 7
integer, parameter :: SEB = 8

end module
