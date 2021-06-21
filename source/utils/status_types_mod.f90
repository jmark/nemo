module status_types_mod

integer, parameter :: STATUS_OK     = 0
integer, parameter :: STATUS_ERROR  = 1
integer, parameter :: STARTUP_ERROR = -1

type status_t

    logical                 :: ok       = .true.
    integer                 :: errcode  = STATUS_OK
    character(len=4096)     :: message  = ''
    
end type

end module
