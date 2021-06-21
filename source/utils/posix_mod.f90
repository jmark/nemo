module posix_mod

use iso_c_binding, only: c_int, c_char, c_int16_t

interface
    function c_mkdir(path,mode) bind(c,name="mkdir")
        use iso_c_binding, only: c_int, c_char, c_int16_t
        character(kind=c_char)      :: path(*)
        integer(c_int16_t), value   :: mode
        integer(c_int)              :: c_mkdir
    end function
end interface

!! int rename(const char *old_filename, const char *new_filename)
interface
    function c_rename(old_filename, new_filename) bind(c,name="rename")
        use iso_c_binding, only: c_int, c_char
        character(kind=c_char)      :: old_filename(*)
        character(kind=c_char)      :: new_filename(*)
        integer(c_int)              :: c_rename
    end function
end interface

contains

function posix_mkdir(path)
    use iso_c_binding, only: c_int16_t !! pmake: ignore
    character(len=*) :: path
    integer          :: posix_mkdir

    posix_mkdir = c_mkdir(trim(path) // char(0), int(o'750',c_int16_t))
end function

function posix_rename(old_filename,new_filename)
    character(len=*) :: old_filename
    character(len=*) :: new_filename
    integer          :: posix_rename

    posix_rename = c_rename(trim(old_filename) // char(0), trim(new_filename) // char(0))
end function

end module
