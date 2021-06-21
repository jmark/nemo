module strings_mod

contains

function upper(s1)  result (s2)
    character(*)       :: s1
    character(len(s1)) :: s2
    character          :: ch
    integer,parameter  :: duc = iachar('A') - iachar('a')
    integer            :: i

    do i = 1,len(s1)
       ch = s1(i:i)
       if (ch >= 'a'.and.ch <= 'z') ch = char(iachar(ch)+duc)
       s2(i:i) = ch
    end do
end function

function lower(s1)  result (s2)
    character(*)       :: s1
    character(len(s1)) :: s2
    character          :: ch
    integer,parameter  :: duc = iachar('A') - iachar('a')
    integer            :: i

    do i = 1,len(s1)
       ch = s1(i:i)
       if (ch >= 'A'.and.ch <= 'Z') ch = char(iachar(ch)-duc)
       s2(i:i) = ch
    end do
end function

function trim2(s1)  result (s2)
    character(*)       :: s1
    character(len(s1)) :: s2

    s2 = trim(adjustL(s1))
end function

function strcomp(sA, sB) result(r)

    character(len=*), intent(in) :: sA
    character(len=*), intent(in) :: sB

    logical :: r
    integer :: i

    r = .true.
    do i = 1,MERGE(len(sA),len(sB),len(sA) < len(sB))     
        if (sA(i:i) /= sB(i:i)) then
            r = .false.
            exit
        end if
    end do

end function

function mkrankstr(irank) result(rankstr)
    integer, intent(in) :: irank
    character(len=10)   :: rankstr

    write(rankstr,'(I0.4)') irank
end function

end module
