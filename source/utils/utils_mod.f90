module utils_mod

use constants_mod

interface kronproduct
    procedure kronproduct_matmat
    procedure kronproduct_diagMat
end interface

interface pp
    procedure vp
    procedure mp
    procedure bp
    procedure tp
    procedure tp5
end interface

contains

!! pure elemental function logmean(uL,uR) result(r)
!!     real(dp), intent(in)    :: uL,uR
!!     real(dp)                :: r
!! 
!!     r = (uL-uR+TOL)/(log(uL)-log(uR)+TOL)
!! end function

pure function calc_dist(a,b) result(dist)
            
    real(dp), intent(in)    :: a(N_DIMS),b(N_DIMS)
    real(dp)                :: dist
    real(dp)                :: v(N_DIMS)

    v = b-a
    dist = SQRT(SUM(v*v))

end function

pure elemental function logmean(aL,aR) result(r)

    real, intent(in)         :: aL, aR
    real(dp)                     :: x,u,r
    real(dp),parameter           :: eps = 1.0e-02_dp

    x = al/ar
    u = (x*(x-2.0_dp)+1.0_dp)/(x*(x+2.0_dp)+1.0_dp)
    r = merge((aL+aR)*52.50_dp/(105.0_dp + u*(35.0_dp + u*(21.0_dp +u*15.0_dp))), (aL-aR)/log(x), u < eps)

end function

!! subroutine write_headers(outputfp, buflen, buffer)
!! 
!!     character(len=*), intent(in)    :: outputfp
!!     integer,intent(in)                 :: buflen
!!     character(len=*), intent(in)    :: buffer(:)
!! 
!!     integer :: ioUnit
!!     integer :: openStat
!! 
!!     character(range(buflen)+2) :: countStr
!! 
!!     OPEN(NEWUNIT  = ioUnit             , &
!!          FILE     = TRIM(outputfp)     , &
!!          FORM     = 'FORMATTED'        , &
!!          STATUS   = 'UNKNOWN'          , &
!!          POSITION = 'REWIND'           , &
!!          RECL     = 50000              , &
!!          IOSTAT   = openStat             )
!! 
!!     IF(openStat.NE.0) THEN
!!         CALL abort('ERROR: cannot open '// TRIM(outputfp))
!!     END IF
!! 
!!     !! write nr. of columns to string
!!     write(countStr,'(i0)') buflen+1
!!     WRITE(ioUnit, '('// trim(countStr) //'(A))') "# ", buffer(1:buflen)
!!     CLOSE(ioUnit)
!! 
!! end subroutine
!! 

! This is a simple function to search for an available unit.
! LUN_MIN and LUN_MAX define the range of possible LUNs to check.
! The UNIT value is returned by the function, and also by the optional
! argument. This allows the function to be used directly in an OPEN
! statement, and optionally save the result in a local variable.
! If no units are available, -1 is returned.
function utils_get_newunit(unit) result(newunit)

    integer, intent(out), optional :: unit

    integer, parameter :: LUN_MIN=10, LUN_MAX=1000
    logical :: opened
    integer :: lun

    integer :: newunit

    newunit=-1

    do lun=LUN_MIN,LUN_MAX
        inquire(unit=lun,opened=opened)
        if (.not. opened) then
            newunit=lun
            exit
        end if
    end do

    if (present(unit)) unit = newunit

end function

subroutine write2file(outputfp, buflen, buffer)

    character(len=*), intent(in)    :: outputfp
    integer,intent(in) :: buflen
    real(dp),intent(in) :: buffer(:)

    integer :: ioUnit
    integer :: openStat

    character(range(buflen)+2) :: countStr

    !OPEN(NEWUNIT  = ioUnit             , &
    OPEN(UNIT     = utils_get_newunit(ioUnit), &
         FILE     = TRIM(outputfp)     , &
         FORM     = 'FORMATTED'        , &
         STATUS   = 'UNKNOWN'          , &
         POSITION = 'APPEND'           , &
         RECL     = 50000              , &
         IOSTAT   = openStat             )

    IF(openStat.NE.0) THEN
        !CALL abort('ERROR: cannot open '// TRIM(outputfp))
        write (*,*) 'ERROR: cannot open '// TRIM(outputfp)
        return
    END IF

    !! write nr. of columns to string
    write(countStr,'(i0)') buflen
    WRITE(ioUnit, '('// trim(countStr) //'(ES32.23,1x))') buffer(1:buflen)
    CLOSE(ioUnit)

end subroutine write2file

function incr(n)

    implicit none

    integer, intent(inout) :: n
    integer :: incr

    n = n + 1
    incr = n

end function

pure function krondelta(i,j) result(k)

    integer, intent(in) :: i,j
    integer :: k

    if (i == j) then
        k = 1
    else
        k = 0
    end if

end function

pure function isintegerdivisible(a,b) result(c)

    integer, intent(in) :: a,b
    integer :: c

    if (mod(a,b) == 0) then
        c = 1
    else
        c = 0
    end if

end function

function kronproduct_matmat(A,B) result(C)

    real(dp), intent(in) :: A(:,:)
    real(dp), intent(in) :: B(:,:)

    real(dp) :: C(size(A,1)*size(B,1),size(A,2)*size(B,2))

    integer     :: i,j
    integer     :: bn(2),cn(2)

    bn = shape(B)
    cn = shape(A)*shape(B)

    do i = 1,cn(1); do j = 1,cn(2)
        C(i,j) = A((i-1)/bn(1)+1,(j-1)/bn(2)+1) * B(mod(i-1,bn(1))+1, mod(j-1,bn(2))+1)
    end do; end do

end function

function kronproduct_diagMat(A,B) result(C)

    real(dp), intent(in) :: A(:)
    real(dp), intent(in) :: B(:,:)

    real(dp) :: C(size(A)*size(B,1),size(A)*size(B,2))

    integer     :: i,j, iA, jA
    integer     :: bn(2),cn(2)

    bn = shape(B)
    cn = size(A)*shape(B)

    do i = 1,cn(1); do j = 1,cn(2)
        iA = (i-1)/bn(1)+1
        jA = (j-1)/bn(2)+1

        if (iA == jA) then
            C(i,j) = A(iA) * B(mod(i-1,bn(1))+1, mod(j-1,bn(2))+1)
        else
            C(i,j) = 0.0_dp
        end if
    end do; end do

end function

function eyeproduct(N,B) result(C)

    integer, intent(in)     :: N
    real(dp), intent(in)    :: B(:,:)

    real(dp) :: C(N*size(B,1),N*size(B,2))

    integer     :: i,j, iA, jA
    integer     :: bn(2),cn(2)

    bn = shape(B)
    cn = N*shape(B)

    C = 0.0_dp

    do i = 1,cn(1); do j = 1,cn(2)
        iA = (i-1)/bn(1)+1
        jA = (j-1)/bn(2)+1

        if (iA == jA) then
            C(i,j) = B(mod(i-1,bn(1))+1, mod(j-1,bn(2))+1)
        end if
    end do; end do

end function

pure function calc_index_window(nelems,nranks,rankid) result(window)

    integer(kind=8), intent(in) :: nelems
    integer(kind=4), intent(in) :: nranks,rankid
    integer(kind=8)             :: window(2)

    integer(kind=8) :: a,b,NelemsPerRank,NelemsResidual

    NelemsPerRank = nelems / nranks
    NelemsResidual = mod(nelems,INT(nranks,kind=8))

    if (rankid < NelemsResidual) then
        a = rankid * NelemsPerRank + mod(INT(rankid,kind=8),NelemsResidual) + 1
        b = a + NelemsPerRank
    else
        a = rankid * NelemsPerRank + NelemsResidual + 1
        b = a + NelemsPerRank - 1
    end if

    window = (/a,b/)

end function

pure elemental function vmin(a,b) result(c)

    real(dp), intent(in)    :: a,b
    real(dp)                :: c

    c = min(a,b)

end function

pure elemental function vmax(a,b) result(c)

    real(dp), intent(in)    :: a,b
    real(dp)                :: c

    c = max(a,b)

end function

function tostr(num) result(str)

    integer, intent(in) :: num
    character(len=16)   :: str

    write(str,'(i0)') num

end function

subroutine tp5(title,blk)

    character(len=*), intent(in)    :: title
    real(dp), intent(in)            :: blk(:,:,:,:,:)

    integer :: i,k,l,q

    do q = 1,size(blk,5)
        do l = 1,size(blk,4)
            do k = 1,size(blk,3)
                write (*,'((a12),(a6),(i3),(a6),(i3),(a6),(i3))') trim(title),' | q = ', q, ' | w = ', l, ' | z = ', k
                do i = 1,size(blk,1)
                    write (*,'(99(ES12.4))') blk(i,:,k,l,q)
                end do
                write (*,*)
            end do
            write (*,*)
        end do
        write (*,*)
    end do
    write (*,*)

end subroutine

subroutine tp(title,blk)

    character(len=*), intent(in)    :: title
    real(dp), intent(in)            :: blk(:,:,:,:)

    integer :: i,k,l

    do l = 1,size(blk,4)
        do k = 1,size(blk,3)
            write (*,'((a12),(a6),(i3),(a6),(i3))') trim(title),' | w = ', l, ' | z = ', k
            do i = 1,size(blk,1)
                write (*,'(99(ES12.4))') blk(i,:,k,l)
            end do
            write (*,*)
        end do
        write (*,*)
    end do
    write (*,*)

end subroutine

subroutine bp(title,blk)

    character(len=*), intent(in)    :: title
    real(dp), intent(in)            :: blk(:,:,:)

    integer :: i,k

    do k = 1,size(blk,3)
        write (*,*) trim(title) // ' | z = ' // tostr(k)
        do i = 1,size(blk,1)
            write (*,'(99(ES12.4))') blk(i,:,k)
        end do
        write (*,*)
    end do
    write (*,*)

end subroutine

subroutine mp(title,mat)

    character(len=*), intent(in)    :: title
    real(dp), intent(in)            :: mat(:,:)

    integer :: i

    write (*,*) trim(title)
    do i = 1,size(mat,1)
        write (*,'(99(ES12.4))') mat(i,:)
    end do
    write (*,*)

end subroutine

subroutine vp(title,vec)

    character(len=*), intent(in)    :: title
    real(dp), intent(in)            :: vec(:)

    integer :: i

    write (*,*) trim(title)
    write (*,'(99(ES12.4))') vec
    write (*,*)

end subroutine

end module
