module tensor_mod

use constants_mod, only: dp

interface blockmulXYZ
    procedure blockmulXYZ_single_op
    procedure blockmulXYZ_triple_op
end interface

interface blockmulXY
    procedure blockmulXY_single_op
    procedure blockmulXY_double_op
end interface

interface blockmulXZ
    procedure blockmulXZ_single_op
    procedure blockmulXZ_double_op
end interface

interface blockmulYZ
    procedure blockmulYZ_single_op
    procedure blockmulYZ_double_op
end interface

contains

pure function blockdotX(vect,cube) result(prod)

    real(dp), intent(in) :: vect(:)
    real(dp), intent(in) :: cube(:,:,:)

    real(dp) :: prod(size(cube,2),size(cube,3))

    integer :: i,j,q

    do i = 1,size(cube,2); do j = 1,size(cube,3);
        prod(i,j) = sum(vect*cube(:,i,j))
    end do; end do;
    
end function

pure function blockdotY(vect,cube) result(prod)

    real(dp), intent(in) :: vect(:)
    real(dp), intent(in) :: cube(:,:,:)

    real(dp) :: prod(size(cube,1),size(cube,3))

    integer :: i,j,q

    do i = 1,size(cube,1); do j = 1,size(cube,3);
        prod(i,j) = sum(vect*cube(i,:,j))
    end do; end do;
    
end function

pure function blockdotZ(vect,cube) result(prod)

    real(dp), intent(in) :: vect(:)
    real(dp), intent(in) :: cube(:,:,:)

    real(dp) :: prod(size(cube,1),size(cube,2))

    integer :: i,j,q

    do i = 1,size(cube,1); do j = 1,size(cube,2);
        prod(i,j) = sum(vect*cube(i,j,:))
    end do; end do;
    
end function

!! ========================================================================== !!

pure function blockmulX(matrix,incube) result(excube)

    real(dp), intent(in) :: matrix(:,:)
    real(dp), intent(in) :: incube(:,:,:)

    real(dp) :: excube(size(matrix,1),size(incube,2),size(incube,3))

    integer :: i,j,q

    do q = 1,size(matrix,1); do i = 1,size(incube,2); do j = 1,size(incube,3);
        excube(q,i,j) = sum(matrix(q,:) * incube(:,i,j))
    end do; end do; end do;

end function

pure function blockmulY(matrix,incube) result(excube)

    real(dp), intent(in) :: matrix(:,:)
    real(dp), intent(in) :: incube(:,:,:)

    real(dp) :: excube(size(incube,1),size(matrix,1),size(incube,3))

    integer :: i,j,q

    do q = 1,size(matrix,1); do i = 1,size(incube,1); do j = 1,size(incube,3);
        excube(i,q,j) = sum(matrix(q,:) * incube(i,:,j))
    end do; end do; end do;

end function

pure function blockmulZ(matrix,incube) result(excube)

    real(dp), intent(in) :: matrix(:,:)
    real(dp), intent(in) :: incube(:,:,:)

    real(dp) :: excube(size(incube,1),size(incube,2),size(matrix,1))

    integer :: i,j,q

    do q = 1,size(matrix,1); do i = 1,size(incube,1); do j = 1,size(incube,2);
        excube(i,j,q) = sum(matrix(q,:) * incube(i,j,:))
    end do; end do; end do;

end function

!! ========================================================================== !!

pure function blockmulXYZ_single_op(matrix,incube) result(excube)

    real(dp), intent(in) :: matrix(:,:)
    real(dp), intent(in) :: incube(:,:,:)

    real(dp) :: excube(size(matrix,1),size(matrix,1),size(matrix,1))

    excube = blockmulXYZ(matrix,matrix,matrix,incube) 

end function

pure function blockmulXYZ_triple_op(matrA,matrB,matrC,cube0) result(cubeC)

    real(dp), intent(in) :: matrA(:,:)
    real(dp), intent(in) :: matrB(:,:)
    real(dp), intent(in) :: matrC(:,:)
    real(dp), intent(in) :: cube0(:,:,:)

    real(dp) :: cubeA(size(matrA,1),size(cube0,2),size(cube0,3))
    real(dp) :: cubeB(size(matrA,1),size(matrB,1),size(cube0,3))
    real(dp) :: cubeC(size(matrA,1),size(matrB,1),size(matrC,1))

    cubeA = blockmulX(matrA,cube0) 
    cubeB = blockmulY(matrB,cubeA) 
    cubeC = blockmulZ(matrC,cubeB) 

end function

!! ========================================================================== !!

pure function blockmulXY_double_op(matrA,matrB,cube0) result(cubeB)

    real(dp), intent(in) :: matrA(:,:)
    real(dp), intent(in) :: matrB(:,:)
    real(dp), intent(in) :: cube0(:,:,:)

    real(dp) :: cubeA(size(matrA,1),size(cube0,2),size(cube0,3))
    real(dp) :: cubeB(size(matrA,1),size(matrB,1),size(cube0,3))

    cubeA = blockmulX(matrA,cube0) 
    cubeB = blockmulY(matrB,cubeA) 

end function

pure function blockmulXZ_double_op(matrA,matrB,cube0) result(cubeB)

    real(dp), intent(in) :: matrA(:,:)
    real(dp), intent(in) :: matrB(:,:)
    real(dp), intent(in) :: cube0(:,:,:)

    real(dp) :: cubeA(size(matrA,1),size(cube0,2),size(cube0,3))
    real(dp) :: cubeB(size(matrA,1),size(cube0,2),size(matrB,1))

    cubeA = blockmulX(matrA,cube0) 
    cubeB = blockmulZ(matrB,cubeA) 

end function

pure function blockmulYZ_double_op(matrA,matrB,cube0) result(cubeB)

    real(dp), intent(in) :: matrA(:,:)
    real(dp), intent(in) :: matrB(:,:)
    real(dp), intent(in) :: cube0(:,:,:)

    real(dp) :: cubeA(size(cube0,1),size(matrA,1),size(cube0,3))
    real(dp) :: cubeB(size(cube0,1),size(matrA,1),size(matrB,1))

    cubeA = blockmulY(matrA,cube0) 
    cubeB = blockmulZ(matrB,cubeA) 

end function

!! ========================================================================== !!

pure function blockmulXY_single_op(matrA,cube0) result(cubeA)

    real(dp), intent(in) :: matrA(:,:)
    real(dp), intent(in) :: cube0(:,:,:)

    real(dp) :: cubeA(size(matrA,1),size(matrA,1),size(cube0,3))

    cubeA = blockmulXY_double_op(matrA,matrA,cube0) 

end function

pure function blockmulXZ_single_op(matrA,cube0) result(cubeA)

    real(dp), intent(in) :: matrA(:,:)
    real(dp), intent(in) :: cube0(:,:,:)

    real(dp) :: cubeA(size(matrA,1),size(cube0,2),size(matrA,1))

    cubeA = blockmulXZ_double_op(matrA,matrA,cube0) 

end function

pure function blockmulYZ_single_op(matrA,cube0) result(cubeA)

    real(dp), intent(in) :: matrA(:,:)
    real(dp), intent(in) :: cube0(:,:,:)

    real(dp) :: cubeA(size(cube0,1),size(matrA,1),size(matrA,1))

    cubeA = blockmulYZ_double_op(matrA,matrA,cube0) 

end function

end module
