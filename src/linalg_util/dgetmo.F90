!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine DGETMO(A,ldA,M,N,B,ldB)

! TRANSPOSE A REGULAR MATRIX (OUT-OF-PLACE)

real*8 A(ldA,*), B(ldB,*)

if (M <= 0) then
  write(6,*)
  write(6,*) '  *** Error in subroutine DGETMO ***'
  write(6,*) '  Invalid dimension of matrix A :'
  write(6,*) '  The number of columns, M, must be larger than zero'
  write(6,*)
end if
if (N <= 0) then
  write(6,*)
  write(6,*) '  *** Error in subroutine DGETMO ***'
  write(6,*) '  Invalid leading dimension of matrix B :'
  write(6,*) '  The number of rows, N, must be larger than zero'
  write(6,*)
end if
if ((ldA <= 0) .or. (ldA < M)) then
  write(6,*)
  write(6,*) '  *** Error in subroutine DGETMO ***'
  write(6,*) '  Invalid leading dimension of matrix A :'
  write(6,*) '  ldA must be larger than 0 and larger than M'
  write(6,*)
end if
if ((ldB <= 0) .or. (ldB < N)) then
  write(6,*)
  write(6,*) '  *** Error in subroutine DGETMO ***'
  write(6,*) '  Invalid leading dimension of matrix B :'
  write(6,*) '  ldB must be larger than 0 and larger than N'
  write(6,*)
end if
INC = 8
do j=1,M,INC
  jj = min(M-j+1,INC)
  Go To(1,2,3,4,5,6,7,8) jj
  write(6,*) 'Error in DGETMO!'
1 continue
  do i=1,N
    B(i,j) = A(j,i)
  end do
  Go To 99
2 continue
  do i=1,N
    B(i,j) = A(j,i)
    B(i,j+1) = A(j+1,i)
  end do
  Go To 99
3 continue
  do i=1,N
    B(i,j) = A(j,i)
    B(i,j+1) = A(j+1,i)
    B(i,j+2) = A(j+2,i)
  end do
  Go To 99
4 continue
  do i=1,N
    B(i,j) = A(j,i)
    B(i,j+1) = A(j+1,i)
    B(i,j+2) = A(j+2,i)
    B(i,j+3) = A(j+3,i)
  end do
  Go To 99
5 continue
  do i=1,N
    B(i,j) = A(j,i)
    B(i,j+1) = A(j+1,i)
    B(i,j+2) = A(j+2,i)
    B(i,j+3) = A(j+3,i)
    B(i,j+4) = A(j+4,i)
  end do
  Go To 99
6 continue
  do i=1,N
    B(i,j) = A(j,i)
    B(i,j+1) = A(j+1,i)
    B(i,j+2) = A(j+2,i)
    B(i,j+3) = A(j+3,i)
    B(i,j+4) = A(j+4,i)
    B(i,j+5) = A(j+5,i)
  end do
  Go To 99
7 continue
  do i=1,N
    B(i,j) = A(j,i)
    B(i,j+1) = A(j+1,i)
    B(i,j+2) = A(j+2,i)
    B(i,j+3) = A(j+3,i)
    B(i,j+4) = A(j+4,i)
    B(i,j+5) = A(j+5,i)
    B(i,j+6) = A(j+6,i)
  end do
  Go To 99
8 continue
  do i=1,N
    B(i,j) = A(j,i)
    B(i,j+1) = A(j+1,i)
    B(i,j+2) = A(j+2,i)
    B(i,j+3) = A(j+3,i)
    B(i,j+4) = A(j+4,i)
    B(i,j+5) = A(j+5,i)
    B(i,j+6) = A(j+6,i)
    B(i,j+7) = A(j+7,i)
  end do
99 continue
end do

return

end subroutine DGETMO
