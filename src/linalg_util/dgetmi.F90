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

subroutine DGETMI(A,ldA,N)

! TRANSPOSE A SQUARE MATRIX (IN-PLACE)

real*8 A(ldA,*), Temp

if (N <= 0) then
  write(6,*)
  write(6,*) '  *** Error in subroutine DGETMI ***'
  write(6,*) '  Invalid dimension of matrix A :'
  write(6,*) '  The number of rows/columns, N, must be larger than zero'
  write(6,*)
end if
if ((ldA <= 0) .or. (ldA < N)) then
  write(6,*)
  write(6,*) '  *** Error in subroutine DGETMI ***'
  write(6,*) '  Invalid leading dimension of matrix A :'
  write(6,*) '  ldA must be larger than 0 and larger than N'
  write(6,*)
end if

do i=1,N
  do j=1,i-1
    Temp = A(j,i)
    A(j,i) = A(i,j)
    A(i,j) = Temp
  end do
end do

return

end subroutine DGETMI
