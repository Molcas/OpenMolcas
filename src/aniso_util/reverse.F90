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

subroutine REVERSE(A_dir,A_inv,det)
! THIS ROUTINE CALCULATES THE INVERSE OF A SQUARE 3x3 MATRIX, AND ITS DETERMINANT.

implicit none
integer, parameter :: wp = kind(0.d0)
real(kind=8) :: A_dir(3,3)
real(kind=8) :: A_inv(3,3)
real(kind=8) :: A(3,3)
real(kind=8) :: B(3,3)
real(kind=8) :: det
real(kind=8) :: FindDetR
external :: FindDetR

det = 0.0_wp
call dcopy_(3*3,[0.0_wp],0,A,1)
call dcopy_(3*3,[0.0_wp],0,B,1)
call dcopy_(3*3,[0.0_wp],0,A_inv,1)
call dcopy_(3*3,A_dir,1,A,1)

det = FindDetR(A,3)

B(1,1) = A(2,2)*A(3,3)-A(3,2)*A(2,3)
B(1,2) = -A(1,2)*A(3,3)+A(3,2)*A(1,3)
B(1,3) = A(1,2)*A(2,3)-A(2,2)*A(1,3)
B(2,1) = -A(2,1)*A(3,3)+A(3,1)*A(2,3)
B(2,2) = A(1,1)*A(3,3)-A(3,1)*A(1,3)
B(2,3) = -A(1,1)*A(2,3)+A(2,1)*A(1,3)
B(3,1) = A(2,1)*A(3,2)-A(3,1)*A(2,2)
B(3,2) = -A(1,1)*A(3,2)+A(3,1)*A(1,2)
B(3,3) = A(1,1)*A(2,2)-A(2,1)*A(1,2)

call dscal_(3*3,1.0_wp/det,B,1)
call dcopy_(3*3,B,1,A_inv,1)

return

end subroutine REVERSE
