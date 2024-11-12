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

use Definitions, only: wp

implicit none
real(kind=wp), intent(in) :: A_dir(3,3)
real(kind=wp), intent(out) :: A_inv(3,3), det
real(kind=wp) :: A(3,3)
real(kind=wp), external :: FindDetR

A(:,:) = A_dir(:,:)

det = FindDetR(A,3)

A_inv(1,1) = A(2,2)*A(3,3)-A(3,2)*A(2,3)
A_inv(1,2) = -A(1,2)*A(3,3)+A(3,2)*A(1,3)
A_inv(1,3) = A(1,2)*A(2,3)-A(2,2)*A(1,3)
A_inv(2,1) = -A(2,1)*A(3,3)+A(3,1)*A(2,3)
A_inv(2,2) = A(1,1)*A(3,3)-A(3,1)*A(1,3)
A_inv(2,3) = -A(1,1)*A(2,3)+A(2,1)*A(1,3)
A_inv(3,1) = A(2,1)*A(3,2)-A(3,1)*A(2,2)
A_inv(3,2) = -A(1,1)*A(3,2)+A(3,1)*A(1,2)
A_inv(3,3) = A(1,1)*A(2,2)-A(2,1)*A(1,2)

A_inv(:,:) = A_inv(:,:)/det

return

end subroutine REVERSE
