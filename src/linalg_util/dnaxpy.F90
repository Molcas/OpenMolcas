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

subroutine DNAXPY(N,M,A,INCA,X,INCXI,INCXO,Y,INCYI,INCYO)
! MULTIPLY A VECTOR, X, BY A SCALAR, ADD TO A VECTOR, Y, AND
! STORE THE RESULT IN THE VECTOR Y. REPEAT THIS N TIMES.

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: N, M, INCA, INCXI, INCXO, INCYI, INCYO
real(kind=wp), intent(in) :: A((N-1)*INCA+1), X(((M-1)*INCXI+1)*((N-1)*INCXO+1))
real(kind=wp), intent(inout) :: Y(((M-1)*INCYI+1)*((N-1)*INCYO+1))
integer(kind=iwp) :: I

do I=1,N
  call DAXPY_(M,A(1+(I-1)*INCA),X(1+(I-1)*INCXO),INCXI,Y(1+(I-1)*INCYO),INCYI)
end do

return

end subroutine DNAXPY
