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

subroutine DYAX(N,ALPHA,X,INCX,Y,INCY)
! MULTIPLY A VECTOR, X, BY A SCALAR AND STORE THE RESULT IN THE VECTOR Y.

#include "intent.fh"

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: N, INCX, INCY
real(kind=wp), intent(in) :: ALPHA, X(1+(N-1)*INCX)
real(kind=wp), intent(_OUT_) :: Y(1+(N-1)*INCY)
integer(kind=iwp) :: I

do I=1,N
  Y(1+(I-1)*INCY) = ALPHA*X(1+(I-1)*INCX)
end do

return

end subroutine DYAX
