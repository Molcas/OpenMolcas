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

subroutine ICOPY(n,x,incX,y,incY)
!***********************************************************************
!                                                                      *
!     Copy vector X into vector Y                                      *
!                                                                      *
!     calling arguments:                                               *
!     N       : Integer, input.                                        *
!               Number of input elements.                              *
!     X       : Array of Integer, input.                               *
!               Input vector, X.                                       *
!     incX    : Integer, input.                                        *
!               Stride of vector X.                                    *
!     Y       : Array of Integer, output.                              *
!               Output vector, Y.                                      *
!     incY    : Integer, input.                                        *
!               Stride of vector Y.                                    *
!                                                                      *
!***********************************************************************

use Definitions, only: iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: n, x(*), incx, incy
integer(kind=iwp), intent(_OUT_) :: y(*)
integer(kind=iwp) :: i, ix, iy, m, mp1

! copies integer vector, x, to integer vector, y.
! uses unrolled loops for increments equal to one.
if (n <= 0) return
if ((incx == 1) .and. (incy == 1)) then

  ! code for both increments equal to 1

  ! clean-up loop

  m = mod(n,7)
  if (m /= 0) y(1:m) = x(1:m)
  if ((m == 0) .or. (n >= 7)) then
    mp1 = m+1
    do i=mp1,n,7
      y(i:i+6) = x(i:i+6)
    end do
  end if

else

  ! code for unequal increments or equal increments not equal to 1

  ix = 1
  iy = 1
  if (incx < 0) ix = (-n+1)*incx+1
  if (incy < 0) iy = (-n+1)*incy+1
  do i=1,n
    y(iy) = x(ix)
    ix = ix+incx
    iy = iy+incy
  end do

end if

return

end subroutine ICOPY
