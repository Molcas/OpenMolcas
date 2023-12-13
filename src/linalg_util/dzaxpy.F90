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

subroutine DZAXPY(N,DA,DX,INCX,DY,INCY,DZ,INCZ)
! MULTIPLY A VECTOR, X, BY A SCALAR, ADD TO A VECTOR, Y, AND
! STORE THE RESULT IN THE VECTOR Z.

#include "intent.fh"

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: N, INCX, INCY, INCZ
real(kind=wp), intent(in) :: DA, DX(*), DY(*)
real(kind=wp), intent(_OUT_) :: DZ(*)
integer(kind=iwp) :: I, IX, IY, IZ, M, MP1

if (N <= 0) return
if ((INCX /= 1) .or. (INCY /= 1)) then

  ! CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL TO 1

  IY = 1
  IZ = 1
  if (INCY < 0) IY = (-N+1)*INCY+1
  if (INCZ < 0) IZ = (-N+1)*INCZ+1
  if (DA /= ZERO) then
    IX = 1
    if (INCX < 0) IX = (-N+1)*INCX+1
    do I=1,N
      DZ(IY) = DY(IY)+DA*DX(IX)
      IX = IX+INCX
      IY = IY+INCY
      IZ = IZ+INCZ
    end do
  else
    do I=1,N
      DZ(IY) = DY(IY)
      IY = IY+INCY
      IZ = IZ+INCZ
    end do
  end if

else

  ! CODE FOR BOTH INCREMENTS EQUAL TO 1

  ! CLEAN-UP LOOP

  M = mod(N,4)
  if (DA /= ZERO) then
    if (M /= 0) then
      do I=1,M
        DZ(I) = DY(I)+DA*DX(I)
      end do
      if (N < 4) return
    end if
    MP1 = M+1
    do I=MP1,N,4
      DZ(I) = DY(I)+DA*DX(I)
      DZ(I+1) = DY(I+1)+DA*DX(I+1)
      DZ(I+2) = DY(I+2)+DA*DX(I+2)
      DZ(I+3) = DY(I+3)+DA*DX(I+3)
    end do
  else
    if (M /= 0) then
      do I=1,M
        DZ(I) = DY(I)
      end do
      if (N < 4) return
    end if
    MP1 = M+1
    do I=MP1,N,4
      DZ(I) = DY(I)
      DZ(I+1) = DY(I+1)
      DZ(I+2) = DY(I+2)
      DZ(I+3) = DY(I+3)
    end do
  end if

end if

return

end subroutine DZAXPY
