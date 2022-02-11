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

subroutine FMUL2(A,B,C,NROW,NCOL,N)

use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: NROW, NCOL, N
real(kind=wp) :: A(NROW,N), B(NCOL,N), C(NROW,NCOL)
#include "warnings.h"
integer(kind=iwp) :: I, J, J1, K
real(kind=wp) :: CJ(1000), FAC !IFG

if (nRow > 1000) then
  write(u6,*)
  call XFLUSH(u6)
  write(u6,*) ' *** Error in Subroutine FMUL2 ***'
  call XFLUSH(u6)
  write(u6,*) ' row dimension exceeds local buffer size'
  call XFLUSH(u6)
  write(u6,*)
  call XFLUSH(u6)
  call Quit(_RC_INTERNAL_ERROR_)
end if

do J=1,NCOL
  do I=1,NROW
    CJ(I) = Zero
  end do
  if (J == NCOL) GO TO 16
  J1 = J+1
  do K=1,N
    FAC = B(J,K)
    if (FAC == Zero) GO TO 20
    do I=J1,NROW
      CJ(I) = CJ(I)+FAC*A(I,K)
    end do
20  continue
  end do
16 do I=1,NROW
    C(I,J) = CJ(I)
  end do
end do

return

end subroutine FMUL2
