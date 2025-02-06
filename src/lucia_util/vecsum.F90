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

subroutine VECSUM(C,A,B,FACA,FACB,NDIM)
! CALCULATE THE VECTOR C(I)=FACA*A(I)+FACB*B(I)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: C(*), A(*), B(*), FACA, FACB
integer(kind=iwp) :: NDIM
integer(kind=iwp) :: I
real(kind=wp) :: S

if ((FACA /= Zero) .and. (FACB /= Zero)) then
  do I=1,NDIM
    S = FACA*A(I)+FACB*B(I)
    C(I) = S
  end do

else if ((FACA == Zero) .and. (FACB /= Zero)) then
  do I=1,NDIM
    S = FACB*B(I)
    C(I) = S
  end do

else if ((FACA /= Zero) .and. (FACB == Zero)) then
  do I=1,NDIM
    S = FACA*A(I)
    C(I) = S
  end do

else if ((FACA == Zero) .and. (FACB == Zero)) then
  do I=1,NDIM
    C(I) = Zero
  end do

end if

end subroutine VECSUM
