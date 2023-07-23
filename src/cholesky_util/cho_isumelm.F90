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

function CHO_ISUMELM(IVEC,N)
!
! Purpose: sum elements of integer vector.

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: CHO_ISUMELM
integer(kind=iwp) :: IVEC(*), N
integer(Kind=iwp) :: I, ISUM

if (N > 0) then
  ISUM = IVEC(1)
  do I=2,N
    ISUM = ISUM+IVEC(I)
  end do
else
  ISUM = 0
end if

CHO_ISUMELM = ISUM

end function CHO_ISUMELM
