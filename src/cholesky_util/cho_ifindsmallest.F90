!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2007, Thomas Bondo Pedersen                            *
!***********************************************************************

function Cho_iFindSmallest(iVec,n)
!
! Thomas Bondo Pedersen, June 2007.
!
! Return index of smallest element in array iVec.

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: Cho_iFindSmallest
integer(kind=iwp), intent(in) :: n, iVec(n)
integer(kind=iwp) :: i, imin

if (n < 1) then
  imin = 0
else
  imin = n
  do i=n-1,1,-1
    if (iVec(i) < iVec(imin)) imin = i
  end do
end if

Cho_iFindSmallest = imin

end function Cho_iFindSmallest
