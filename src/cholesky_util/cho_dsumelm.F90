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

function CHO_DSUMELM(VEC,N)
!
! Purpose: sum elements of double precision vector.

implicit none
real*8 CHO_DSUMELM
real*8 VEC(*), DSUM
integer I, N

if (N > 0) then
  DSUM = VEC(1)
  do I=2,N
    DSUM = DSUM+VEC(I)
  end do
else
  DSUM = 0.0d0
end if

CHO_DSUMELM = DSUM

end function CHO_DSUMELM
