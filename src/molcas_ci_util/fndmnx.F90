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

function FNDMNX(VECTOR,NDIM,MINMAX)
! FIND SMALLEST(MINMAX=1) OR LARGEST(MINMAX=2)
! ABSOLUTE VALUE OF ELEMENTS IN VECTOR

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: FNDMNX
real(kind=wp) :: VECTOR(*)
integer(kind=iwp) :: NDIM, MINMAX
integer(kind=iwp) :: I
real(kind=wp) :: res

! jwk-cleanup
res = Zero
if (MINMAX == 1) then
  res = abs(VECTOR(1))
  do I=2,NDIM
    res = min(res,abs(VECTOR(I)))
  end do
end if

if (MINMAX == 2) then
  res = abs(VECTOR(1))
  do I=2,NDIM
    res = max(res,abs(VECTOR(I)))
  end do
end if

FNDMNX = res

return

end function FNDMNX
