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
integer(kind=iwp), intent(in) :: NDIM, MINMAX
real(kind=wp), intent(in) :: VECTOR(NDIM)
integer(kind=iwp) :: I
real(kind=wp) :: res

res = Zero
if (NDIM > 0) then

  if (MINMAX == 1) then
    res = huge(res)
    do I=1,NDIM
      res = min(res,abs(VECTOR(I)))
    end do
  else if (MINMAX == 2) then
    do I=1,NDIM
      res = max(res,abs(VECTOR(I)))
    end do
  end if

end if

FNDMNX = res

return

end function FNDMNX
