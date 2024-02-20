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

function distance(N,C1,C2)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: distance
integer(kind=iwp), intent(in) :: N
real(kind=wp), intent(in) :: C1(N), C2(N)
integer(kind=iwp) :: i
real(kind=wp) :: R, X

X = Zero
do i=1,N
  R = C1(i)-C2(i)
  X = X+R*R
end do
distance = sqrt(X)

return

end function distance
