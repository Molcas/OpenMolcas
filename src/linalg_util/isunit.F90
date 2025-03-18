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

function isUnit(A,N,ldA,thrs)

! Find out if a block of A is close enough to a unit matrix

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
logical(kind=iwp) :: isunit
integer(kind=iwp), intent(in) :: N, ldA
real(kind=wp), intent(in) :: A(ldA,N), thrs
integer(kind=iwp) :: i

do i=1,N
  if (abs(A(i,i)-One) > thrs) exit
  if (any(abs(A(1:i-1,i)) > thrs)) exit
  if (any(abs(A(i+1:N,i)) > thrs)) exit
end do
isunit = (i > N)

return

end function isUnit
