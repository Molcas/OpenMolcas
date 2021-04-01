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

subroutine vec(threshold,u,j,k,iErr)

use ZMatConv_Mod, only: Coords
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: threshold
real(kind=wp), intent(out) :: u(3)
integer(kind=iwp), intent(in) :: j, k
integer(kind=iwp), intent(out) :: iErr
integer(kind=iwp) :: i
real(kind=wp) :: r(3), r2

iErr = 0
r2 = Zero

do i=1,3
  r(i) = Coords(i,j)-Coords(i,k)
  r2 = r2+r(i)*r(i)
end do
r2 = sqrt(r2)
if (r2 < threshold) then
  iErr = 1
  return
end if
do i=1,3
  u(i) = r(i)/r2
end do

return

end subroutine vec
