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
! Copyright (C) 1991, Roland Lindh                                     *
!***********************************************************************

function iDeg(Coor)
!***********************************************************************
!                                                                      *
! Object: to compute the degeneracy of a coordinate.                   *
!                                                                      *
!     Author: Roland Lindh, Dep. of Theoretical Chemistry,             *
!             University of Lund, SWEDEN                               *
!             March '91                                                *
!***********************************************************************

use Symmetry_Info, only: iOper, nIrrep
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: iDeg
real(kind=wp), intent(in) :: Coor(3)
integer(kind=iwp) :: i, j
real(kind=wp) :: Cx(3,8), r(3), x, y, z
logical(kind=iwp) :: New

iDeg = 1
Cx(:,1) = Coor(:)
do i=1,nIrrep-1
  r(:) = One
  if (btest(iOper(i),0)) r(1) = -One
  if (btest(iOper(i),1)) r(2) = -One
  if (btest(iOper(i),2)) r(3) = -One
  x = r(1)*Coor(1)
  y = r(2)*Coor(2)
  z = r(3)*Coor(3)
  New = .true.
  do j=1,iDeg
    if (New .and. (x == Cx(1,j)) .and. (y == Cx(2,j)) .and. (z == Cx(3,j))) New = .false.
  end do
  if (New) then
    iDeg = iDeg+1
    Cx(1,iDeg) = x
    Cx(2,iDeg) = y
    Cx(3,iDeg) = z
  end if
end do

return

end function iDeg
