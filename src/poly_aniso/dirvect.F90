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

subroutine dirvect(P1,R1,P2,R2,vec,dist)
! this Subroutine computes the directional vector of the origins of two points P1 and P2
! P1 : coords of the first point
! P2 : coords of the second point
! R1 : rot. matrix of the  first point to the general coordinate system
! R2 : rot. matrix of the second point to the general coordinate system

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: P1(3), R1(3,3), P2(3), R2(3,3)
real(kind=wp), intent(out) :: vec(3), dist
integer(kind=iwp) :: j
real(kind=wp) :: C1(3), C2(3)
real(kind=wp), external :: distance

C1(:) = Zero
C2(:) = Zero

do j=1,3
  C1(:) = C1(:)+P1(j)*R1(:,j)
  C2(:) = C2(:)+P2(j)*R2(:,j)
end do
dist = distance(3,C1,C2)
vec(:) = (C1(:)-C2(:))/dist

return

end subroutine dirvect
