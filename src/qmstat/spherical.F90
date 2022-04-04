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

! Take higher multipole into spherical representation.
subroutine Spherical(dMul)

use Index_Functions, only: nTri_Elem1
use Constants, only: Three, Half
use Definitions, only: wp

implicit none
real(kind=wp), intent(inout) :: dMul(nTri_Elem1(2))
real(kind=wp) :: x2, xy, xz, y2, yz, z2
real(kind=wp), parameter :: d3 = sqrt(Three)

x2 = dMul(1)
y2 = dMul(4)
z2 = dMul(6)
xy = dMul(2)
xz = dMul(3)
yz = dMul(5)
dMul(1) = d3*xy
dMul(2) = d3*xz
dMul(3) = z2-Half*(x2+y2)
dMul(4) = d3*yz
dMul(5) = Half*d3*(x2-y2)

return

end subroutine Spherical
