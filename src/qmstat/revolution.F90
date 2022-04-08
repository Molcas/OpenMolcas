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

! Construct rotation matrix.
subroutine Revolution(v,Rinv,Rotte)

use Definitions, only: wp

implicit none
real(kind=wp), intent(in) :: v(3), Rinv
real(kind=wp), intent(out) :: Rotte(3,3)
real(kind=wp) :: u(3), w(3), t(3)

! Obtain base-vectors for the plane to which v in the normal vector.

call PlaneVectors(u,w,v,Rinv)

! Normalize v.

t(:) = Rinv*v

! Assemble rotation matrix

Rotte(1,:) = u
Rotte(2,:) = w
Rotte(3,:) = t

return

end subroutine Revolution
