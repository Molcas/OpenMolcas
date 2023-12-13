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

! Routine to give base vectors of the plane with v as normal.
subroutine PlaneVectors(u,w,v,Rinv)

use Constants, only: Zero, One, Half
use Definitions, only: wp

implicit none
real(kind=wp), intent(out) :: u(3), w(3)
real(kind=wp), intent(in) :: v(3), Rinv
real(kind=wp) :: const, dLu, p(3), Scal, Shitx, Shity, Shitz

! Construct an arbitrary normalized vector orthogonal to the v-vector.

const = Zero
Shitx = One
Shity = Zero
Shitz = Zero
do
  p(1) = Shitx+const
  p(2) = Shity+Half*const
  p(3) = Shitz-const
  Scal = p(1)*v(1)+p(2)*v(2)+p(3)*v(3)
  u(:) = p-Scal*Rinv**2*v
  if ((abs(u(1)) >= 1.0e-6_wp) .or. (abs(u(2)) >= 1.0e-6_wp) .or. (abs(u(3)) >= 1.0-6_wp)) exit
  const = const+One
end do
dLu = sqrt(u(1)**2+u(2)**2+u(3)**2)
u(:) = u/dLu

! Construct the final pi-vector, which is orthogonal to the v-vector
! and the recently constructed pi-vector.

w(1) = Rinv*(u(2)*v(3)-u(3)*v(2))
w(2) = Rinv*(u(3)*v(1)-u(1)*v(3))
w(3) = Rinv*(u(1)*v(2)-u(2)*v(1))

return

end subroutine PlaneVectors
