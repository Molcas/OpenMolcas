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

implicit real*8(a-h,o-z)
dimension u(3), w(3), v(3), p(3)

! Construct an arbitrary normalized vector orthogonal to the v-vector.

const = 0.0d0
Shitx = 1.0d0
Shity = 0.0d0
Shitz = 0.0d0
do
  p(1) = Shitx+1.0d0*const
  p(2) = Shity+0.5d0*const
  p(3) = Shitz-1.0d0*const
  Scal = p(1)*v(1)+p(2)*v(2)+p(3)*v(3)
  u(1) = p(1)-Scal*Rinv**2*v(1)
  u(2) = p(2)-Scal*Rinv**2*v(2)
  u(3) = p(3)-Scal*Rinv**2*v(3)
  if ((abs(u(1)) >= 1d-6) .or. (abs(u(2)) >= 1d-6) .or. (abs(u(3)) >= 1d-6)) exit
  const = const+1.0d0
end do
dLu = sqrt(u(1)**2+u(2)**2+u(3)**2)
u(1) = u(1)/dLu
u(2) = u(2)/dLu
u(3) = u(3)/dLu

! Construct the final pi-vector, which is orthogonal to the v-vector
! and the recently constructed pi-vector.

w(1) = Rinv*(u(2)*v(3)-u(3)*v(2))
w(2) = Rinv*(u(3)*v(1)-u(1)*v(3))
w(3) = Rinv*(u(1)*v(2)-u(2)*v(1))

return

end subroutine PlaneVectors
