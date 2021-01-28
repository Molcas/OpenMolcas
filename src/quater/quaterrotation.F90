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
! Copyright (C) Yannick Carissan                                       *
!***********************************************************************
!  quaterRotation
!
!> @brief
!>   Performs the rotation of \p U in \p V via the quaternion \p Q
!> @author Y. Carissan
!>
!> @details
!> Performs the rotation of \p U in \p V via the quaternion \p Q.
!>
!> @param[in]  Q Quaternion used for the rotation
!> @param[in]  U Vector to be rotated
!> @param[out] V Vector rotated
!***********************************************************************

subroutine quaterRotation(Q,U,V)

use Constants, only: One,Two
use Definitions, only: wp, r8

implicit none
real(kind=wp), intent(in) :: Q(0:3), U(3)
real(kind=wp), intent(out) :: V(3)
real(kind=wp) :: T(3), C1, C2, C3
real(kind=r8), external :: ddot_

call CheckQuater(Q)
call Cross(Q(1:3),U,T)         ! T = Q x U

C1 = Two*Q(0)**2-One           ! C1 = 2 * Q(0)^2 - 1
C2 = Two*Q(0)                  ! C2 = 2 * Q(0)
C3 = Two*ddot_(3,Q(1:3),1,U,1) ! C3 = 2 * Q.U
V(:) = C1*U(:) - C2*T(:) + C3*Q(1:3)

return

end subroutine quaterRotation
