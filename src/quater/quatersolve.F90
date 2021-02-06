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
!  quatersolve
!
!> @brief
!>   Computes the quaternion which corresponds to the best rotation that transforms \p U1 in \p U2 and \p V1 in \p V2
!> @author Y. Carissan
!>
!> @details
!> Input vectors will be normalized and \p V2 changed.
!> Computes the quaternion which corresponds to the
!> rotation that transforms \p U1 in \p U2 and \p V1 in \p V2.
!> If such a rotation does not exists, i.e. if
!> \f[ U_1 \cdot V_1 \ne U_2 \cdot V_2 \f]
!> the quaternion contains the rotation
!> which transforms \p U1 into \p U2 and \p V1 into \p V2' such that
!> \f[ U_1 \cdot V_1 = U_2 \cdot V_2' \f]
!> and
!> \f[ U_1 \times V_1 = b (U_2 \times V_2') \f] with \f$ b \f$ real.
!>
!> @param[in,out] U1 vector \f$ U_1 \f$
!> @param[in,out] U2 vector \f$ U_2 \f$
!> @param[in,out] V1 vector \f$ V_1 \f$
!> @param[in,out] V2 vector \f$ V_2 \f$
!> @param[out]    Q  quaternion \f$ Q \f$
!***********************************************************************

subroutine quatersolve(U1,U2,V1,V2,Q)

use Quater_globals, only: debug
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, r8

implicit none
real(kind=wp), intent(inout) :: U1(3), U2(3), V1(3), V2(3)
real(kind=wp), intent(out) :: Q(0:3)
real(kind=wp) :: U3(3), V3(3), Uref(3), Vref(3), K(3), Vtmp(3), c
logical(kind=iwp) :: skip
real(kind=wp), parameter :: thrs = 1.0e-3_wp
real(kind=r8), external :: ddot_

if (debug) then
  call RecPrt('IN SOLVE U1',' ',U1,3,1)
  call RecPrt('IN SOLVE V1',' ',V1,3,1)
  call RecPrt('IN SOLVE U2',' ',U2,3,1)
  call RecPrt('IN SOLVE V2',' ',V2,3,1)
end if
call quatersetup(U1,U2,V1,V2)
if (debug) call RecPrt('new V2',' ',V2,3,1)

Uref(:) = U1(:)
Vref(:) = V1(:)

call GetKandC(U1,U2,V1,V2,K,C)

skip = .false.
if (C < thrs) then
  call cross(U1,U2,U3)
  call cross(V1,V2,V3)
  call getKandC(U1,V1,U3,V3,K,C)
  if (C < thrs) then
    call getKandC(U2,V2,U3,V3,K,C)
    if (C < thrs) then
      skip = .true.
    else
      Uref(:) = U2(:)
      Vref(:) = V2(:)
    end if
  end if
end if

if (skip) then
  Q(:) = [One,Zero,Zero,Zero]
else
  Q(1:3) = Half*K(:)/sqrt(C)

  call Cross(Uref,Q(1:3),Vtmp) ! Vtmp = Uref x Q

  Q(0) = Half*ddot_(3,Vref,1,Vtmp,1)/ddot_(3,Vtmp,1,Vtmp,1) ! Q(0)=0.5 * Vref.(UrefxQ)/(UrefxQ)^2
end if

call CheckQuater(Q)
call setMatrix(Q)
if (debug) call RecPrt('Quaternion',' ',Q,4,1)

return

end subroutine quatersolve
