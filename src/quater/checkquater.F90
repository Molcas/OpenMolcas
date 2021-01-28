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
!  CheckQuater
!
!> @brief
!>   Check whether the quaternion represents a rotation
!> @author Y. Carissan
!>
!> @details
!> Check whether the quaternion represents a rotation.
!>
!> @param[in] Q The quaternion to be checked
!***********************************************************************

subroutine CheckQuater(Q)

use Quater_globals, only: debug
use Constants, only: One, Two, pi
use Definitions, only: wp, r8, u6

implicit none
real(kind=wp), intent(in) :: Q(0:3)
real(kind=wp) :: res, angle, axis(3)
real(kind=wp), parameter :: thrs=1.0e-6_wp
real(kind=wp), external :: modangle
real(kind=r8), external :: ddot_

res = ddot_(4,Q,1,Q,1)

if (abs(res-One) > thrs) then
  call RecPrt('Quaternion tested',' ',Q,4,1)
  call SysAbendMsg('CheckQuater','Quaternion does not represent a rotation','')
end if

angle = modangle(Two*acos(Q(0)),Two*pi)

if (debug) then
  call RecPrt('Quaternion',' ',Q(0),4,1)
  write(u6,'(a8,f10.6,a3,f10.2,a3)') 'Angle = ',angle,'Rad',180*angle/pi,'Deg'
end if
axis(:) = Q(1:3)
call normalizeVec(axis)
if (debug) then
  call RecPrt('Axis',' ',axis,3,1)
end if

return

end subroutine CheckQuater
