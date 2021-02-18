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
!  quater
!
!> @brief
!>   Driver for quater
!> @author Y. Carissan
!>
!> @details
!> Driver for quater.
!>
!> @param[out] ireturn return code
!***********************************************************************

subroutine quater(ireturn)

use Quater_globals, only: debug
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: ireturn
real(kind=wp) :: U1(3), U2(3), V1(3), V2(3), V1best(3), V2best(3), Q(0:3)

debug = .false.
call quaterinit()

call RdInput_Quater(U1,U2,V1,V2)

if (debug) then
  write(u6,*) 'Reference axis'
  call RecPrt('U1',' ',U1,3,1)
  call RecPrt('U2',' ',U2,3,1)
  write(u6,*) 'New axis'
  call RecPrt('V1',' ',V1,3,1)
  call RecPrt('V2',' ',V2,3,1)
end if

call QuaterSolve(U1,U2,V1,V2,Q)

! for test
call Add_Info('Quaternion',Q,4,8)

if (debug) then
  write(u6,*) 'Normalized Reference axis'
  call RecPrt('U1',' ',U1,3,1)
  call RecPrt('U2',' ',U2,3,1)
  write(u6,*) 'Normalized New axis'
  call RecPrt('V1',' ',V1,3,1)
  call RecPrt('V2',' ',V2,3,1)
  call QuaterRotation(Q,U1,V1best)
  call QuaterRotation(Q,U2,V2best)
  call RecPrt('Best V1',' ',V1best,3,1)
  call RecPrt('Best V2',' ',V2best,3,1)
end if

call GenerateGeoms(Q)

call QuaterReport(V1best,V2best,V1,V2)

call QuaterFinish()

ireturn = 0

return

end subroutine quater
