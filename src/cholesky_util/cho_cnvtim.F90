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

subroutine CHO_CNVTIM(TIM,IHR,IMN,SEC)
!
! Purpose: convert TIM to hours/minutes/seconds

use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: TIM
integer(kind=iwp), intent(out) :: IHR, IMN
real(kind=wp), intent(out) :: SEC
real(kind=wp) :: XHR, XMN
real(kind=wp), parameter :: hour_s = 3600.0_wp, min_s = 60.0_wp

XHR = TIM/hour_s
IHR = int(XHR)

XMN = (TIM-hour_s*real(IHR,kind=wp))/min_s
IMN = int(XMN)

SEC = TIM-hour_s*real(IHR,kind=wp)-min_s*real(IMN,kind=wp)

end subroutine CHO_CNVTIM
