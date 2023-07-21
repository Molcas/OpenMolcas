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

implicit none
real*8 TIM
integer IHR, IMN
real*8 XHR, XMN, SEC

XHR = TIM/3.6d3
IHR = int(XHR)

XMN = (TIM-3.6d3*dble(IHR))/6.0d1
IMN = int(XMN)

SEC = TIM-3.6d3*dble(IHR)-6.0d1*dble(IMN)

end subroutine CHO_CNVTIM
