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

subroutine DecideonESPF(Do_ESPF)

implicit real*8(a-h,o-z)
logical Do_ESPF

call Get_iScalar('System BitSwitch',iOption)

! this is a logical variable (logical function)
Do_ESPF = btest(iOption,11)

return

end subroutine DecideonESPF
