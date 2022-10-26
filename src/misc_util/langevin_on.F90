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

logical function Langevin_on()

implicit real*8(a-h,o-z)

!call Get_iOption(iOption)
call Get_iScalar('System BitSwitch',iOption)

Langevin_on = iand(iOption,8) == 8

return

end function Langevin_on
