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

function Langevin_on()

use Definitions, only: iwp

implicit none
logical(kind=iwp) :: Langevin_on
integer(kind=iwp) :: iOption

!call Get_iOption(iOption)
call Get_iScalar('System BitSwitch',iOption)

Langevin_on = btest(iOption,3)

return

end function Langevin_on
