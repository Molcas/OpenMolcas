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

function RF_On()

use Definitions, only: iwp

implicit none
logical(kind=iwp) :: RF_On
integer(kind=iwp) :: iOption = -99

!call Get_iOption(iOption)
if (iOption == -99) call Get_iScalar('System BitSwitch',iOption)

RF_On = btest(iOption,2)

return

end function RF_On
