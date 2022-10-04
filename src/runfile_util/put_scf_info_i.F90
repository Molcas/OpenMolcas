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

subroutine Put_SCF_Info_I(iCase,Arr,nArr)

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: iCase, nArr, Arr(nArr)
character(len=24) :: Label

Label = 'SCFInfoI'
if (iCase == 1) Label = 'SCFInfoI_ab'
call Put_iArray(Label,Arr,nArr)

return

end subroutine Put_SCF_Info_I
