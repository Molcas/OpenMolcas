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

subroutine Put_PLMO(PLMO,nPLMO)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nPLMO
real(kind=wp) :: PLMO(nPLMO)
character(len=24) :: Label

Label = 'PLMO'
call Put_dArray(Label,PLMO,nPLMO)

return

end subroutine Put_PLMO
