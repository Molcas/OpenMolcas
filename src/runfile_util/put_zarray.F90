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

! complex number in runfile
subroutine Put_zArray(Label,zData,nData)

use Definitions, only: wp, iwp

implicit none
character(len=*), intent(in) :: Label
integer(kind=iwp), intent(in) :: nData
complex(kind=wp), intent(in) :: zData(nData)

call Put_dArray('R'//Label,real(zData),nData)
call Put_dArray('I'//Label,aimag(zData),nData)

return

end subroutine Put_zArray
