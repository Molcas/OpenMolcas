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
subroutine Put_zArray(Label,data,nData)

implicit none
character*(*) Label
integer nData
complex*16 data(nData)

call Put_dArray('R'//Label,real(data),nData)
call Put_dArray('I'//Label,aimag(data),nData)

return

end subroutine Put_zArray
