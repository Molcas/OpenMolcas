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
subroutine Get_zArray(Label,data,nData)

implicit none
character*(*) Label
integer nData
real*8 RData(nData), IData(nData)
complex*16 data(nData)

call Get_dArray('R'//Label,RData,nData)
call Get_dArray('I'//Label,IData,nData)

data = RData+(0.0,1.0)*IData

return

end subroutine Get_zArray
