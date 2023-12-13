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
subroutine Get_zArray(Label,zData,nData)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Onei
use Definitions, only: wp, iwp

implicit none
character(len=*), intent(in) :: Label
integer(kind=iwp), intent(in) :: nData
complex(kind=wp), intent(out) :: zData(nData)
real(kind=wp), allocatable :: rData(:)

call mma_allocate(rData,nData,label='rData')
call Get_dArray('R'//Label,rData,nData)
zData(:) = rData
call Get_dArray('I'//Label,rData,nData)
zData(:) = zData+Onei*rData
call mma_deallocate(rData)

return

end subroutine Get_zArray
