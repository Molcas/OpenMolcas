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

! This subroutine should be in a module, to avoid explicit interfaces
#ifndef _IN_MODULE_
#error "This file must be compiled inside a module"
#endif

subroutine Get_PC_Coord_New(Coord,nData)

use stdalloc, only: mma_allocate
use Definitions, only: wp, iwp

implicit none
real(kind=wp), allocatable, intent(out) :: Coord(:)
integer(kind=iwp), intent(out) :: nData
logical(kind=iwp) :: Found
character(len=*), parameter :: Label = 'GeoNewPC'

call qpg_dArray(Label,Found,nData)
if ((.not. Found) .or. (nData == 0)) return
call mma_allocate(Coord,nData,label='Coord')
call Get_dArray(Label,Coord,nData)

return

end subroutine Get_PC_Coord_New
