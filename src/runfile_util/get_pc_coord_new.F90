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

subroutine Get_PC_Coord_New(Coord,nData)

implicit real*8(a-h,o-z)
#include "stdalloc.fh"
character*24 Label
logical Found
real*8, dimension(:), allocatable :: Coord

Label = 'GeoNewPC'
call qpg_dArray(Label,Found,nData)
if ((.not. Found) .or. (nData == 0)) return
call mma_allocate(Coord,nData,label='Coord')
call Get_dArray(Label,Coord,nData)

return

end subroutine Get_PC_Coord_New
