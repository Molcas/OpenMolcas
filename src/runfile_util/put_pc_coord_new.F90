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

subroutine Put_PC_Coord_New(Coord,nAtoms)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAtoms
real(kind=wp), intent(in) :: Coord(7,nAtoms)
character(len=*), parameter :: Label = 'GeoNewPC'

call Put_dArray(Label,Coord,7*nAtoms)

return

end subroutine Put_PC_Coord_New
