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

subroutine Put_Coord_Full(Coord,nAtoms)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAtoms
real(kind=wp), intent(in) :: Coord(3,nAtoms)
integer(kind=iwp) :: nAtoms_All

call Get_nAtoms_All(nAtoms_All)
call Put_Coord_New(Coord,nAtoms_All)
call Put_dArray('MMO Coords',Coord(:,nAtoms_All+1:),3*(nAtoms-nAtoms_All))

return

end subroutine Put_Coord_Full
