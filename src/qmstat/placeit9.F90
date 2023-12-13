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

!----------------------------------------------------------------------*
! With this function we wish to place the QM-molecule properly when we *
! run with solvetn configurations from the sampfile. This we do by     *
! making the center-of-masses to coincide.                             *
!----------------------------------------------------------------------*
subroutine PlaceIt9(Coord,Cordst,info_atom,iQ_Atoms)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iQ_Atoms, info_atom(iQ_Atoms)
real(kind=wp), intent(in) :: Coord(3,iQ_Atoms)
real(kind=wp), intent(out) :: Cordst(3,iQ_Atoms)
integer(kind=iwp) :: i
real(kind=wp) :: CMSamx, CMSamy, CMSamz, CMSewx, CMSewy, CMSewz, Tx, Ty, Tz, Wtot

CMSewx = Zero
CMSewy = Zero
CMSewz = Zero
CMSamx = Zero
CMSamy = Zero
CMSamz = Zero
Wtot = Zero
do i=1,iQ_Atoms
  CMSewx = CMSewx+Coord(1,i)*info_atom(i)
  CMSewy = CMSewy+Coord(2,i)*info_atom(i)
  CMSewz = CMSewz+Coord(3,i)*info_atom(i)
  CMSamx = CMSamx+Cordst(1,i)*info_atom(i)
  CMSamy = CMSamy+Cordst(2,i)*info_atom(i)
  CMSamz = CMSamz+Cordst(3,i)*info_atom(i)
  Wtot = Wtot+real(info_atom(i),kind=wp)
end do
CMSewx = CMSewx/Wtot
CMSewy = CMSewy/Wtot
CMSewz = CMSewz/Wtot
CMSamx = CMSamx/Wtot
CMSamy = CMSamy/Wtot
CMSamz = CMSamz/Wtot
Tx = CMSewx-CMSamx
Ty = CMSewy-CMSamy
Tz = CMSewz-CMSamz
Cordst(1,:) = Coord(1,:)-Tx
Cordst(2,:) = Coord(2,:)-Ty
Cordst(3,:) = Coord(3,:)-Tz

return

end subroutine PlaceIt9
