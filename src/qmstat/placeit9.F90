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
#include "maxi.fh"
real(kind=wp) :: Coord(MxAt*3), Cordst(MxCen*MxPut,3)
integer(kind=iwp) :: info_atom(MxAt), iQ_Atoms
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
  CMSewx = CMSewx+Coord((i-1)*3+1)*info_atom(i)
  CMSewy = CMSewy+Coord((i-1)*3+2)*info_atom(i)
  CMSewz = CMSewz+Coord((i-1)*3+3)*info_atom(i)
  CMSamx = CMSamx+Cordst(i,1)*info_atom(i)
  CMSamy = CMSamy+Cordst(i,2)*info_atom(i)
  CMSamz = CMSamz+Cordst(i,3)*info_atom(i)
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
do i=1,iQ_Atoms
  Cordst(i,1) = Coord((i-1)*3+1)-Tx
  Cordst(i,2) = Coord((i-1)*3+2)-Ty
  Cordst(i,3) = Coord((i-1)*3+3)-Tz
end do

return

end subroutine PlaceIt9
