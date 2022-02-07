!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1991, Roland Lindh                                     *
!***********************************************************************

subroutine CoW(Coor,CoF,W,nAtom,T)
!***********************************************************************
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             January '91                                              *
!***********************************************************************

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAtom
real(kind=wp), intent(in) :: Coor(3,nAtom), W(nAtom)
real(kind=wp), intent(out) :: CoF(3), T
#include "print.fh"
integer(kind=iwp) :: iAtom, iPrint, iRout

iRout = 140
iPrint = nPrint(iRout)

if (iPrint >= 99) then
  call RecPrt(' In CoW: Coor',' ',Coor,3,nAtom)
  call RecPrt(' In CoW: W',' ',W,nAtom,1)
end if
T = Zero
do iAtom=1,nAtom
  T = T+W(iAtom)
end do
CoF(:) = Zero
do iAtom=1,nAtom
  CoF(:) = CoF(:)+Coor(:,iAtom)*W(iAtom)
end do
if (T /= Zero) then
  CoF(:) = CoF(:)/T
else
  CoF(:) = Zero
end if
if (iPrint >= 99) then
  call RecPrt(' In CoW: CoF',' ',CoF,1,3)
  call RecPrt(' In CoW: T',' ',[T],1,1)
end if

return

end subroutine CoW
