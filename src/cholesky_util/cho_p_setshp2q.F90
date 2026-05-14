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

subroutine Cho_P_SetShP2Q(irc,iLoc,iShlAB,nAB)
!
! Purpose: set mapping from shell pair to qualified.
!          Global reduced set index arrays are needed, so we swap
!          local and global index arrays before (and after) calling
!          the original serial routine.

use Cholesky, only: Cho_Real_Par
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp), intent(in) :: iLoc, iShlAB, nAB(8)

if (Cho_Real_Par) then
  call Cho_P_IndxSwp()
  call Cho_SetShP2Q(irc,iLoc,iShlAB,nAB)
  call Cho_P_IndxSwp()
else
  call Cho_SetShP2Q(irc,iLoc,iShlAB,nAB)
end if

end subroutine Cho_P_SetShP2Q
