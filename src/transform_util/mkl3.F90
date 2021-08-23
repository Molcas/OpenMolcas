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
! Copyright (C) 2005, Giovanni Ghigo                                   *
!***********************************************************************

subroutine MkL3(iSymA,iSymI,iI,numV,LyType,iJy,AddLx0,SameLx)
!***********************************************************************
! Author :  Giovanni Ghigo                                             *
!           Lund University, Sweden & Torino University, Italy         *
!           February 2005                                              *
!----------------------------------------------------------------------*
! Purpuse:  Generation of the Cholesky matrix of Secondary(iSymA) for  *
!           occupied iI(iSymI) for numV vectors.                       *
!***********************************************************************

use Cho_Tra, only: nIsh, nSsh, TCVX
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: iSymA, iSymI, iI, numV
integer(kind=iwp), intent(inout) :: LyType, iJy
real(kind=wp), intent(_OUT_) :: AddLx0(*)
logical(kind=iwp), intent(inout) :: SameLx
integer(kind=iwp) :: iAddLx, iAddTCVX, iIx, iV, LxType

! Build Lx
if (iI <= nIsh(iSymI)) then
  LxType = 3
  iIx = iI
else
  LxType = 5
  iIx = iI-nIsh(iSymI)
end if

if (.not. SameLx) then
  LyType = LxType
  iJy = iIx
else
  if ((LyType == LxType) .and. (iIx == iJy)) then
    return
  else
    SameLx = .false.
  end if
end if

iAddLx = 1
iAddTCVX = 1+nSsh(iSymA)*(iIx-1)
do iV=1,numV
  call dCopy_(nSsh(iSymA),TCVX(LxType,iSymA,iSymI)%A(iAddTCVX,iV),1,AddLx0(iAddLx),1)
  iAddLx = iAddLx+nSsh(iSymA)
end do

return

end subroutine MkL3
