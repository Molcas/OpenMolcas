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

subroutine MkLij(iSymI,iSymJ,iI,iJ,numV,Lij)
!***********************************************************************
! Author :  Giovanni Ghigo                                             *
!           Torino University, Italy                                   *
!           July 2005                                                  *
!----------------------------------------------------------------------*
! Purpuse:  Generation of the Cholesky vector for occupied iI(iSymI),  *
!           iJ(iSymJ) for numV vectors.                                *
!***********************************************************************

use Cho_Tra

implicit real*8(a-h,o-z)
implicit integer(i-n)
integer iSymI, iSymJ, iI, iJ, numV
real*8 Lij(NumV)
#include "rasdim.fh"
#include "SysDef.fh"

if (iI <= nIsh(iSymI)) then
  iIx = iI
  nIx = nIsh(iSymI)
  if (iJ <= nIsh(iSymJ)) then
    LijType = 1
    iJy = iJ
    nJy = nIsh(iSymJ)
  else
    LijType = 7
    iJy = iJ-nIsh(iSymJ)
    nJy = nAsh(iSymJ)
  end if
else
  iIx = iI-nIsh(iSymI)
  nIx = nAsh(iSymI)
  if (iJ <= nIsh(iSymJ)) then
    LijType = 2
    iJy = iJ
    nJy = nIsh(iSymJ)
  else
    LijType = 4
    iJy = iJ-nIsh(iSymJ)
    nJy = nAsh(iSymJ)
  end if
end if

!-----------------------------------------------------------------------
if (IfTest) then
  write(6,*) '     Cho_MkLij: TCVx(',LijType,': ',iSymI,',',iSymJ,')'
  call XFlush(6)
end if
!-----------------------------------------------------------------------
iAddTCVX = iIx+nIx*(iJy-1)
call dCopy_(numV,TCVX(LijType,iSymI,iSymJ)%A(iAddTCVX,1),nIx*nJy,Lij,1)

return

end subroutine MkLij
