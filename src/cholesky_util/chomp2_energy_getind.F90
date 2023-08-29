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
! Copyright (C) 2004,2005, Thomas Bondo Pedersen                       *
!***********************************************************************

subroutine ChoMP2_Energy_GetInd(LnT2am,LiT2am,iBatch,jBatch)
!
! Thomas Bondo Pedersen, Dec. 2004 / Feb. 2005.
!
! Purpose: setup (ai|bj) index arrays for batch i,j.
!          For iBatch=jBatch and ChoAlg=2, (ai|bj) is stored as
!          the matrix M(ab,ij) with i<=j.

use Cholesky, only: nSym
use ChoMP2, only: ChoAlg, LnMatij, LnT1am, nMatab
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: LnT2am, LiT2am(8)
integer(kind=iwp), intent(in) :: iBatch, jBatch
integer(kind=iwp) :: iSym
character(len=14) :: String
character(len=*), parameter :: SecNam = 'ChoMP2_Energy_GetInd'

if (iBatch == jBatch) then
  LnT2am = 0
  if (ChoAlg == 1) then
    do iSym=1,nSym
      LiT2am(iSym) = LnT2am
      LnT2am = LnT2am+LnT1am(iSym,iBatch)*(LnT1am(iSym,iBatch)+1)/2
    end do
  else if (ChoAlg == 2) then
    do iSym=1,nSym
      LiT2am(iSym) = LnT2am
      LnT2am = LnT2am+nMatab(iSym)*LnMatij(iSym,iBatch)
    end do
  else
    write(String,'(A8,I6)') 'ChoAlg =',ChoAlg
    call SysAbendMsg(SecNam,'ChoAlg out-of-bounds error!',String)
  end if
else
  LnT2am = 0
  do iSym=1,nSym
    LiT2am(iSym) = LnT2am
    LnT2am = LnT2am+LnT1am(iSym,iBatch)*LnT1am(iSym,jBatch)
  end do
end if

end subroutine ChoMP2_Energy_GetInd
