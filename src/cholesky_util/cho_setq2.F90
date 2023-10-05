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

subroutine Cho_SetQ2(iQuAB2,LstSP,nSP,iSym,jLoc,iLoc)
!
! Purpose: set mapping from qualified diagonals of symmetry iSym to
!          reduced set indexed by arrays at location iLoc>1,
!          counting only shell pairs that contain qualified
!          diagonals (in the order in which they were qualified).
!          The qualified index array iQuAB (pointer in Cholesky)
!          is assumed to refer to index arrays at location jLoc>1.

use Cholesky, only: iiBstR, iiBstRSh, IndRed, IndRSh, iQuAB, nnBstRSh, nQual
use Definitions, only: iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(_OUT_) :: iQuAB2(*)
integer(kind=iwp), intent(in) :: nSP, LstSP(nSP), iSym, jLoc, iLoc
integer(kind=iwp) :: iAB, iAB0, iAB1, iC, iQ, iSP, jAB, jAB1, jShlAB, jShlAB_Ref, kShlAB
integer(kind=iwp), external :: Cho_F2SP, Cho_P_LocalSP

iC = 0
jShlAB_Ref = -1
do iQ=1,nQual(iSym)
  jAB = iQuAB(iQ,iSym) ! addr in rs jLoc
  jAB1 = IndRed(jAB,jLoc) ! addr in 1st rs
  jShlAB = Cho_P_LocalSP(Cho_F2SP(IndRSh(jAB1))) ! local SP
  if (jShlAB /= jShlAB_Ref) then
    iC = 0
    iSP = 0
    do while (iSP < nSP)
      iSP = iSP+1
      kShlAB = Cho_P_LocalSP(LstSP(iSP)) ! local SP
      if (kShlAB == jShlAB) then
        iSP = nSP ! break while loop
      else
        iC = iC+nnBstRSh(iSym,kShlAB,iLoc)
      end if
    end do
    jShlAB_Ref = jShlAB
  end if
  iAB0 = iiBstR(iSym,iLoc)+iiBstRSh(iSym,jShlAB,iLoc)
  iAB = 0
  do while (iAB < nnBstRSh(iSym,jShlAB,iLoc))
    iAB = iAB+1
    iAB1 = IndRed(iAB0+iAB,iLoc) ! addr in 1st rs
    if (iAB1 == jAB1) then
      iQuAB2(iQ) = iC+iAB
      iAB = nnBstRSh(iSym,jShlAB,iLoc) ! break while loop
    end if
  end do
end do

end subroutine Cho_SetQ2
