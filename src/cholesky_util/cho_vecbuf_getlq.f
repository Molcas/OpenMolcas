************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SubRoutine Cho_VecBuf_GetLQ(QVec,l_QVec)
C
C     Purpose: extract elements corresponding to qualified diagonals
C              from vectors in buffer.
C
      use ChoSwp, only: iQuAB
      use ChoVecBuf
      Implicit Real*8 (a-h,o-z)
      Real*8, Target::  QVec(l_QVec)
#include "cholesky.fh"
      Real*8, Pointer:: BVec(:,:), Q(:,:)

      Integer nVecTot(8)
      Integer iS, iE, lRow, lCol

C     Check if there is any buffer at all.
C     ------------------------------------

      If (.NOT.Allocated(CHVBUF)) Return

C     Extract from each symmetry block.
C     ---------------------------------

      Call Cho_P_GetGV(nVecTot,nSym)

      kOffQ = 0
      Do iSym = 1,nSym

         lRow = nnBstR(iSym,2)
         lCol = nVec_in_Buf(iSym)
         iS = ip_ChVBuf_Sym(iSym)
         iE = iS - 1 + lRow*lCol
         BVec(1:lRow,1:lCol) => CHVBUF(iS:iE)

         lRow = nQual(iSym)
         lCol = nVec_in_Buf(iSym)
         iS = kOffQ + 1
         iE = iS - 1 + lRow*lCol
         Q(1:lRow,1:lCol) => QVec(iS:iE)

         If (nQual(iSym) .gt. 0) Then
            Do iVec = 1,nVec_in_Buf(iSym)
*              kQ = kOffQ + nQual(iSym)*(iVec-1)
               Do iQ = 1,nQual(iSym)
                  iAB = iQuAB(iQ,iSym) - iiBstR(iSym,2)
*                 QVec(kQ+iQ) = BVec(iAB,iVec)
                  Q(iQ,iVec) = BVec(iAB,iVec)
               End Do
            End Do
            kOffQ = kOffQ + nQual(iSym)*nVecTot(iSym)
         End If
      End Do

      BVec => Null()

      End
