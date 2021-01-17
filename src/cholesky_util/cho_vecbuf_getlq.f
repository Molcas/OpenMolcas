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
      Implicit Real*8 (a-h,o-z)
      Real*8  QVec(l_QVec)
#include "cholesky.fh"
#include "chovecbuf.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

      Integer nVecTot(8)

      BVec(i,j,k)=Work(ip_ChVBuf_Sym(k)-1+nnBstR(k,2)*(j-1)+i)

C     Check if there is any buffer at all.
C     ------------------------------------

      If (l_ChVBuf .lt. 1) Return

C     Extract from each symmetry block.
C     ---------------------------------

      Call Cho_P_GetGV(nVecTot,nSym)

      kOffQ = 0
      Do iSym = 1,nSym
         If (nQual(iSym) .gt. 0) Then
            Do iVec = 1,nVec_in_Buf(iSym)
               kQ = kOffQ + nQual(iSym)*(iVec-1)
               Do iQ = 1,nQual(iSym)
                  iAB = iQuAB(iQ,iSym) - iiBstR(iSym,2)
                  QVec(kQ+iQ) = BVec(iAB,iVec,iSym)
               End Do
            End Do
            kOffQ = kOffQ + nQual(iSym)*nVecTot(iSym)
         End If
      End Do

      End
