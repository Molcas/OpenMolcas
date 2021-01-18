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
      SubRoutine Cho_SetQ2(iQuAB2,LstSP,nSP,iSym,jLoc,iLoc)
C
C     Purpose: set mapping from qualified diagonals of symmetry iSym to
C              reduced set indexed by arrays at location iLoc>1,
C              counting only shell pairs that contain qualified
C              diagonals (in the order in which they were qualified).
C              The qualified index array iQuAB (pointer in choptr.fh)
C              is assumed to refer to index arrays at location jLoc>1.
C
      use ChoSwp, only: iQuAB, nnBstRSh, iiBstRSh, IndRSh
      Implicit Real*8 (a-h,o-z)
      Integer iQuAB2(*)
      Integer LstSP(nSP)
#include "cholesky.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

      Integer  Cho_P_LocalSP, Cho_F2SP
      External Cho_P_LocalSP, Cho_F2SP

      IndRed(i,j)=iWork(ip_IndRed-1+mmBstRT*(j-1)+i)

      iC = 0
      jShlAB_Ref = -1
      Do iQ = 1,nQual(iSym)
         jAB = iQuAB(iQ,iSym) ! addr in rs jLoc
         jAB1 = IndRed(jAB,jLoc) ! addr in 1st rs
         jShlAB = Cho_P_LocalSP(Cho_F2SP(IndRSh(jAB1))) ! local SP
         If (jShlAB .ne. jShlAB_Ref) Then
            iC = 0
            iSP = 0
            Do While (iSP .lt. nSP)
               iSP = iSP + 1
               kShlAB = Cho_P_LocalSP(LstSP(iSP)) ! local SP
               If (kShlAB .eq. jShlAB) Then
                  iSP = nSP ! break while loop
               Else
                  iC = iC + nnBstRSh(iSym,kShlAB,iLoc)
               End If
            End Do
            jShlAB_Ref = jShlAB
         End If
         iAB0 = iiBstR(iSym,iLoc) + iiBstRSh(iSym,jShlAB,iLoc)
         iAB = 0
         Do While (iAB .lt. nnBstRSh(iSym,jShlAB,iLoc))
            iAB = iAB + 1
            iAB1 = IndRed(iAB0+iAB,iLoc) ! addr in 1st rs
            If (iAB1 .eq. jAB1) Then
               iQuAB2(iQ) = iC + iAB
               iAB = nnBstRSh(iSym,jShlAB,iLoc) ! break while loop
            End If
         End Do
      End Do

      End
