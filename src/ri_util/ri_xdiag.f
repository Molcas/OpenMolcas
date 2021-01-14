************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2007, Thomas Bondo Pedersen                            *
************************************************************************
      SubRoutine RI_XDiag(Diag,nDiag)
C
C     Thomas Bondo Pedersen, Jan. 2007.
C
C     Purpose: compute exact integral diagonal.
C
      use ChoArr, only: iSP2F
      Implicit None
      Integer nDiag
      Real*8  Diag(nDiag)
#include "cholesky.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

      Integer ip0, ip_Scr, l_Scr
      Integer ip_SewMem, l_SewMem
      Integer ID
      Integer iSAB, iShlA, iShlB
      Integer NumAB
      Integer iSym, i1, i2

      Logical  Rsv_Tsk
      External Rsv_Tsk

      Integer i, j, k
      Integer nBstSh, IndRed, iiBstRSh, nnBstRSh
      nBstSh(i)=iWork(ip_nBstSh-1+i)
      IndRed(i,j)=iWork(ip_IndRed-1+mmBstRT*(j-1)+i)
      iiBstRSh(i,j,k)=iWork(ip_iiBstRSh-1+nSym*nnShl*(k-1)+nSym*(j-1)+i)
      nnBstRSh(i,j,k)=iWork(ip_nnBstRSh-1+nSym*nnShl*(k-1)+nSym*(j-1)+i)

C     Allocate memory.
C     ----------------

      Call Init_Tsk(ID,nnShl)

      l_Scr = Mx2Sh
      Call GetMem('xDiaScr','Allo','Real',ip_Scr,l_Scr)
      Call GetMem('getMx','Max ','Real',ip_SewMem,l_SewMem)
      Call xSetMem_Ints(l_SewMem)

C     Initialize diagonal array.
C     --------------------------

      Call fZero(Diag,nnBstRT(1))

C     Parallel loop over shell pairs in first red. set.
C     -------------------------------------------------

      ip0 = ip_Scr - 1
      Do While (Rsv_Tsk(ID,iSAB))

C        Get shells.
C        -----------

         Call Cho_InvPck(iSP2F(iSAB),iShlA,iShlB,.True.)

C        Compute (AB|AB).
C        ----------------

         If (iShlA .eq. iShlB) Then
            NumAB = nBstSh(iShlA)*(nBstSh(iShlA)+1)/2
         Else
            NumAB = nBstSh(iShlA)*nBstSh(iShlB)
         End If
         ShA = iShlA
         ShB = iShlB
         Call Cho_MCA_DiagInt(iShlA,iShlB,Work(ip_Scr),NumAB)

C        Extract diagonal elements.
C        --------------------------

         Do iSym = 1,nSym
            i1 = iiBstR(iSym,1) + iiBstRSh(iSym,iSAB,1) + 1
            i2 = i1 + nnBstRSh(iSym,iSAB,1) - 1
            Do i = i1,i2
               Diag(i) = Work(ip0+IndRed(i,1))
            End Do
         End Do

      End Do
      Call Free_Tsk(ID)

C     Sync diagonal.
C     --------------

      Call GAdGOP(Diag,nnBstRT(1),'+')

C     Deallocate memory.
C     ------------------

      Call xRlsMem_Ints
      Call GetMem('xDiaScr','Free','Real',ip_Scr,l_Scr)

      End
