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
      SubRoutine Cho_P_SyncDiag(Diag,iLoc)
C
C     Purpose: synchronize global diagonal. Diag is the local diagonal
C              and iLoc tells which memory location to use for reduced
C              set index arrays.
C
      Implicit None
      Real*8  Diag(*)
      Integer iLoc
#include "cho_para_info.fh"
#include "cholesky.fh"
#include "choglob.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

      Integer kDG, i, j
      Real*8  c1, c2, w1, w2

      Integer iL2G, IndRed
      iL2G(i)=iWork(ip_iL2G-1+i)
      IndRed(i,j)=iWork(ip_IndRed-1+mmBstRT*(j-1)+i)

C     Skip if serial run.
C     -------------------

      If (.not.Cho_Real_Par) Return
      Call Cho_Timer(c1,w1)

C     Zero all entries in global diagonal.
C     ------------------------------------

      Call Cho_dZero(Work(ip_Diag_G),nnBstRT_G(1))

C     Copy elements from local to global diagonal.
C     --------------------------------------------

      kDG = ip_Diag_G - 1
      If (iLoc .eq. 1) Then
         Do i = 1,nnBstRT(1)
            Work(kDG+iL2G(i)) = Diag(i)
         End Do
      Else
         Do j = 1,nnBstRT(iLoc)
            i = IndRed(j,iLoc)
            Work(kDG+iL2G(i)) = Diag(i)
         End Do
      End If

C     Synchronize global diagonal.
C     ----------------------------

      Call Cho_GADGop(Work(ip_Diag_G),nnBstRT_G(1),'+')

      Call Cho_Timer(c2,w2)
      tMisc(1,4)=tMisc(1,4)+c2-c1
      tMisc(2,4)=tMisc(2,4)+w2-w1

      End
