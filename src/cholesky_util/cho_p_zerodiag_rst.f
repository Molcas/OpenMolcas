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
      SubRoutine Cho_P_ZeroDiag_Rst(Diag,iSym,iABG)
C
C     Purpose: zero diagonal element iABG (in global diagonal, rs1).
C              For serial runs, this is trivial. For parallel runs, we
C              need first to figure out if the treated diagonal element
C              is in fact present in the local diagonal.
C
      Implicit None
      Real*8  Diag(*)
      Integer iSym, iABG
#include "cho_para_info.fh"
#include "cholesky.fh"
#include "choptr.fh"
#include "choglob.fh"
#include "WrkSpc.fh"

      Integer iAB1, iAB2, iAB, jAB, kAB

      Integer i, j
      Integer IndRed, iL2G
      IndRed(i,j)=iWork(ip_IndRed-1+mmBstRT*(j-1)+i)
      iL2G(i)=iWork(ip_iL2G-1+i)

      If (Cho_Real_Par) Then
         iAB1 = iiBstR(iSym,2) + 1
         iAB2 = iAB1 + nnBstR(iSym,2) - 1
         Do iAB = iAB1,iAB2
            jAB = IndRed(iAB,2)  ! addr in local rs1
            kAB = iL2G(jAB)      ! addr in global rs1
            If (kAB .eq. iABG) Then ! found...
               Diag(jAB) = 0.0d0    ! now zero local diagonal elm.
               Return
            End If
         End Do
      Else
         Diag(iABG) = 0.0d0
      End If

      End
