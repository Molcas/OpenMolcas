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
      SubRoutine Cho_P_ZeroDiag(Diag,iSym,iABG)
C
C     Purpose: zero diagonal element iABG (in global diagonal, rs1).
C              For serial runs, this is trivial. For parallel runs, we
C              need first to figure out if the treated diagonal element
C              is in fact present among the qualified in the local
C              diagonal.
C
C     NB! If you wish to test the entire local diagonal (i.e. not just
C         the qualified), use Cho_P_ZeroDiag_Rst instead.
C
      Implicit None
      Real*8  Diag(*)
      Integer iSym, iABG
#include "cho_para_info.fh"
#include "cholesky.fh"
#include "choptr.fh"
#include "choglob.fh"
#include "cholq.fh"
#include "WrkSpc.fh"

      Integer iQ, iAB, jAB, kAB

      Integer i, j
      Integer iQuAB_L, IndRed, iL2G
      iQuAB_L(i,j)=iWork(ip_iQuAB_L-1+MaxQual*(j-1)+i)
      IndRed(i,j)=iWork(ip_IndRed-1+mmBstRT*(j-1)+i)
      iL2G(i)=iWork(ip_iL2G-1+i)

      If (Cho_Real_Par) Then
         Do iQ = 1,nQual_L(iSym)
            iAB = iQuAB_L(iQ,iSym) ! addr in local current rs
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
