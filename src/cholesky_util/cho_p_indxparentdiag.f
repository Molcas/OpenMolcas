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
      Integer Function Cho_P_IndxParentDiag(iQ,iSym)
C
C     Purpose: return index in global 1st reduced set of qualified iQ,
C              sym. iSym.
C
      Implicit None
      Integer iQ, iSym
#include "cho_para_info.fh"
      Integer  Cho_IndxParentDiag_S
      External Cho_IndxParentDiag_S
      Integer  Cho_IndxParentDiag_P
      External Cho_IndxParentDiag_P

      If (Cho_Real_Par) Then
         Cho_P_IndxParentDiag = Cho_IndxParentDiag_P(iQ,iSym)
      Else
         Cho_P_IndxParentDiag = Cho_IndxParentDiag_S(iQ,iSym)
      End If

      End
      Integer Function Cho_IndxParentDiag_P(iQ,iSym)
      Implicit None
      Integer iQ, iSym
#include "cholesky.fh"
#include "choptr.fh"
#include "WrkSpc.fh"
#include "choglob.fh"

      Integer i, j, iQuAB, IndRed
      IndRed(i,j)=iWork(ip_IndRed_G-1+mmBstRT_G*(j-1)+i)
      iQuAB(i,j)=iWork(ip_iQuAB-1+MaxQual*(j-1)+i)

      Cho_IndxParentDiag_P = IndRed(iQuAB(iQ,iSym),2)

      End
      Integer Function Cho_IndxParentDiag_S(iQ,iSym)
      Implicit None
      Integer iQ, iSym
#include "cholesky.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

      Integer i, j, iQuAB, IndRed
      IndRed(i,j)=iWork(ip_IndRed-1+mmBstRT*(j-1)+i)
      iQuAB(i,j)=iWork(ip_iQuAB-1+MaxQual*(j-1)+i)

      Cho_IndxParentDiag_S = IndRed(iQuAB(iQ,iSym),2)

      End
