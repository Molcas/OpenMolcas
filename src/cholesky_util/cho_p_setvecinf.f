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
      SubRoutine Cho_P_SetVecInf(nVec,iSym,iPass)
C
C     Purpose: set global and local info for vectors.
C
      use ChoSwp, only: iQuAB
      Implicit None
      Integer nVec, iSym, iPass
#include "cholesky.fh"
#include "choptr.fh"
#include "WrkSpc.fh"
#include "cho_para_info.fh"
#include "choglob.fh"
      Integer  Cho_P_IndxParentDiag
      External Cho_P_IndxParentDiag

      Integer iV, iVec, iAB

      Integer i, j, IndRed
      IndRed(i,j)=iWork(ip_IndRed-1+mmBstRT*(j-1)+i)

      If (Cho_Real_Par) Then
C Set global vector information (by swapping index arrays)
         Call Cho_P_IndxSwp()
         Do iV = 1,nVec
            iVec = NumCho_G(iSym) + iV
            iAB = IndRed(iQuAB(iV,iSym),2)
            Call Cho_SetVecInf(iWork(ip_InfVec),MaxVec,InfVec_N2,nSym,
     &                         iVec,iSym,iAB,iPass,2)
         End Do
         Call Cho_P_IndxSwp()
C Set local vector information
         Do iV = 1,nVec
            iVec = NumCho_G(iSym) + iV
            iAB = Cho_P_IndxParentDiag(iV,iSym)
            Call Cho_SetVecInf(iWork(ip_InfVec),MaxVec,InfVec_N2,nSym,
     &                         iVec,iSym,iAB,iPass,2)
         End Do
      Else
         Do iV = 1,nVec
            iVec = NumCho(iSym) + iV
            iAB = IndRed(iQuAB(iV,iSym),2)
            Call Cho_SetVecInf(iWork(ip_InfVec),MaxVec,InfVec_N2,nSym,
     &                         iVec,iSym,iAB,iPass,2)
         End Do
      End If

      End
