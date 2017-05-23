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
      SubRoutine Cho_P_ChkInt(xInt,Diag,iSym,nErr,Tol,Report)
C
C     Purpose: check diagonals in qualified integral columns.
C
      Implicit None
      Real*8  xInt(*), Diag(*)
      Integer iSym, nErr
      Real*8  Tol
      Logical Report
#include "cho_para_info.fh"

      If (Cho_Real_Par) Then
         Call Cho_P_QualSwp()
         Call Cho_ChkInt(xInt,Diag,iSym,nErr,Tol,Report)
         Call Cho_P_QualSwp()
      Else
         Call Cho_ChkInt(xInt,Diag,iSym,nErr,Tol,Report)
      End If

      End
