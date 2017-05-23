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
      SubRoutine Cho_P_WrDiag()
C
C     Purpose: store global diagonal on disk (Parallel only).
C              NB: on exit, initial global diagonal is stored!
C
      Implicit None
#include "WrkSpc.fh"
#include "cholesky.fh"
#include "choglob.fh"
#include "cho_para_info.fh"
      Integer ip_Diag_L
      Integer l_Diag_L

      If (Cho_Real_Par) Then
         l_Diag_L = nnBstRT(1)
         Call GetMem('DiagL','Allo','Real',ip_Diag_L,l_Diag_L)
         Call Cho_IODiag(Work(ip_Diag_L),2)
         Call Cho_P_SyncDiag(Work(ip_Diag_L),1)
         Call Cho_P_IndxSwp()
         Call Cho_IODiag(Work(ip_Diag_G),1)
         Call Cho_P_IndxSwp()
         Call GetMem('DiagL','Free','Real',ip_Diag_L,l_Diag_L)
      End If

      End
