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
      SubRoutine Cho_P_ReoQual(iQScr,IDK,nK)
      use ChoSwp, only: iQuAB
      Implicit None
      Integer iQScr(*), IDK(*), nK(*)
#include "cholesky.fh"
#include "choptr.fh"
#include "cho_para_info.fh"

      Call Cho_ReoQual(iQuAB,MaxQual,nSym,iQScr,IDK,nK,nQual)
      If (Cho_Real_Par) Then
         Call Cho_P_QualSwp()
         Call Cho_ReoQual(iQuAB,MaxQual,nSym,iQScr,IDK,nK,nQual)
         Call Cho_P_QualSwp()
      End If

      End
