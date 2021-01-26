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
      SubRoutine Cho_P_IniLQ(MaxQual,nSym)
      use ChoSwp, only: iQuAB_L, iQuAB_L_Hidden
      use ChoArr, only: iQL2G
      Implicit None
      Integer MaxQual, nSym
#include "cholq.fh"
#include "cho_para_info.fh"
#include "stdalloc.fh"

      If (Cho_Real_Par) Then
         Call mma_allocate(iQuAB_L_Hidden,MaxQual,nSym,
     &                     Label='iQuAB_L_Hidden')
         iQuAB_L => iQuAB_L_Hidden
         Call mma_allocate(iQL2G,MaxQual,nSym,Label='iQL2G')
      End If

      Call Cho_iZero(nQual_L,8)

      End
