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
      Implicit None
      Integer MaxQual, nSym
#include "cholq.fh"
#include "cho_para_info.fh"

      If (Cho_Real_Par) Then
         l_iQuAB_L = MaxQual*nSym
         l_iQL2G = l_iQuAB_L
         Call GetMem('iQuAB_L','Allo','Inte',ip_iQuAB_L,l_iQuAB_L)
         Call GetMem('iQL2G','Allo','Inte',ip_iQL2G,l_iQL2G)
      Else
         ip_iQuAB_L = -999999
         l_iQuAB_L = 0
         ip_iQL2G = -999999
         l_iQL2G = 0
      End If

      Call Cho_iZero(nQual_L,8)
      Call Cho_iZero(ip_LQ_Sym,8)
      Call Cho_iZero(l_LQ_Sym,8)
      Call Cho_iZero(ldLQ,8)
      ip_LQ = -999999
      l_LQ = 0

      End
