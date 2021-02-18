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
      Integer Function Cho_P_LocalSP(iShlAB)
C
C     Purpose: return local shell pair corresponding to global shell
C              pair iShlAB (returns 0 if not found).
C
      use ChoArr, only: MySP, n_MySP
      Implicit None
      Integer iShlAB
#include "cho_para_info.fh"

      Integer iSP

      If (Cho_Real_Par) Then
         Cho_P_LocalSP = 0
         Do iSP = 1,n_mySP
            If (mySP(iSP) .eq. iShlAB) Then
               Cho_P_LocalSP = iSP
               Return
            End If
         End Do
      Else
         Cho_P_LocalSP = iShlAB
      End If

      End
