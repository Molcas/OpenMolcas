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
      SubRoutine Cho_P_GetMaxShl(DiaSh,Smax,iShlAB)
C
C     Purpose: get max shell pair and update DiaSh.
C
      Implicit None
      Real*8  DiaSh(*)
      Real*8  Smax
      Integer iShlAB
#include "cho_para_info.fh"

      If (Cho_Real_Par) Then
C        Swap local and global reduced set index arrays and use original
C        serial routine to get max shell pair.
         Call Cho_P_IndxSwp()
         Call Cho_GetMaxShl(DiaSh,Smax,iShlAB)
         Call Cho_P_IndxSwp()
      Else
         Call Cho_GetMaxShl(DiaSh,Smax,iShlAB)
      End If

      End
