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
      SubRoutine Cho_SetGlob()
C
C     Purpose: define entries in choglob.fh
C
      Implicit None
#include "choglob.fh"

      Integer N, iSym

      Integer iLarge
      Parameter (iLarge = 999999)

      ip_Diag_G = -iLarge
      ip_iL2G = -iLarge
      ip_iiBstRSh_G = -iLarge
      ip_nnBstRSh_G = -iLarge
      ip_IndRed_G = -iLarge
      ip_IndRSh_G = -iLarge
      ip_InfRed_G = -iLarge
      ip_InfVec_G = -iLarge
      l_Diag_G = 0
      l_iL2G = 0
      l_iiBstRSh_G = 0
      l_nnBstRSh_G = 0
      l_IndRed_G = 0
      l_IndRSh_G = 0
      l_InfRed_G = 0
      l_InfVec_G = 0

      nnShl_G = 0
      mmBstRT_G = 0

      N = 8*nLoc_G
      Call Cho_iZero(iiBstR_G,N)
      Call Cho_iZero(nnBstR_G,N)
      Call Cho_iZero(nnBstRT_G,nLoc_G)
      Call Cho_iZero(NumCho_G,8)
      NumChT_G = 0

      Do iSym = 1,8
         LuCho_G(iSym) = -iLarge
      End Do
      LuRed_G = -iLarge
      LuRst_G = -iLarge

      End
