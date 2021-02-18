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
      SubRoutine Cho_P_Stat()
      Implicit None
#include "cholesky.fh"
#include "choglob.fh"
#include "cho_para_info.fh"
      Integer iTmp, jTmp

      If (Cho_Real_Par) Then
         Call Cho_P_IndxSwp()
         Call iSwap(nSym,NumCho,1,NumCho_G,1)
         jTmp = NumChT
         NumChT = NumChT_G
         iTmp = LuRed
         LuRed = LuRed_G
         Call Cho_Stat()
         LuRed = iTmp
         NumChT = jTmp
         Call iSwap(nSym,NumCho,1,NumCho_G,1)
         Call Cho_P_IndxSwp()
      Else
         Call Cho_Stat()
      End If

      End
