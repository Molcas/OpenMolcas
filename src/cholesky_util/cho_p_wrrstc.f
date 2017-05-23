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
      SubRoutine Cho_P_WrRstC(iPass)
C
C     Purpose: write global restart info to disk.
C
      Implicit None
      Integer iPass
#include "cholesky.fh"
#include "choglob.fh"
#include "cho_para_info.fh"

      Integer iTmp

      Real*8 c1, c2, w1, w2

      Call Cho_Timer(c1,w1)
      If (Cho_Real_Par) Then
         Call Cho_P_IndxSwp()
         Call iSwap(nSym,NumCho,1,NumCho_G,1)
         iTmp = LuRst
         LuRst = LuRst_G
         Call Cho_WrRstC(iPass)
         LuRst = iTmp
         Call iSwap(nSym,NumCho,1,NumCho_G,1)
         Call Cho_P_IndxSwp()
      Else
         Call Cho_WrRstC(iPass)
      End If
      Call Cho_Timer(c2,w2)
      tMisc(1,3)=tMisc(1,3)+c2-c1
      tMisc(2,3)=tMisc(2,3)+w2-w1

      End
