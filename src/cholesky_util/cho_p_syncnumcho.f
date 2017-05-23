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
      SubRoutine Cho_P_SyncNumCho(NumCho,nSym)
C
C     Purpose: sync global NumCho_G vector counter. On entry, NumCho is
C              the local counter (unchanged).
C
      Implicit None
      Integer nSym
      Integer NumCho(nSym)
#include "choglob.fh"
#include "cho_para_info.fh"

      Integer iSym
      Real*8  c1, c2, w1, w2

      If (Cho_Real_Par) Then
         Call Cho_Timer(c1,w1)
         Call iCopy(nSym,NumCho,1,NumCho_G,1)
         Call Cho_GAIGop(NumCho_G,nSym,'max')
         NumChT_G = NumCho_G(1)
         Do iSym = 2,nSym
            NumChT_G = NumChT_G + NumCho_G(iSym)
         End Do
         Call Cho_Timer(c2,w2)
         Call Cho_P_SyncNumCho_Time(c2-c1,w2-w1)
      End If

      End
      SubRoutine Cho_P_SyncNumCho_Time(c,w)
      Implicit None
      Real*8 c, w
#include "cholesky.fh"
      tMisc(1,5)=tMisc(1,5)+c
      tMisc(2,5)=tMisc(2,5)+w
      End
