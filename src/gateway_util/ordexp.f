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
      SubRoutine OrdExp(nExp,Exp,nCntrc,Cff)
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
      Real*8 Exp(nExp), Cff(nExp,nCntrc)
*
*     Order exponents
*     Make the subsequent change in the contraction
*     matrix
*
      Do iExp = 1, nExp-1
         Exp1 = Exp(iExp)
         kExp = iExp
         Do jExp = iExp+1, nExp
            Exp2 = Exp(jExp)
            If (Exp2.gt.Exp1) Then
               Exp1 = Exp2
               kExp = jExp
            End If
         End Do
         If (kExp.ne.iExp) Then
            Call DSwap_(1,Exp(iExp),1,Exp(kExp),1)
            Call DSwap_(nCntrc,Cff(iExp,1),nExp,Cff(kExp,1),nExp)
         End If
      End Do
*
      Return
      End
