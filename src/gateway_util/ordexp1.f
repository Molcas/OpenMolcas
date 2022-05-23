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
      SubRoutine OrdExp1(nExp,Exp,nCntrc,Cff)
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
*     Move any exponent which are already decontracted to the bottom
*
      iBottom=nExp
      Do iCntrc = nCntrc, 1, -1
         mExp=0
         jExp=-1
         Do iExp = 1, nExp
            If (Cff(iExp,iCntrc).ne.Zero) Then
               jExp=iExp
               mExp=mExp+1
            End If
         End Do
         If (mExp.eq.1) Then
            Call DSwap_(1,Exp(jExp),1,Exp(iBottom),1)
            Call DSwap_(nCntrc,Cff(jExp,1),nExp,Cff(iBottom,1),nExp)
            iBottom=iBottom-1
         End If
      End Do
*
      Return
      End
