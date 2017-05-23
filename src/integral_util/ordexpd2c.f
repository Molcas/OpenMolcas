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
      SubRoutine OrdExpD2C(nExp,Exp,nCntrc,Cff)
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
      Real*8 Exp(nExp), Cff(nExp,nCntrc)
*
*     Order exponents diffuse to compact
*     Make the subsequent change in the contraction
*     matrix
*
      Do iExp = 1, nExp-1
         Exp1 = Exp(iExp)
         kExp = iExp
         Do jExp = iExp+1, nExp
            Exp2 = Exp(jExp)
            If (Exp2.lt.Exp1) Then
               Exp1 = Exp2
               kExp = jExp
            End If
         End Do
         If (kExp.ne.iExp) Then
            Call DSwap_(1,Exp(iExp),1,Exp(kExp),1)
            Call DSwap_(nCntrc,Cff(iExp,1),nExp,Cff(kExp,1),nExp)
         End If
      End Do
#ifdef _ORDER_BAS_
*
*     Now order the contracted basis functions diffuse to compact
*
      Do iCntrc = 1, nCntrc-1
         Bas1=Abs(Cff(1,iCntrc))
         kCntrc = iCntrc
         Do jCntrc = iCntrc+1, nCntrc
            Bas2=Abs(Cff(1,jCntrc))
            If (Bas2.lt.Bas1) Then
               Bas1 = Bas2
               kCntrc = jCntrc
            End If
         End Do
         If (kCntrc.ne.iCntrc) Then
            Call DSwap_(nExp,Cff(1,iCntrc),1,Cff(1,kCntrc),1)
         End If
      End Do
#endif
*
      Return
      End
