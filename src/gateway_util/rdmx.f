************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1991, Roland Lindh                                     *
************************************************************************
      SubRoutine RdMx(RadMax,Exp,nExp,Cff,nCff,cdMax,EtMax)
************************************************************************
*                                                                      *
* Object:                                                              *
*                                                                      *
* Called from: Input                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             August '91                                               *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
      Real*8 Exp(nExp), Cff(nExp,nCff)
#include "real.fh"
*define _DEBUG_
#ifdef _DEBUG_
#include "print.fh"
*
      iRout = 201
      iPrint = nPrint(iRout)
      Call qEnter('RdMx')
*
      Call RecPrt('Exp',' ',Exp,nExp,1)
      Call RecPrt('Cff',' ',Cff,nExp,nCff)
#endif
      Do iExp = 1, nExp
*
         cc = DDot_(nCff,Cff(iExp,1),nExp,Cff(iExp,1),nExp)
         c = Sqrt(cc)
*
         Alpha=Exp(iExp)
         Beta =Exp(iExp)
         Zeta=Alpha+Beta
         If (Zeta.gt.Zero) Then
            Eta=Alpha+Beta
            Rho=(Zeta*Eta)/(Zeta+Eta)
*
            ssss = c**4 * Two*Sqrt(Rho/Pi)*(Pi/Zeta)**(Three/Two)
     &                                    *(Pi/Eta )**(Three/Two)
            If (Sqrt(ssss).gt.RadMax) Then
               RadMax = Sqrt(ssss)
               EtMax  = Eta
               cdMax  = Sqrt(ssss)
            End If
         End If
*
      End Do
*
#ifdef _DEBUG_
      Call qExit('RdMx')
#endif
      Return
      End
